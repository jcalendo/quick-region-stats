import std/[strutils, tables, math, algorithm]
import lapper
import hts/files
import cligen

type
  RegionIv* = object
    startPos*, stopPos*: int
    region_id*: string

  RegionStats = object
    evenness, minCov, maxCov, meanCov, medianCov, fractionCovered, totalCov, cv: float
    length: int
    breadths: seq[float]

proc start*(iv: RegionIv): int {.inline.} = iv.startPos
proc stop*(iv: RegionIv): int {.inline.} = iv.stopPos

proc buildRegionIndex(regionsBed: string): tuple[trees: Table[string, Lapper[RegionIv]], regionLengths: Table[string, int]] =
  var 
    regionDict = initTable[string, seq[RegionIv]]()
    regionLengths = initTable[string, int]()
    f: HTSFile

  try:
    discard f.open(regionsBed, "r")
  except Exception:
    quit("Error: Could not open regions BED file " & regionsBed)

  var line = newStringOfCap(2048)

  while f.readLine(line):
    if line.startsWith("#") or line.strip() == "": continue
    let cols = line.split('\t')
    if cols.len < 3: continue

    let chrom = cols[0]
    let start = parseInt(cols[1])
    let stop = parseInt(cols[2])
    
    let regionId = if cols.len > 3 and cols[3].strip() != "": cols[3].strip() else: chrom & ":" & $start & "-" & $stop
    
    if regionId.len > 0:
      let length = stop - start
      regionLengths[regionId] = regionLengths.getOrDefault(regionId, 0) + length
      
      if not regionDict.hasKey(chrom):
        regionDict[chrom] = @[]
      regionDict[chrom].add(RegionIv(startPos: start, stopPos: stop, region_id: regionId))
  
  f.close()

  var trees = initTable[string, Lapper[RegionIv]]()
  for chrom, regions in regionDict.mpairs:
    trees[chrom] = lapify(regions)
    
  return (trees, regionLengths)


proc parseBedCoverage(bedFile: string, trees: var Table[string, Lapper[RegionIv]]): tuple[totalCov: Table[string, float], logSum: Table[string, float], covBlocks: Table[string, seq[tuple[c: float, w: int]]]] =
  var 
    regionTotalCov = initTable[string, float]()
    regionLogSum = initTable[string, float]()
    regionCovBlocks = initTable[string, seq[tuple[c: float, w: int]]]()
    f: HTSFile

  try:
    discard f.open(bedFile, "r")
  except Exception:
    quit("Error: Could not open coverage BED file " & bedFile)

  var line = newStringOfCap(2048)
  while f.readLine(line):
    let cols = line.split('\t')
    if cols.len < 4: continue
    
    let chrom = cols[0]
    let bStart = parseInt(cols[1])
    let bStop = parseInt(cols[2])
    let cov = parseFloat(cols[3])

    if cov == 0.0 or not trees.hasKey(chrom): 
      continue
      
    for iv in trees[chrom].find(bStart, bStop):
      let rId = iv.region_id
      let overlapStart = max(bStart, iv.startPos)
      let overlapStop = min(bStop, iv.stopPos)
      let width = float(overlapStop - overlapStart)

      if width > 0.0:
        regionTotalCov[rId] = regionTotalCov.getOrDefault(rId, 0.0) + (width * cov)
        regionLogSum[rId] = regionLogSum.getOrDefault(rId, 0.0) + (width * cov * ln(cov))
        
        if not regionCovBlocks.hasKey(rId): regionCovBlocks[rId] = @[]
        regionCovBlocks[rId].add((cov, int(width)))

  f.close()
  return (regionTotalCov, regionLogSum, regionCovBlocks)

proc computeRegionStats(regionLengths: Table[string, int], regionTotalCov: Table[string, float], regionLogSum: Table[string, float], regionCovBlocks: Table[string, seq[tuple[c: float, w: int]]], thresholds: seq[int]): Table[string, RegionStats] =
  var results = initTable[string, RegionStats]()
  
  for rId, L in regionLengths:
    let C = regionTotalCov.getOrDefault(rId, 0.0)
    let logSum = regionLogSum.getOrDefault(rId, 0.0)
    let blocks = regionCovBlocks.getOrDefault(rId, @[])
    
    var stat = RegionStats(length: L, evenness: 0.0, minCov: 0.0, maxCov: 0.0, meanCov: 0.0, medianCov: 0.0, fractionCovered: 0.0, totalCov: C, cv: 0.0)
    stat.breadths = newSeq[float](thresholds.len)
    
    if L > 0:
      stat.meanCov = C / float(L)
      
      var fullBlocks = blocks
      var covered = 0
      for b in blocks: covered += b.w
      stat.fractionCovered = float(covered) / float(L)
      
      # Add 0-coverage block if region isn't fully covered
      if covered < L:
        fullBlocks.add((0.0, L - covered))
        
      # Sort blocks for median/min/max
      fullBlocks.sort(proc(x, y: tuple[c: float, w: int]): int = cmp(x.c, y.c))
      
      if fullBlocks.len > 0:
        stat.minCov = fullBlocks[0].c
        stat.maxCov = fullBlocks[^1].c
        
        let mid1 = (L - 1) div 2
        let mid2 = L div 2
        var current = 0
        var med1, med2 = -1.0
        
        # Calculate Median, CV, and Threshold Breadths simultaneously over fullBlocks
        var sumSqDiff = 0.0
        var thresholdCounts = newSeq[int](thresholds.len)

        for b in fullBlocks:
          # Median logic
          if med1 < 0.0 and current + b.w > mid1: med1 = b.c
          if med2 < 0.0 and current + b.w > mid2: med2 = b.c
          current += b.w
          
          # Variance logic
          let diff = b.c - stat.meanCov
          sumSqDiff += (diff * diff) * float(b.w)

          # Breadth thresholds logic
          for i, t in thresholds:
            if b.c >= float(t):
              thresholdCounts[i] += b.w
          
        stat.medianCov = (med1 + med2) / 2.0
        
        # Finalize CV
        let variance = sumSqDiff / float(L)
        if stat.meanCov > 0.0:
          stat.cv = sqrt(variance) / stat.meanCov
          
        # Finalize Breadths
        for i in 0 ..< thresholds.len:
          stat.breadths[i] = float(thresholdCounts[i]) / float(L)

      # Calculate Evenness
      if C > 0.0:
        var entropy = ln(C) - (logSum / C)
        if entropy < 0.0: entropy = 0.0 
        stat.evenness = exp(entropy) / float(L)

    results[rId] = stat
    
  return results

proc main(bed: string, regions: string, thresholds: seq[int] = @[1, 10, 100, 1000], output: string = "") =
  var (trees, regionLengths) = buildRegionIndex(regions)
  let (regionTotalCov, regionLogSum, regionCovBlocks) = parseBedCoverage(bed, trees)
  let stats = computeRegionStats(regionLengths, regionTotalCov, regionLogSum, regionCovBlocks, thresholds)

  var outStream: File
  if output.len > 0:
    if not open(outStream, output, fmWrite):
      quit("Error: Could not open output file " & output & " for writing.")
  else:
    outStream = stdout 

  # Build dynamic header
  var header = "region_id\tlength\tfraction_covered\ttotal_cov\tmin_cov\tmax_cov\tmean_cov\tmedian_cov\tcv"
  for t in thresholds:
    header &= "\tF" & $t
  header &= "\tevenness"
  
  outStream.writeLine(header)
  
  # Write results
  for rId, s in stats:
    outStream.write(rId, "\t", 
                    s.length, "\t", 
                    formatFloat(s.fractionCovered, ffDecimal, 4), "\t",
                    formatFloat(s.totalCov, ffDecimal, 2), "\t",
                    formatFloat(s.minCov, ffDecimal, 2), "\t",
                    formatFloat(s.maxCov, ffDecimal, 2), "\t",
                    formatFloat(s.meanCov, ffDecimal, 2), "\t",
                    formatFloat(s.medianCov, ffDecimal, 2), "\t",
                    formatFloat(s.cv, ffDecimal, 4))
    
    for b in s.breadths:
      outStream.write("\t", formatFloat(b, ffDecimal, 4))
      
    outStream.writeLine("\t", formatFloat(s.evenness, ffDecimal, 2))

  if output.len > 0:
    outStream.close()

when isMainModule:
  dispatch main, help = {
    "bed": "Path to the mosdepth per-base bed.gz file",
    "regions": "Path to the input BED regions file to build interval trees over",
    "thresholds": "Comma-separated list of depth thresholds (e.g. 1,10,100,1000) to compute breadth fractions for",
    "output": "Optional path to save the TSV output. Defaults to stdout."
  }