import std/[strutils, strformat, tables, math, algorithm]
import lapper
import hts/files
import cligen

type
  RegionIv = object
    ## An interval in the lapper tree. `idx` points into the `seq[Region]`
    ## so the hot loop avoids string hashing and string copies.
    startPos, stopPos: int
    idx: int

  DepthBlock = tuple[depth: float, width: int]

  Region = object
    ## Per-region accumulator, filled while streaming the per-base BED.
    id: string
    length: int
    totalDepth: float
    logSum: float            # sum of width * depth * ln(depth), used for evenness
    blocks: seq[DepthBlock]

  RegionStats = object
    length: int
    fractionCovered, totalDepth, minDepth, maxDepth: float
    meanDepth, medianDepth, cv, evenness: float
    breadths: seq[float]

  Trees = Table[string, Lapper[RegionIv]]

proc start(iv: RegionIv): int {.inline.} = iv.startPos
proc stop(iv: RegionIv): int {.inline.} = iv.stopPos

proc openOrQuit(path, what: string): HTSFile =
  if not result.open(path, "r"):
    quit(&"Error: Could not open {what} {path}")

proc mergeIntervals(ivs: var seq[(string, int, int)]): seq[(string, int, int)] =
  ## Sort and merge overlapping/abutting intervals so bases are never counted twice.
  ivs.sort()
  for (chrom, s, e) in ivs:
    if result.len > 0 and result[^1][0] == chrom and s <= result[^1][2]:
      result[^1][2] = max(result[^1][2], e)
    else:
      result.add((chrom, s, e))

proc buildRegionIndex(regionsBed: string): (Trees, seq[Region]) =
  ## Read the regions BED. Records sharing a name are aggregated into one region,
  ## and output order follows first appearance in the file.
  var
    f = openOrQuit(regionsBed, "regions BED file")
    regions: seq[Region]
    idxById: Table[string, int]
    rawIvs: seq[seq[(string, int, int)]]   # per region: (chrom, start, stop)
    line = newStringOfCap(2048)
  defer: f.close()

  while f.readLine(line):
    let line = line.strip()
    if line.len == 0 or line.startsWith("#") or line.startsWith("track") or
       line.startsWith("browser"):
      continue
    let cols = line.split('\t')
    if cols.len < 3: continue

    let
      chrom = cols[0]
      start = parseInt(cols[1])
      stop = parseInt(cols[2])
    if stop <= start:
      quit(&"Error: invalid interval {chrom}:{start}-{stop} in {regionsBed}")

    let id = if cols.len > 3 and cols[3].strip().len > 0: cols[3].strip()
             else: &"{chrom}:{start}-{stop}"

    let idx = idxById.mgetOrPut(id, regions.len)
    if idx == regions.len:
      regions.add(Region(id: id))
      rawIvs.add(@[])
    rawIvs[idx].add((chrom, start, stop))

  var ivsByChrom: Table[string, seq[RegionIv]]
  for idx, ivs in rawIvs.mpairs:
    for (chrom, s, e) in mergeIntervals(ivs):
      regions[idx].length += e - s
      ivsByChrom.mgetOrPut(chrom, @[]).add(RegionIv(startPos: s, stopPos: e, idx: idx))

  var trees: Trees
  for chrom, ivs in ivsByChrom.mpairs:
    trees[chrom] = lapify(ivs)
  (trees, regions)

proc accumulateDepth(bedFile: string, trees: var Trees, regions: var seq[Region]) =
  ## Stream the per-base BED and add each block's overlap to its region(s).
  var
    f = openOrQuit(bedFile, "coverage BED file")
    line = newStringOfCap(256)
  defer: f.close()

  while f.readLine(line):
    let cols = line.split('\t', maxsplit = 4)
    if cols.len < 4: continue

    let depth = parseFloat(cols[3].strip())
    if depth == 0.0: continue
    trees.withValue(cols[0], tree):
      let
        bStart = parseInt(cols[1])
        bStop = parseInt(cols[2])
        logDepth = ln(depth)
      for iv in tree[].find(bStart, bStop):
        let width = min(bStop, iv.stopPos) - max(bStart, iv.startPos)
        if width > 0:
          let w = float(width)
          regions[iv.idx].totalDepth += w * depth
          regions[iv.idx].logSum += w * depth * logDepth
          regions[iv.idx].blocks.add((depth, width))

proc computeStats(r: Region, thresholds: seq[int]): RegionStats =
  let L = r.length
  result = RegionStats(length: L, totalDepth: r.totalDepth,
                       breadths: newSeq[float](thresholds.len))
  if L <= 0: return

  let fL = float(L)
  result.meanDepth = r.totalDepth / fL

  var blocks = r.blocks
  var covered = 0
  for b in blocks: covered += b.width
  result.fractionCovered = float(covered) / fL
  if covered < L:
    blocks.add((0.0, L - covered))   # uncovered bases have depth 0

  blocks.sort(proc(x, y: DepthBlock): int = cmp(x.depth, y.depth))
  result.minDepth = blocks[0].depth
  result.maxDepth = blocks[^1].depth

  # Single pass over sorted blocks for median, variance and threshold breadths.
  let (mid1, mid2) = ((L - 1) div 2, L div 2)
  var
    seen = 0
    med1, med2 = NaN
    sumSqDiff = 0.0
    thresholdCounts = newSeq[int](thresholds.len)
  for b in blocks:
    if med1.isNaN and seen + b.width > mid1: med1 = b.depth
    if med2.isNaN and seen + b.width > mid2: med2 = b.depth
    seen += b.width

    sumSqDiff += (b.depth - result.meanDepth) ^ 2 * float(b.width)

    for i, t in thresholds:
      if b.depth >= float(t): thresholdCounts[i] += b.width

  result.medianDepth = (med1 + med2) / 2.0
  if result.meanDepth > 0.0:
    result.cv = sqrt(sumSqDiff / fL) / result.meanDepth
  for i, n in thresholdCounts:
    result.breadths[i] = float(n) / fL

  # Evenness: exp(Shannon entropy of depth distribution) / length
  if r.totalDepth > 0.0:
    let entropy = max(0.0, ln(r.totalDepth) - r.logSum / r.totalDepth)
    result.evenness = exp(entropy) / fL

proc main(bed: string, regions: string, thresholds: seq[int] = @[1, 10, 100, 1000],
          output: string = "") =
  var (trees, regs) = buildRegionIndex(regions)
  accumulateDepth(bed, trees, regs)

  let outStream = if output.len > 0: open(output, fmWrite) else: stdout
  defer:
    if outStream != stdout: outStream.close()

  var header = @["region_id", "length", "fraction_covered", "total_depth", "min_depth",
                 "max_depth", "mean_depth", "median_depth", "cv"]
  for t in thresholds: header.add(&"F{t}")
  header.add("evenness")
  outStream.writeLine(header.join("\t"))

  for r in regs:
    let s = computeStats(r, thresholds)
    var row = @[r.id, $s.length, &"{s.fractionCovered:.4f}", &"{s.totalDepth:.2f}",
                &"{s.minDepth:.2f}", &"{s.maxDepth:.2f}", &"{s.meanDepth:.2f}",
                &"{s.medianDepth:.2f}", &"{s.cv:.4f}"]
    for b in s.breadths: row.add(&"{b:.4f}")
    row.add(&"{s.evenness:.2f}")
    outStream.writeLine(row.join("\t"))

when isMainModule:
  dispatch main, help = {
    "bed": "Path to the mosdepth per-base bed.gz file",
    "regions": "Path to the input BED regions file to build interval trees over",
    "thresholds": "Comma-separated list of depth thresholds (e.g. 1,10,100,1000) to compute breadth fractions for",
    "output": "Optional path to save the TSV output. Defaults to stdout."
  }
