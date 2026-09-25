import std/[unittest, math, os, tables]
import qrs

const
  dataDir = currentSourcePath.parentDir
  thresholds = @[1, 10, 100, 1000]

proc runRegions(regionsBed, bed: string): seq[Region] =
  var (trees, regions) = buildRegionIndex(regionsBed)
  accumulateDepth(bed, trees, regions)
  regions

proc statsById(regions: seq[Region]): Table[string, RegionStats] =
  for r in regions: result[r.id] = computeStats(r, thresholds)

proc writeTemp(name, content: string): string =
  result = getTempDir() / name
  writeFile(result, content)

proc checkBreadths(actual, expected: seq[float]) =
  check actual.len == expected.len
  for i in 0 ..< expected.len:
    check actual[i].almostEqual(expected[i])

suite "test regions":
  let regions = runRegions(dataDir / "regions.bed", dataDir / "per-base.bed")
  let stats = statsById(regions)

  test "rows follow regions BED order":
    var ids: seq[string]
    for r in regions: ids.add(r.id)
    check ids == @["reg_uniform", "reg_variable", "reg_partial", "reg_empty"]

  test "reg_uniform: 100bp at 50x":
    let s = stats["reg_uniform"]
    check s.length == 100
    check s.fractionCovered == 1.0
    check s.totalDepth == 5000.0
    check s.minDepth == 50.0
    check s.maxDepth == 50.0
    check s.meanDepth == 50.0
    check s.medianDepth == 50.0
    check s.cv == 0.0
    check s.evenness.almostEqual(1.0)
    checkBreadths(s.breadths, @[1.0, 1.0, 0.0, 0.0])

  test "reg_variable: 50bp at 10x, 50bp at 50x":
    let s = stats["reg_variable"]
    check s.length == 100
    check s.fractionCovered == 1.0
    check s.totalDepth == 3000.0
    check s.minDepth == 10.0
    check s.maxDepth == 50.0
    check s.meanDepth == 30.0
    check s.medianDepth == 30.0   # even length: mean of depths at bases 49 and 50
    check s.cv.almostEqual(20.0 / 30.0)
    let expEven = exp(ln(3000.0) - (500.0 * ln(10.0) + 2500.0 * ln(50.0)) / 3000.0) / 100.0
    check s.evenness.almostEqual(expEven)
    checkBreadths(s.breadths, @[1.0, 1.0, 0.0, 0.0])

  test "reg_partial: 40bp at 5x, 40bp uncovered, 20bp at 100x":
    let s = stats["reg_partial"]
    check s.length == 100
    check s.fractionCovered.almostEqual(0.6)
    check s.totalDepth == 2200.0
    check s.minDepth == 0.0
    check s.maxDepth == 100.0
    check s.meanDepth == 22.0
    check s.medianDepth == 5.0
    let variance = (40.0 * 22.0^2 + 40.0 * 17.0^2 + 20.0 * 78.0^2) / 100.0
    check s.cv.almostEqual(sqrt(variance) / 22.0)
    let expEven = exp(ln(2200.0) - (200.0 * ln(5.0) + 2000.0 * ln(100.0)) / 2200.0) / 100.0
    check s.evenness.almostEqual(expEven)
    checkBreadths(s.breadths, @[0.6, 0.2, 0.2, 0.0])

  test "reg_empty: chromosome absent from per-base BED":
    let s = stats["reg_empty"]
    check s.length == 100
    check s.fractionCovered == 0.0
    check s.totalDepth == 0.0
    check s.minDepth == 0.0
    check s.maxDepth == 0.0
    check s.meanDepth == 0.0
    check s.medianDepth == 0.0
    check s.cv == 0.0
    check s.evenness == 0.0
    checkBreadths(s.breadths, @[0.0, 0.0, 0.0, 0.0])

suite "region parsing":
  test "overlapping intervals with the same name stay separate":
    let regionsBed = writeTemp("qrs_overlap_regions.bed",
      "chr1\t1000\t1100\tgeneA\nchr1\t1050\t1150\tgeneA\n")
    let bed = writeTemp("qrs_overlap_cov.bed", "chr1\t1000\t1120\t10\n")
    let regions = runRegions(regionsBed, bed)
    check regions.len == 2
    check (regions[0].startPos, regions[0].stopPos) == (1000, 1100)
    check (regions[1].startPos, regions[1].stopPos) == (1050, 1150)
    let (s0, s1) = (computeStats(regions[0], thresholds), computeStats(regions[1], thresholds))
    check s0.length == 100
    check s0.totalDepth == 1000.0
    check s0.fractionCovered == 1.0
    check s1.length == 100
    check s1.totalDepth == 700.0
    check s1.fractionCovered.almostEqual(0.7)

  test "same-name intervals on separate chromosomes stay separate":
    let regionsBed = writeTemp("qrs_multi_regions.bed",
      "chr1\t0\t10\tgeneB\nchr2\t0\t30\tgeneB\n")
    let bed = writeTemp("qrs_multi_cov.bed", "chr1\t0\t10\t4\nchr2\t0\t10\t8\n")
    let regions = runRegions(regionsBed, bed)
    check regions.len == 2
    check regions[0].chrom == "chr1"
    check regions[1].chrom == "chr2"
    let (s0, s1) = (computeStats(regions[0], thresholds), computeStats(regions[1], thresholds))
    check s0.length == 10
    check s0.totalDepth == 40.0
    check s0.fractionCovered == 1.0
    check s1.length == 30
    check s1.totalDepth == 80.0
    check s1.fractionCovered.almostEqual(1.0 / 3.0)

  test "records with invalid coordinates are skipped":
    let regionsBed = writeTemp("qrs_invalid_regions.bed",
      "chr1\t0\t10\tok1\n" &
      "chr1\t100\t50\tbackwards\n" &
      "chr1\t20\t20\tempty\n" &
      "chr1\t-5\t10\tnegative\n" &
      "chr1\tabc\t10\tnotanumber\n" &
      "chr1\t30\t40\tok2\n")
    let (_, regions) = buildRegionIndex(regionsBed)
    var ids: seq[string]
    for r in regions: ids.add(r.id)
    check ids == @["ok1", "ok2"]

  test "missing name column falls back to chrom:start-stop; headers skipped":
    let regionsBed = writeTemp("qrs_noname_regions.bed",
      "track name=x\n# comment\n\nchr1\t1000\t1100\n")
    let (_, regions) = buildRegionIndex(regionsBed)
    check regions.len == 1
    check regions[0].id == "chr1:1000-1100"
    check regions[0].length == 100
