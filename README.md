## Quick region stats from per-base BED files (qrs)

A simple program for computing depth statistics over genomic regions given a per-base BED file
and a BED file specifying regions for which to compute stats over. 

## Usage

```{bash}
Usage:
  main [REQUIRED,optional-params] 
Options:
  -h, --help                                 print this cligen-erated help
  --help-syntax                              advanced: prepend,plurals,..
  -b=, --bed=         string  REQUIRED       Path to the mosdepth per-base bed.gz file
  -r=, --regions=     string  REQUIRED       Path to the input BED regions file to build interval trees over
  -t=, --thresholds=  ints    1,10,100,1000  Comma-separated list of depth thresholds (e.g. 1,10,100,1000) to compute breadth fractions for
  -o=, --output=      string  ""             Optional path to save the TSV output. Defaults to stdout.
```

The program will output one line for each valid record in the regions.bed file, in the same order
as the input. Each record is treated as its own region, even if several records share a name.
Records with invalid coordinates (non-integer, negative start, or end <= start) are skipped with a
warning on stderr. The output has the following columns:

- **chrom**, **start**, **end**: Coordinates of the region, as given in the BED file
- **region_id**: Name of the record in the BED file (if 'Name' column is present), otherwise `chrom:start-end`
- **length**: Length of the region
- **fraction_covered**: Fraction of bases in the region with depth > 0
- **total_depth**: Total depth of the region
- **min_depth**: Minimum depth observed
- **max_depth**: maximum depth observed
- **mean_depth**: Average depth
- **median_depth**: Median depth
- **cv**: Coefficient of variation of the depth across the region
- **F1,10,100,...**: Fraction of the region covered at 1x, 10x, 100x, etc.
- **evenness**: Normalized observed entropy

### Example

```{bash}
qrs -b per-base.bed.gz -r regions.bed -o stats.tsv 
```

## Installation

Check out the [releases page](https://github.com/jcalendo/quick-region-stats/releases/) for a 
pre-compiled linux binary. 

To build the script, a system installation of [htslib](https://www.htslib.org/download/) and of 
course [Nim](https://nim-lang.org/install.html) are required. Once these dependencies are met:

```{bash}
git clone https://github.com/jcalendo/quick-region-stats.git
cd quick-region-stats

# Fetch required Nim packages
nimble install -d 

# Compile with speed optimizations
nim c -d:danger -d:release --opt:speed src/qrs.nim
```
