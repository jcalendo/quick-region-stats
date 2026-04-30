## Quick region stats from per-base BED files (qrs)

A simple program for computing coverage statistics over genomic regions given a per-base BED file
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
  -t=, --thresholds=  ints    1,10,100,1000  Comma-separated list of coverage thresholds (e.g. 1,10,100,1000) to compute breadth fractions for
  -o=, --output=      string  ""             Optional path to save the TSV output. Defaults to stdout.
```

The program will output one line for each region in the regions.bed file with the following columns:

- **region_id**: Name of the record in the BED file (if 'Name' column is present) 
- **length**: Length of the region measured
- **fraction_covered**: Fraction of bases in the region with coverage > 0
- **total_cov**: Total coverage of the region
- **min_cov**: Minimum coverage observed
- **max_cov**: maximum coverage observed
- **mean_cov**: Average coverage
- **median_cov**: Median coverage
- **cv**: Coeffieicnt of variation of the coverage across the region
- **evenness**: Normalized observed entropy
- **F1,10,100,...**: Fraction of the region covered at 1x, 10x, 100x, etc.
