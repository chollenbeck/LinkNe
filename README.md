# LinkNe

### Overview

LinkNe is a program for estimating effective population size from LD at loci with known linkage relationships. A detailed description of the program can be found in:

Hollenbeck CM, Portnoy DS, Gold JR (2016) A method for detecting recent changes in contemporary effective population size from linkage disequilibrium at linked and unlinked loci. Heredity, 117, 207–216. [[link](http://www.nature.com/hdy/journal/v117/n4/abs/hdy201630a.html)]

### Input

Requirements:

1) Genepop file of diploid genotypes
    - the file should contain only one population
    - only the three-character allele format is currently supported
    - missing genotypes are indicated with '000000'

2) Recombination rate estimates for at least a subset of the genetic markers. Two formats are allowed:

- A genetic map with the following format (note that map positions are in cM and the header is required):

```
    locus chromosome  position
    loc1  1 5.0
    loc2  1 15.0
    loc3  2 50.5
```

- A square matrix of recombination rates (note that recombination rates are in Morgans). Example (equivalent to the map above):

```
    loc1  0 0.05 0.50
    loc2  0.05 0 0.50
    loc3  0.50 0.50 0
```

### Output

#### Ne.out

The columns in the Ne.out file correspond to the following:

- MIDPOINT_C: The midpoint recombination rate of the bin
- MEAN_C: The mean recombination rate of the bin
- NE: The point estimate of effective population size
- PARA_95_LOW: Lower bound of 95% parametric confidence interval
- PARA_95_HIGH: Upper bound of 95% parametric confidence interval
- PAIRWISE: Number of pairwise locus comparisons used in the estimate
- S: Mean sample size
- CV: Coefficient of variation - from Hill (1981)
- EMP_95_LOW: Lower bound of 95% confidence interval based on CV from Hill (1981)
- EMP_95_HIGH: Upper bound of 95% confidence interval based on CV from Hill (1981)
- MEAN_R_SQ: Mean total r-square
- MEAN_EXP_R_SQ: Mean expected r-square
- R_SQ_DRIFT: R-square attributable to drift
- R_SQ_DRIFT_LOW: Lower bound of 95% confidence interval for r-square drift
- R_SQ_DRIFT_HIGH: Upper bound of 95% confidence interval for r-square drift

### Usage


    perl LinkNe.pl -i input_file -map linkage_map [options]


     Options:

         -i, --infile
                 Genepop input file - currently only supports three-character alleles where '000000' indicates a missing genotype

         -map    Linkage map in tab-separated format (with header):

                <locus> <chromosome>  <position>

                Example:

                locus chromosome  position
                loc1  1 5.0
                loc2  1 15.0
                loc3  2 50.5

                Position values should be in centiMorgans

         -m, --matfile
                 Matrix of recombination frequencies for loci

         -o, --outfile
                 Outfile

         -b, --binsize
                 Size of bins (in Morgans) for estimating Ne [Default: 0.05]

         -t, --timebin
                 Bin pairwise estimates of LD by generations (measured by 1 / 2c)
                 rather than by recombination frequency (c) itself

         -a, --allele_cutoff
                 Cutoff frequency for excluding rare alleles from the analysis
                 [Default: 0.05]

         -e, --rec_cutoff
                 exclude locus pairs below a specified recombination faction
                 [Default: No cutoff]

         -c, --no_bias_corr
                 Turn off correction of expected r2 values using the bias
                 correction of Waples (2006) - not recommended except for
                 experimental use

         -v, --moving_avg
                 Compute a moving average for effective size relative to
                 recombination rate. Must specify a window size and interval for
                 moving the window.

         -w, --window
                 Size of sliding window for calculating moving average (in
                 Morgans) [Default: 0.01]

         -n, --interval
                 Distance the sliding window moves for each calculation (in
                 Morgans) [Default: 0.005]

         -s, --save
                 Saves relevant pairwise data to a file for later recalculation
                 of the moving average


### Dependencies

The following perl modules are required for running LinkNe:

Data::Dumper<br />
Getopt::Long<br />
Pod::Usage<br />
Statistics::Distributions<br />
