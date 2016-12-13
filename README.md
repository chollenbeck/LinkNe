# LinkNe
#### A program for estimating effective population size from LD at loci with known linkage relationships

### Synopsis

Usage:
     perl LinkNe.pl -i <inputfile> -m <recombination_matrix> [options]


     Options:

         -i, --infile
                 Genepop input file - currently supports three-character alleles

         -map    Linkage map in tab-separated format: <locus> <chromosome>
                 <position>

                 Example: locus chromosome position loc1 1 5.0 loc2 1 15.0 loc3 2
                 50.5

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
