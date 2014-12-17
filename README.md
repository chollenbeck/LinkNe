# LinkNe
#### A program for estimating effective population size from LD at loci with known linkage relationships

### Synopsis

Usage:
     perl LinkNe.pl -i <inputfile> -m <recombination_matrix> [options]


Options:
    -i, --infile
            Genepop input file

    -m, --matfile
            Square matrix of recombination frequencies for loci

    -o, --outfile
            Name of outfile

    -b, --binsize
            Size of bins (in Morgans) for estimating Ne [Default: 0.05]

    -a, --allele_cutoff
            Cutoff frequency for excluding rare alleles from the analysis
            [Default: 0.05]

    -c, --correct_bias
            Correct expected r2 values using the bias correction of Waples
            (2006) - recommended

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

    -r, --recalculate
            Recalculates the moving average based on previously saved data
            in order to avoid the overhead of recalculating all pairwise
            values



### Dependencies

The following perl modules are required for running LinkNe:

Data::Dumper<br />
Getopt::Long<br />
Pod::Usage<br />
Statistics::Distributions<br />