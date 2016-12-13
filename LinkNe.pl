#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Statistics::Distributions;

my $version = '1.0.9';

my $command = 'LinkNe.pl ' . join(" ", @ARGV);

pod2usage(-verbose => 1) if @ARGV == 0;

my $infile = '';
my $outfile = 'Ne.out';
my $matfile = '';
my $mapfile = '';
my $binsize = 0.05;
my $timebin = '';
my $allele_cutoff = 0.05;
my $rec_cutoff = 0;
my $no_bias_corr = '';
my $moving_avg = '';
my $window_size = 0.05;
my $interval = 0.01;
my $save_data = '';
my $reanalyze = '';
my $opt_version = '';
my $debug = '';

GetOptions(	'infile|i=s' => \$infile,
			'outfile|o=s' => \$outfile,
			'matfile|m=s' => \$matfile,
			'map=s' => \$mapfile,
			'bins|b=s' => \$binsize,
			'timebin|t' => \$timebin,
			'allele_cutoff|a=s' => \$allele_cutoff,
			'rec_cutoff|e=s' => \$rec_cutoff,
			'no_bias_corr|c' => \$no_bias_corr,
			'moving_avg|v=s' => \$moving_avg,
			'window|w=s' => \$window_size,
			'interval|n=s' => \$interval,
			'save|s' => \$save_data,
			'reanalyze|r' => \$reanalyze,
			'version' => \$opt_version,
			'debug' => \$debug,
);

if ($opt_version) {
	die "Version ",  $version, "\n";
}

if ($reanalyze) {
	reanalyze($infile, $window_size, $interval, $outfile);
}


open(OUT, ">", $outfile) or die $!;
print OUT join("\t", "MIDPOINT_C", "MEAN_C", "NE", "PARA_95_LOW", "PARA_95_HIGH", "PAIRWISE", "S", "CV", "EMP_95_LOW", "EMP_95_HIGH", "MEAN_R_SQ", "MEAN_EXP_R_SQ", "R_SQ_DRIFT", "R_SQ_DRIFT_LOW", "R_SQ_DRIFT_HIGH"), "\n";


open(R2, ">", $outfile . '.R2.log') or die $!;
open(LOG, ">", $outfile . '.log') or die $!;

if ($debug) {
	open(DUMP, ">", 'dumper.out') or die $!;
}

print R2 join("\t", "loc1", "loc2", "c", "r2"), "\n";

if ($save_data) {
	open(SAV, ">", $infile . '.calc.log') or die $!;
}

# Read in the recombination matrix

my %mat_index;

my @rec_matrix;
if ($mapfile) {
	my ($matrix, $matrix_index) = map_to_matrix($mapfile);
	@rec_matrix = @{$matrix};
	%mat_index = %{$matrix_index};
} elsif ($matfile) {
	open(MAT, "<", $matfile) or die $!;
	my $loc_count = 0;
	while(<MAT>) {
		last if $_ =~ /^\s/;
		my ($locus, @recs) = split;
		push @rec_matrix, \@recs;
		$mat_index{$locus} = $loc_count;
		$loc_count++;
	}
	close MAT;
} else {
	die "You need to provide a recombination matrix (-m) or map file (-map)\n"
}

my ($pops, $pop_loci) = read_genepop($infile, \%mat_index);
#print DUMP Dumper($pops);
my @loci = keys %mat_index;
my @pop_loci = @{$pop_loci};
my @shared_loci;
foreach my $locus (@pop_loci) {
	if (defined $mat_index{$locus}) {
		push @shared_loci, $locus;
	}
}

my @bins; # The data structure that will hold the data for the regular estimates
my @all_bin; # The data structure that will hold the data for the moving average
my @bin_range;


for (my $i = 1; $i <= 0.5 / $binsize; $i++) {
	my $binmax = $i * $binsize;
	my $binmin = $binmax - $binsize;
	my $binmid = ($binmax + $binmin) / 2;
	my @range = ($binmin, $binmid, $binmax);
	push @bin_range, \@range;
}

# Combine bins based on the time period that they represent

if ($timebin) {
	my @new_bins;
	my $prev_time = 0;
	for (my $i = 0; $i < scalar(@bin_range); $i++) {
		my @bin = @{$bin_range[$i]};
		my $binmid = $bin[1];
		my $time = 1 / (2 * $binmid); # Convert to time based on equation from Hayes et al. 2003
		my $rounded_time;
		if ($time =~ /(\d+)\.\d+/) {
			$rounded_time = $1;
		} elsif ($time =~ /^(\d+)$/) {
			$rounded_time = $1;
		} else {
			die "Can't round time to nearest integer\n";
		}
		my $curr_time = $rounded_time;
		if ($i == 0) { # The first bin in the array
			$prev_time = $curr_time;
			push @new_bins, \@bin;
			next;
		} else { # The subsequent bins in the array
			if (abs($curr_time - $prev_time) < 2) {
				# Combine the bins
				my $prev_bin = pop(@new_bins);
				my $new_mid = ($prev_bin->[0] + $bin[2]) / 2;
				my @comb_bin = ($prev_bin->[0], $new_mid, $bin[2]);
				push @new_bins, \@comb_bin;
				$prev_time = $curr_time;
			} else {
				push @new_bins, \@bin;
				$prev_time = $curr_time;
			}
		}
	}
	@bin_range = @new_bins;

}

for (my $i = 0; $i < scalar(@bin_range); $i++) {
	$bins[$i] = [ [],[],[],[],[] ];
}

foreach my $pop (@{$pops}) {

	my %freqs;
	my %exclude;
	my %missing;
	my %sample_size;
	my @filt_loci;
	foreach my $locus (@shared_loci) {
		if (! grep($locus, @loci)) {
			print DUMP "Filtered: $locus\n";
			next;
		}
		my ($freq, $locus_n, $locus_miss) = get_freqs($locus, $pop);

		#print DUMP "All freqs\n";
		#print DUMP Dumper($freq);

		my $excluded = 0;
		# Filter out alleles below a specified frequency cutoff
		foreach my $allele (keys %{$freq}) {
			if ($freq->{$allele} < $allele_cutoff) {
				$exclude{$locus}{$allele} = 1;
				$excluded++;
				#next;
			}
			$freqs{$locus}{$allele} = $freq->{$allele};
		}

		$missing{$locus} = $locus_miss;
		$sample_size{$locus} = $locus_n;


		# Skip the locus if there is only one allele after excluding low-frequency alleles


		if (scalar(keys %{$freqs{$locus}}) - $excluded == 1) {
			next;
		}

		push @filt_loci, $locus;
	}

	print DUMP Dumper(\%freqs) if $debug;

	# Define variables to hold values for each locus pair

	my @loc_burrows;
	my @loc_exp_rsqs;

	my $total_loci = scalar(@filt_loci);
	print "$total_loci loci to test\n";
	print LOG "Total Loci: $total_loci\n";

	# Iterate through each pair of loci
	my $pairwise = 0;
	my $binned = 0;
	for(my $y = 0; $y < scalar(@filt_loci) - 1; $y++) {
		for(my $z = $y + 1; $z < scalar(@filt_loci); $z++) {
			my $processed = $y + 2;
			print "Processed $processed of $total_loci loci\r";

			my @exp_rsqs;
			my @burrows;

			my @pair = ($filt_loci[$y], $filt_loci[$z]);
			print DUMP join("\t", @pair), "\n" if $debug;

			# Look up the recombination frequency for the locus pair

			my $c = $rec_matrix[$mat_index{$pair[0]}][$mat_index{$pair[1]}];

			next if $c eq 'NA';
			next if $c < $rec_cutoff;
			$c = 0.00001 if $c == 0;
			#$c = 0.49999 if $c == 0.5;
			$pairwise++;

			my ($geno_counts, $shared_n, $geno_pos_1, $geno_pos_2, $allele_counts) = get_counts($filt_loci[$y], $filt_loci[$z], $pop, \%freqs);
			my @counts = @$geno_counts;
			my %geno_pos_1 = %$geno_pos_1;
			my %geno_pos_2 = %$geno_pos_2;
			my %allele_counts = %$allele_counts;


			die unless $shared_n > 0; # Trap this exception

			#print DUMP "Allele counts\n";
			#print DUMP Dumper(\@gen_counts);
			#print DUMP Dumper(%{$counts});


			my @alleles1 = sort keys %{$freqs{$pair[0]}};
			my @alleles2 = sort keys %{$freqs{$pair[1]}};


			# Calculate frequencies for alleles in shared individuals

			my %allele_freqs;
			foreach my $locus (@pair) {
				my $total_alleles;
				foreach my $allele (keys %{$allele_counts->{$locus}}) {
					$total_alleles += $allele_counts->{$locus}{$allele};
				}
				foreach my $allele (keys %{$allele_counts->{$locus}}) {
					$allele_freqs{$locus}{$allele} = $allele_counts->{$locus}{$allele} / $total_alleles;
				}

			}

			# Skip the locus pair if one of the loci is monomorphic


			if (scalar(keys %{$allele_freqs{$pair[0]}}) < 2 || scalar(keys %{$allele_freqs{$pair[1]}}) < 2) {
				#print DUMP "Skipping pair: $pair[0] $pair[1]\n";
				next;
			}

			print DUMP Dumper(\%allele_freqs) if $debug;

			#my $loc1_n = $sample_size{$pair[0]};
			#my $loc2_n = $sample_size{$pair[1]};
			my $nij = (scalar(keys %{$allele_freqs{$pair[0]}}) - 1) * (scalar(keys %{$allele_freqs{$pair[1]}}) - 1);
			my $Sij = $shared_n;
			my $wij = $nij * $Sij**2;

			#my @alleles1 = sort keys %{$allele_freqs{$pair[0]}};
			#my @alleles2 = sort keys %{$allele_freqs{$pair[1]}};

			# Iterate through all combinations of alleles between the two loci

			for(my $i = 0; $i < scalar(@alleles1); $i++) {
				if ($exclude{$filt_loci[$y]}{$alleles1[$i]}) {
					print LOG "Excluded ", $filt_loci[$y], ' ', $alleles1[$i], "\n";
					next;
				}
				for(my $j = 0; $j < scalar(@alleles2); $j++) {
					next if $exclude{$filt_loci[$z]}{$alleles2[$j]};

					# Make a 3 x 3 genotype matrix for the allele pair, collapsing all other alleles into a single allele

					my @matrix = ( [0, 0, 0], [0, 0, 0], [0, 0, 0] );

					for (my $k = 0; $k < scalar(@counts); $k++) { # $k indexes rows (locus 1 genotypes)
						for (my $l = 0; $l < @{$counts[0]}; $l++) { # $l indexes columns (locus 2 genotypes)
							#my $val = $counts[$k][$l];
							# Determine which row the value should go into
							my $row;
							my $geno1 = $geno_pos_1{$k};
							if ($geno1 eq "$alleles1[$i]$alleles1[$i]") { # a homozygote
								$row = 0;
							} elsif (substr($geno1, 0, length($alleles1[$i])) eq $alleles1[$i] || substr($geno1, length($alleles1[$i]), length($alleles1[$i])) eq $alleles1[$i]) { # heterozygote
								$row = 1;
							} else { # the allele is not in the genotype
								$row = 2;
							}
							# Determine which column the value should go into
							my $col;
							my $geno2 = $geno_pos_2{$l};
							if ($geno2 eq "$alleles2[$j]$alleles2[$j]") { # a homozygote
								$col = 0;
							} elsif (substr($geno2, 0, length($alleles2[$j])) eq $alleles2[$j] || substr($geno2, length($alleles2[$j]), length($alleles2[$j])) eq $alleles2[$j]) { # heterozygote
								$col = 1;
							} else { # the allele is not in the genotype
								$col = 2;
							}
							#print COUNT "Placing value in coordinates ($row, $col)\n";
							$matrix[$row][$col] += $counts[$k][$l];

						}
					}


					my $hi = ($matrix[0][0] + $matrix[0][1] + $matrix[0][2]) / $shared_n;
					my $hj = ($matrix[0][0] + $matrix[1][0] + $matrix[2][0]) / $shared_n;
					my $hets = (2 * $matrix[0][0]) + $matrix[0][1] + $matrix[1][0] + ($matrix[1][1] / 2);


					my $p = $allele_counts->{$pair[0]}{$alleles1[$i]} / (2 * $shared_n);
					my $q = $allele_counts->{$pair[1]}{$alleles2[$j]} / (2 * $shared_n);

					if ($debug) {
						print DUMP "$alleles1[$i] $alleles2[$j]\n";
						print DUMP "Hets: $hets\n";
						print DUMP "n: $shared_n\n";
						print DUMP "p: $p\n";
						print DUMP "q: $q\n";
						print DUMP "hi: $hi\n";
						print DUMP "hj: $hj\n";
						print DUMP "nij: $nij\n";
						print DUMP "Sij: $Sij\n";
					}

					my $D_hat = ($shared_n / ($shared_n - 1)) * (($hets / $shared_n) - (2 * $p * $q));


					# This equation causes a division by zero error under certain circumstances (such as when an allele frequency
					# is exactly 0.5 and there are no corresponding homozygotes). This happens sometimes with small sample sizes.
					# This traps the error and skips the allele pair

					my $r = '';
					eval { $r = $D_hat / ( sqrt( (($p*(1-$p))+($hi-$p**2)) * (($q*(1-$q))+($hj-$q**2)) ) ) };
					next if $@;

					$r = 1 if $r > 1 || $r < -1;



					# if ($@) {

						# print DUMP "$pair[0] $pair[1]\n";
						# print DUMP "$alleles1[$i] $alleles2[$j]\n";
						# print DUMP "Hets: $hets\n";
						# print DUMP "n: $shared_n\n";
						# print DUMP "p: $p\n";
						# print DUMP "q: $q\n";
						# print DUMP "hi: $hi\n";
						# print DUMP "hj: $hj\n";
						# print DUMP "nij: $nij\n";
						# print DUMP "Sij: $Sij\n";

					# }
					my $r_sq = $r**2;

					if ($debug) {
						print DUMP "D: $D_hat\n";
						print DUMP "r: $r\n";
						print DUMP "r^2: $r_sq\n";
					}

					my $exp_rsq;
					unless ($no_bias_corr) {
						if ($Sij >= 30) {
							$exp_rsq = (1 / $Sij) + (3.19 / $Sij**2);
						} else {
							$exp_rsq = 0.0018 + (0.907 / $Sij) + (4.44 / $Sij**2);
						}
					} else {
						$exp_rsq = (1 / $Sij);
					}

					push @burrows, $r_sq;

					push @exp_rsqs, $exp_rsq;

				}
			}
			next if scalar(@burrows) == 0;

			# Calculate the means for the locus pair

			my $total_loc_rsq;
			my $total_loc_exp_rsq;
			for(my $i = 0; $i < scalar(@burrows); $i++) {
				$total_loc_rsq += $burrows[$i];
				$total_loc_exp_rsq += $exp_rsqs[$i];
			}
			my $mean_loc_rsq = $total_loc_rsq / scalar(@burrows);
			my $mean_loc_exp_rsq = $total_loc_exp_rsq / scalar(@exp_rsqs);

			if ($moving_avg) {
				push @all_bin, [ $mean_loc_rsq, $mean_loc_exp_rsq, $nij, $Sij, $wij, $c ];
			}
			if ($save_data) {
				print SAV join("\t", $pair[0], $pair[1], $mean_loc_rsq, $mean_loc_exp_rsq, $nij, $Sij, $wij, $c), "\n";
			}
			if ($debug) {
				my $locus_rsq_drift = $mean_loc_rsq - $mean_loc_exp_rsq;
				print DUMP "Mean Locus r^2: $mean_loc_rsq\n";
				print DUMP "Mean Locus r^2-drift: $locus_rsq_drift\n";
			}

			my $binable = 0;
			for (my $g = 0; $g < scalar(@bins); $g++) {

				my $bin_min = $bin_range[$g][0];
				my $bin_max = $bin_range[$g][2];

				if ($c > $bin_min) {
					if ($c <= $bin_max) {
						push @{$bins[$g][0]}, $mean_loc_rsq;
						push @{$bins[$g][1]}, $mean_loc_exp_rsq;
						push @{$bins[$g][2]}, $nij;
						push @{$bins[$g][3]}, $Sij;
						push @{$bins[$g][4]}, $wij;
						push @{$bins[$g][5]}, $c;
						$binned++;
						$binable = 1;
						print R2 join("\t", $pair[0], $pair[1], $c, $mean_loc_rsq - $mean_loc_exp_rsq), "\n";
						#print "Binned!\n";
						last;
					}
				}
			}

			if ($binable == 0) {
				#print DUMP "Unable to bin: $pair[0], $pair[1]\n";
				#print DUMP "c: $c\n";
				die "Locus pair unable to be binned\n";
			}

		}
	}

	# Calculate the Ne for each bin

	# i indexes the bin, j indexes locus pair within the bin
	open(CVAL, ">", 'cvalues.txt');
	for(my $i = 0; $i < scalar(@bins); $i++) {
		print CVAL "BIN\n";
		# Calculate N (total independent comparisons) and S (weighted harmonic mean of sample size) for each bin
		my $N;
		my $N_over_S;
		for (my $j = 0; $j < scalar(@{$bins[$i][0]}); $j++) {
			$N += $bins[$i][2][$j]; # sum the nij values to compute N
			$N_over_S += $bins[$i][2][$j] / $bins[$i][3][$j]; # sum (nij / Sij) to compute N / S
		}

		my $S = $N / $N_over_S;

		my $total_W;
		my $total_weighted_r_sq;
		my $total_weighted_exp_r_sq;
		my $total_c;

		for(my $j = 0; $j < scalar(@{$bins[$i][0]}); $j++) {
			my $product = $bins[$i][4][$j] * $bins[$i][0][$j]; # Mean locus r_sq * corresponding wij
			$total_weighted_r_sq += $product; # Add the product to the total weighted r_sq

			my $exp_product = $bins[$i][4][$j] * $bins[$i][1][$j]; # Mean locus exp_r_sq * corresponding wij
			$total_weighted_exp_r_sq += $exp_product; # Add the product to the total weighted exp_r_sq

			$total_W += $bins[$i][4][$j]; # Add the wij to the total W
			$total_c += $bins[$i][5][$j]; # Add the cij to the total c
			print CVAL $bins[$i][5][$j], "\n";
		}
		#print join("\n", @wijs), "\n";
		#print join("\n", @exp_rsqs), "\n";

		my $mean_exp_r_sq = $total_weighted_exp_r_sq / $total_W;
		my $mean_r_sq = $total_weighted_r_sq / $total_W;
		my $mean_c = $total_c / scalar(@{$bins[$i][0]});
		#print DUMP join("\t", 'r^2', $mean_r_sq), "\n";

		my $r_sq_drift = $mean_r_sq - $mean_exp_r_sq;

		my $midpoint_c = $bin_range[$i][1];

		# Calculate a gamma value for the bin, where c is the mean of the bin

		my $gamma = ((1 - $mean_c)**2 + $mean_c**2) / (2 * $mean_c * (2 - $mean_c));

		# From Hill 1981 / Waples 2006
		my $Ne = ($gamma / $r_sq_drift);

		# Tenesa et al. 2007 approximation
		#my $Ne = ((1 / $r_sq_drift) - 1) / (4 * $bin_means[$i]);

		# Calculate parametric CIs based on the LDNe approach

		my $chis_low = Statistics::Distributions::chisqrdistr($N, 0.025);
		my $chis_high = Statistics::Distributions::chisqrdistr($N, 0.975);
		my $r_sq_low = ($N * $mean_r_sq) / $chis_low;
		my $r_sq_high = ($N * $mean_r_sq) / $chis_high;
		my $r_sq_drift_low = $r_sq_low - $mean_exp_r_sq;
		my $r_sq_drift_high = $r_sq_high - $mean_exp_r_sq;

		my $Ne_high = $gamma / $r_sq_drift_low;
		my $Ne_low = $gamma / $r_sq_drift_high;

		# Calculate CV from Hill (1981)

		my $CV = (1 + ($Ne / ($gamma * $S))) * sqrt(2 / $N);
		my $SD = $CV * $Ne;
		my $rough_low = $Ne - (1.98 * $SD);
		my $rough_high = $Ne + (1.98 * $SD);

		print OUT join("\t", $midpoint_c, $mean_c, $Ne, $Ne_low, $Ne_high, scalar(@{$bins[$i][0]}), $S, $CV, $rough_low, $rough_high, $mean_r_sq, $mean_exp_r_sq, $r_sq_drift, $r_sq_drift_low, $r_sq_drift_high), "\n";
		#print OUT join("\t", $mean_c, $Ne, $Ne_low, $Ne_high, scalar(@{$bins[$i][0]}), $S, $CV, $rough_low, $rough_high, $mean_r_sq, $mean_exp_r_sq, $r_sq_drift, $r_sq_drift_low, $r_sq_drift_high), "\n";

	}

	if ($moving_avg) {
		calc_moving_avg(\@all_bin, $window_size, $interval, $moving_avg);
	}

}

sub read_genepop {
	my $file = $_[0];
	my %mat_index = %{$_[1]};

	my $commas = 0;
	open (TEST, "<", $file) or die $!;
	<TEST>;
	my $locus_string = <TEST>;
	if ($locus_string =~ /,/) {
		$commas = 1;
	}
	close TEST;

	open (IN, "<", $file) or die;

	my @loci;
	chomp(my $title = <IN>);

	if ($commas == 1) {
		my $locus_string = <IN>;
		@loci = split ",", $locus_string;
		foreach my $locus (@loci) {
			$locus =~ s/\s//g;
		}
		my $first_pop = <IN>;
	} else {
		while (<IN>) {
			if ($_ =~ /pop/i) { last; }
			chomp;
			push @loci, $_;
		}
	}

	my @pops = [];
	my $pop_counter = 0;
	$pops[$pop_counter] = {};
	while (<IN>) {
		if ($_ =~ /^pop/i) {
			$pop_counter++;
			$pops[$pop_counter] = {};
			next;
		}
		my ($ind, $genotypes) = split(/,/, $_);
		$ind =~ s/\s//g;
		$pops[$pop_counter]{$ind} = {};
		my $locus_no = 0;
		foreach my $genotype (split(/\s/, $genotypes)) {
			#print "-$genotype-";
			my $width = (length $genotype) / 2;
			next if $width == 0;
			my $allele1 = substr($genotype, 0, $width);
			my $allele2 = substr($genotype, $width, $width);
			if (! $mat_index{$loci[$locus_no]}) {
				$locus_no++;
				#print "Skipping\n";
				next;
			}
			$pops[$pop_counter]{$ind}{$loci[$locus_no]} = [$allele1, $allele2];
			$locus_no++;
		}
	}
	close IN;
	return (\@pops, \@loci);
}

sub get_freqs {
	my $locus = $_[0];
	my %pop = %{$_[1]};

	my $loc_miss;
	my $total = 0;
	my %freq;
	my %allele_count;
	foreach my $ind (keys %pop) {
		my @alleles = @{$pop{$ind}{$locus}};
		if ($alleles[0] eq '000' || $alleles[1] eq '000') { # This should be updated to account for the two allele case
			$loc_miss++;
			next;
		}
		$allele_count{$alleles[0]}++;
		$allele_count{$alleles[1]}++;
		$total += 2;
	}
	foreach my $allele (keys %allele_count) {
		$freq{$allele} = $allele_count{$allele} / $total;
	}
	my $loc_n = $total / 2;

	return (\%freq, $loc_n, $loc_miss);
}

sub get_counts {
	my $locus1 = $_[0];
	my $locus2 = $_[1];
	my %pop = %{$_[2]};
	my %freqs = %{$_[3]};

	my @loci = ($locus1, $locus2);

	my @alleles1 = sort keys %{$freqs{$locus1}};
	my @alleles2 = sort keys %{$freqs{$locus2}};

	my $num_gts_1 = (scalar(@alleles1) * (1 + scalar(@alleles1))) / 2;
	my $num_gts_2 = (scalar(@alleles2) * (1 + scalar(@alleles2))) / 2;

	# Assign each allele a number, for sorting purposes

	my %allele_1_id;
	my %allele_2_id;
	for(my $i = 0; $i < scalar(@alleles1); $i++) {
		$allele_1_id{$alleles1[$i]} = $i;
	}
	for(my $i = 0; $i < scalar(@alleles2); $i++) {
		$allele_2_id{$alleles2[$i]} = $i;
	}

	# Set up an empty matrix to record genotype counts

	my @counts;
	for(my $i = 0; $i < $num_gts_1; $i++) {
		for(my $j = 0; $j < $num_gts_2; $j++) {
			$counts[$i][$j] = 0;
		}
	}

	# Create a hash that specifies the position of each genotype in the matrix
	my %geno_pos_1;
	my %geno_pos_2;
	my $pos = 0;
	for(my $i = 0; $i < scalar(@alleles1); $i++) {
		for(my $j = $i; $j < scalar(@alleles1); $j++) {
			$geno_pos_1{"$alleles1[$i]$alleles1[$j]"} = $pos;
			$pos++;
		}
	}
	$pos = 0;
	for(my $i = 0; $i < scalar(@alleles2); $i++) {
		for(my $j = $i; $j < scalar(@alleles2); $j++) {
			$geno_pos_2{"$alleles2[$i]$alleles2[$j]"} = $pos;
			$pos++;
		}
	}

	# Iterate through individuals, counting genotypes

	my %allele_counts;
	my $shared_n = 0;
	foreach my $ind (keys %pop)	{
		next if $pop{$ind}{$locus1}[0] eq '000' || $pop{$ind}{$locus2}[0] eq '000';

		my @alleles_1 = sort { $allele_1_id{$a} <=> $allele_1_id{$b} } @{$pop{$ind}{$locus1}};
		my @alleles_2 = sort { $allele_2_id{$a} <=> $allele_2_id{$b} } @{$pop{$ind}{$locus2}};

		#next if $alleles_1[0] eq '000' || $alleles_2[0] eq '000';

		my $geno1 = "$alleles_1[0]$alleles_1[1]";
		my $geno2 = "$alleles_2[0]$alleles_2[1]";

		# Update the allele counts for the two loci
		$allele_counts{$locus1}{$alleles_1[0]}++;
		$allele_counts{$locus1}{$alleles_1[1]}++;
		$allele_counts{$locus2}{$alleles_2[0]}++;
		$allele_counts{$locus2}{$alleles_2[1]}++;

		# Find the number of the row in the matrix where the genotype of the first locus is counted
		my $row = $geno_pos_1{$geno1};

		# Find the number of the column in the matrix where the genotype of the second locus is counted
		my $col = $geno_pos_2{$geno2};

		$counts[$row][$col]++;
		$shared_n++;
	}



	%geno_pos_1 = reverse %geno_pos_1;
	%geno_pos_2 = reverse %geno_pos_2;

	return \@counts, $shared_n, \%geno_pos_1, \%geno_pos_2, \%allele_counts;

}

sub calc_moving_avg {
	my @data = @{$_[0]};
	my $window_size = $_[1];
	my $interval = $_[2];
	my $out_file = $_[3];

	my @sorted = sort { $a->[5] <=> $b->[5] } @data;


	#print TEST Dumper(\@sorted);


	my $max = 0;
	my $last_start = 0;
	my $last_end = 0;

	open(MAVG, ">", $out_file) or die $!;


	for (my $min = 0; $max < 0.5; $min += $interval) {
		$max = $min + $window_size;
		my @bin;
		foreach my $pair (@sorted) {
			if ($pair->[5] > $min) {
				if ($pair->[5] <= $max) {
					push @bin, $pair;
				}
			}
			#print TEST Dumper($pair) unless $pair->[5];
		}

		# Calculate N (total independent comparisons) and S (weighted harmonic mean of sample size) for each bin
		my $N;
		my $N_over_S;

		foreach my $pair (@bin) {
			$N += $pair->[2];
			$N_over_S += ($pair->[2] / $pair->[3]);
		}

		my $S = $N / $N_over_S;


		#print "Calculating for bin $i\n";
		my $total_W;
		my $total_weighted_r_sq;
		my $total_weighted_exp_r_sq;
		my $total_c;

		foreach my $pair (@bin) {
			my $product = $pair->[0] * $pair->[4]; # Mean locus r_sq * corresponding wij
			$total_weighted_r_sq += $product; # Add the product to the total weighted r_sq

			my $exp_product = $pair->[1] * $pair->[4]; # Mean locus exp_r_sq * corresponding wij
			$total_weighted_exp_r_sq += $exp_product; # Add the product to the total weighted exp_r_sq

			$total_W += $pair->[4]; # Add the wij to the total W
			$total_c += $pair->[5]; # Add the cij to the total c
		}


		my $mean_exp_r_sq = $total_weighted_exp_r_sq / $total_W;
		my $mean_r_sq = $total_weighted_r_sq / $total_W;
		#print DUMP join("\t", 'r^2', $mean_r_sq), "\n";

		my $r_sq_drift = $mean_r_sq - $mean_exp_r_sq;

		my $midpoint = ($max + $min) / 2;
		my $mean_c = $total_c / scalar(@bin);

		# Calculate a gamma value for the bin, where c is the midpoint of the bin
		#my $gamma = ((1 - $midpoint)**2 + $midpoint**2) / (2 * $midpoint * (2 - $midpoint));

		my $gamma = ((1 - $mean_c)**2 + $mean_c**2) / (2 * $mean_c * (2 - $mean_c));

		# From Hill 1981 / Waples 2006
		my $Ne = ($gamma / $r_sq_drift);

		# Tenesa et al. 2007 approximation
		#my $Ne = ((1 / $r_sq_drift) - 1) / (4 * $bin_means[$i]);

		# Calculate parametric CIs based on the LDNe approach

		my $chis_low = Statistics::Distributions::chisqrdistr($N, 0.025);
		my $chis_high = Statistics::Distributions::chisqrdistr($N, 0.975);
		my $r_sq_low = ($N * $mean_r_sq) / $chis_low;
		my $r_sq_high = ($N * $mean_r_sq) / $chis_high;
		my $r_sq_drift_low = $r_sq_low - $mean_exp_r_sq;
		my $r_sq_drift_high = $r_sq_high - $mean_exp_r_sq;

		my $Ne_high = $gamma / $r_sq_drift_low;
		my $Ne_low = $gamma / $r_sq_drift_high;

		# Calculate CV from Hill (1981)

		my $CV = (1 + ($Ne / ($gamma * $S))) * sqrt(2 / $N);
		my $SD = $CV * $Ne;
		my $rough_low = $Ne - (2 * $SD);
		my $rough_high = $Ne + (2 * $SD);

		#print MAVG join("\t", $midpoint, $Ne), "\n";
		print MAVG join("\t", $midpoint, $mean_c, $Ne, $Ne_low, $Ne_high, scalar(@bin), $S, $CV, $rough_low, $rough_high, $mean_r_sq, $mean_exp_r_sq, $r_sq_drift, $r_sq_drift_low, $r_sq_drift_high), "\n";


	}
}

sub reanalyze {
	my $file = $_[0];
	my $window_size = $_[1];
	my $interval = $_[2];
	my $out_file = $_[3];

	open(IN, "<", $file) or die $!;
	#open(TEST, ">", 'pairs.txt');

	my @all_data;
	while(<IN>) {
		last if $_ =~ /^\s/;
		chomp;
		my ($loc1, $loc2, $mean_loc_rsq, $mean_loc_exp_rsq, $nij, $Sij, $wij, $c) = split;
		push @all_data, [ $mean_loc_rsq, $mean_loc_exp_rsq, $nij, $Sij, $wij, $c ];
	}
	close IN;
	print "Finished reading data\n";
	#print TEST Dumper(\@all_data);
	calc_moving_avg(\@all_data, $window_size, $interval, $out_file);

	die "Finished reanalyzing data\n";
}

sub map_to_matrix {

	$mapfile = shift;

	open(MAP, "<", $mapfile) or die $!;

	my %map_pos;
	my %map_lg;
	my @loci;
	<MAP>;
	while(<MAP>) {
		chomp;
		my ($locus, $lg, $pos) = split;
		$map_pos{$locus} = $pos;
		$map_lg{$locus} = $lg;
		push @loci, $locus;
	}
	close MAP;

	my @mat;
	for (my $i = 0; $i < scalar(@loci); $i++) {
		for (my $j = 0; $j < scalar(@loci); $j++) {
			my $loc1 = $loci[$i];
			my $loc2 = $loci[$j];

			# Calculate the distance between loci
			if ($map_lg{$loc1} != $map_lg{$loc2}) {
				$mat[$i][$j] = 0.5;
				$mat[$j][$i] = 0.5;
				next;
			} else {
				my $dist = abs($map_pos{$loc1} - $map_pos{$loc2});
				if ($dist >= 50) {
					$mat[$i][$j] = 0.5;
					$mat[$j][$i] = 0.5;
					next;
				} else {
					$mat[$i][$j] = $dist / 100;
					$mat[$j][$i] = $dist / 100;
				}
			}
		}
	}

	my %mat_index;
	my @final_mat;
	for (my $i = 0; $i < scalar(@loci); $i++) {
		my $row = [ $loci[$i], @{$mat[$i]} ];
		push @final_mat, $row;
		$mat_index{$loci[$i]} = $i;
		#print join("\t", @$row), "\n";
	}


	return \@final_mat, \%mat_index;

}




__END__

=head1 NAME

LinkNe.pl

=head1 SYNOPSIS

perl LinkNe.pl -i <inputfile> -m <recombination_matrix> -o <outfile> -b <binsize> [options]

Options:
     -i	<inputfile>		input genepop file

	 -map	<linkage_map>		tab separated linkage map

	 -m	<recombination_matrix>		matrix of recombination frequencies for loci

	 -o	<outfile>		name of output file

	 -b	<binsize>		size of bins (in Morgans) for estimating Ne

	 -t		bin estimates by time (1 / 2c) rather than by recombination frequency

	 -a	<allele_cutoff>		cutoff frequency for excluding rare alleles from the analysis

	 -e	<rec_cutoff>		exclude locus pairs below a specified recombination faction

	 -c		correct expected r2 values using the bias correction of Waples (2006) - recommended

	 -v		compute a moving average for effective size relative to recombination rate

	 -w <window_size>		size of sliding window for calculating moving average

	 -n <interval>		distance the sliding window moves for each calculation

	 -s		save data to a file for later use

	 -r		recalculate a moving average using previously save data




=head1 OPTIONS

=over 8

=item B<-i, --infile>

Genepop input file

=item B<-map>

Linkage map in tab-separated format:
<locus>	<chromosome>	<position>

Example:
locus	chromosome	position
loc1	1	5.0
loc2	1	15.0
loc3	2	50.5

Position values should be in centiMorgans

=item B<-m, --matfile>

Matrix of recombination frequencies for loci

=item B<-o, --outfile>

Outfile

=item B<-b, --binsize>

Size of bins (in Morgans) for estimating Ne [Default: 0.05]

=item B<-t, --timebin>

Bin pairwise estimates of LD by generations (measured by 1 / 2c) rather than by recombination frequency (c) itself

=item B<-a, --allele_cutoff>

Cutoff frequency for excluding rare alleles from the analysis [Default: 0.05]

=item B<-e, --rec_cutoff>

exclude locus pairs below a specified recombination faction [Default: No cutoff]

=item B<-c, --no_bias_corr>

Turn off correction of expected r2 values using the bias correction of Waples (2006) - not recommended except for experimental use

=item B<-v, --moving_avg>

Compute a moving average for effective size relative to recombination rate. Must specify a window size and interval for moving the window.

=item B<-w, --window>

Size of sliding window for calculating moving average (in Morgans) [Default: 0.01]

=item B<-n, --interval>

Distance the sliding window moves for each calculation (in Morgans) [Default: 0.005]

=item B<-s, --save>

Saves relevant pairwise data to a file for later recalculation of the moving average


=back

=head1 DESCRIPTION

B<LinkNe.pl> estimates effective population size from linkage disequilibrium using linked and unlinked loci

=cut
