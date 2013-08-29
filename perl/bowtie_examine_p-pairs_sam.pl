#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 3/21/12


# This scripts compares aligned reads to their p-pair, and verifies that they
# are on the same chromosome, within 1 kb of each other, and facing each other.

# The loops have been arranged such that entries only remain in the hash until
# their pair is found.

################################################################################

my %full_hash;
my %strand_hash;
my %chrom_hash;
my %start_hash;
my $difference;


unless ( @ARGV == 1 ) {
	die "Usage: bowtie_examine_pairs_v1.4.pl default_bowtie_output\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}


unless ( open INPUT, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

unless ( open DIFCHROM, ">", "$ARGV[0]_difchrom" ) {
	die "Cannot open file $ARGV[0]_difchrom: $!";
}

unless ( open TOOFAR, ">", "$ARGV[0]_toofar" ) {
	die "Cannot open file $ARGV[0]_toofar: $!";
}

unless ( open BADORIENT, ">", "$ARGV[0]_badorientation" ) {
	die "Cannot open file $ARGV[0]_badorientation: $!";
}

# OUTPUT is just the 'real' reads, not the |p reads.
unless ( open OUTPUT, ">", "$ARGV[0]_nofarpair" ) {
	die "Cannot open file $ARGV[0]_nofarpair: $!";
}

# BOTHPAIR is both the 'real' reads and the |p reads.
unless ( open BOTHPAIR, ">", "$ARGV[0]_nofarpair_with-p" ) {
	die "Cannot open file $ARGV[0]_nofarpair: $!";
}

# Read in the bowtie file
while (<INPUT>) {

	chomp;
	my $full_line = $_;
	# This version splits on whitespace and discards leading whitespace (just in
	# case)
	my @line = split;
	
	# Put the whole line, and various important parts of the line, in their
	# respective hashes, always using the name of the read as the key.
	$full_hash{$line[0]} = $full_line;
	$strand_hash{$line[0]} = $line[1];
	$chrom_hash{$line[0]} = $line[2];
	$start_hash{$line[0]} = $line[3];
	
	# Bowtie reports alignments beginning at the leftmost base.
	# If the read is -, change the start position to be the rightmost base.
	# That way, the comparisons will be based on the 5' base of each.
	if ($strand_hash{$line[0]} == 16) {
		
		# The '-1' is necessary to make accurate comparisons of position.
		$start_hash{$line[0]} += length ($line[9]) - 1;
		
	}
	
	my $name = $line[0];
	my $pair_name = $name;
	
	# Note this 'if' statement is different than in other scripts; it's set such
	# that $name will end in 1 or 2, and $pair_name will end in |p.
	if ( substr($pair_name, -1) eq '1' ) {
		substr($pair_name, -1) = '2|p';
	} elsif ( substr($pair_name, -1) eq '2' ) {
		substr($pair_name, -1) = '1|p';
	} elsif ( substr($pair_name, -3) eq '1|p') {
		substr($name, -3) = '2';
	} elsif ( substr($pair_name, -3) eq '2|p') {
		substr($name, -3) = '1';
	} else {
		next;
	}
	
	# If both pairs are in the hash, compare them.
	if ( exists $full_hash{$pair_name} && exists $full_hash{$name}) {
		
		# Test if they're on different chromosomes
		if ($chrom_hash{$pair_name} ne $chrom_hash{$name}) {
			
			# If they're on different chromosomes, print to DIFCHROM.
			print DIFCHROM "$name is on $chrom_hash{$name}; $pair_name is on $chrom_hash{$pair_name}\n";
			
			# Remove both entries from the hashes.
			delete $full_hash{$pair_name};
			delete $strand_hash{$pair_name};
			delete $chrom_hash{$pair_name};
			delete $start_hash{$pair_name};
			
			delete $full_hash{$name};
			delete $strand_hash{$name};
			delete $chrom_hash{$name};
			delete $start_hash{$name};
			
			next;
		
		} else {
			
			# If pairs are on the same chromosome, see how far apart they are.
			if ($start_hash{$pair_name} > $start_hash{$name}) {
				$difference = $start_hash{$pair_name} - $start_hash{$name};
			} else {
				$difference = $start_hash{$name} - $start_hash{$pair_name};
			}
			
			# If they're too far apart (1 kb or more), get rid of them.
			if ($difference >= 1000) {
				
				# Note the distance reported is between the 5' bases of  the 
				# aligned portion of each read.
				print TOOFAR "$name and $pair_name are $difference bp apart\n";
				
				delete $full_hash{$pair_name};
				delete $strand_hash{$pair_name};
				delete $chrom_hash{$pair_name};
				delete $start_hash{$pair_name};
				
				delete $full_hash{$name};
				delete $strand_hash{$name};
				delete $chrom_hash{$name};
				delete $start_hash{$name};
				
				next;
			
			} else {
				
				# If the pairs are less than 1 kb apart, check if they point
				# toward each other.
				if ($start_hash{$pair_name} <= $start_hash{$name} && $strand_hash{$pair_name} == 0 && $strand_hash{$name} == 16 ) {
						
					# If the pair came first or the same, is +, and the other 
					# name is -, it's good.
					print OUTPUT "$full_hash{$name}\n";
					print BOTHPAIR "$full_hash{$name}\n$full_hash{$pair_name}\n";
			
					delete $full_hash{$pair_name};
					delete $strand_hash{$pair_name};
					delete $chrom_hash{$pair_name};
					delete $start_hash{$pair_name};
					
					delete $full_hash{$name};
					delete $strand_hash{$name};
					delete $chrom_hash{$name};
					delete $start_hash{$name};	
							
				} elsif ($start_hash{$pair_name} >= $start_hash{$name} && $strand_hash{$name} == 0 && $strand_hash{$pair_name} == 16 ) {
					
					# If the main name is first or the same, is +, and the pair 
					# is -, it's good.
					print OUTPUT "$full_hash{$name}\n";
					print BOTHPAIR "$full_hash{$name}\n$full_hash{$pair_name}\n";
					
					delete $full_hash{$pair_name};
					delete $strand_hash{$pair_name};
					delete $chrom_hash{$pair_name};
					delete $start_hash{$pair_name};
					
					delete $full_hash{$name};
					delete $strand_hash{$name};
					delete $chrom_hash{$name};
					delete $start_hash{$name};
				
				} else {
					
					# Otherwise, it's no good.
					print BADORIENT "$name and $pair_name do not face each other\n";
					
					delete $full_hash{$pair_name};
					delete $strand_hash{$pair_name};
					delete $chrom_hash{$pair_name};
					delete $start_hash{$pair_name};
					
					delete $full_hash{$name};
					delete $strand_hash{$name};
					delete $chrom_hash{$name};
					delete $start_hash{$name};
					
					next;
				
				}
			
			}
			
		}
		
	} else {
		
		# If you can't find the pair in the hash, you just may not have 
		# encountered it yet.  Keep the current line in the hash and proceed.
		next;
			
	}

}

close BOTHPAIR;

close OUTPUT;

close BADORIENT;

close TOOFAR;

close DIFCHROM;

close INPUT;

# There are presumably still things left in the hash; specifically, unpaired
# reads.  You can count them up, here.  Note, however, that if a read does not
# have a p-pair that made it through the alignment step, then that read will
# NOT make it past this step.

my $pairless = keys %full_hash;

print "There are $pairless names without a paired end.\n";



my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
