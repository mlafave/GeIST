#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 5/14/12
# Last revised: 5/18/12

# Input is a list of integration sites & the number of fragments detected at
# that site.
# This script takes that list of entries sorted by the 2nd column (chromosome, 
# presumably), and removes entries within X bp of entries with a higher
# value in the first column (count, presumably).  X is given on the 
# command line; in the main shell script, X=5.

# The way it goes about doing so is somewhat more complicated than the other
# scripts.  In essence, it collects all the entries on a chromosome, sorts them
# by integration count, and starts moving sites to a hash of things to be 
# printed as output.  If it finds something already in the hash that is within
# 5 bp (inclusive) of the site in question, then the site in question is not
# moved to the output hash.  In this way, preference is given to sites with
# more fragments associated with them.  Ties are broken by position on the 
# chromosome, with smaller-numbered positions winning out.  There's no 
# biological reason for breaking ties that way - rather, it's there so, if bias
# is introduced, we at least know what kind of bias there would be.


################################################################################

my $prev_chrom;
my %count_hash;
my %output_hash;
my $flag = 0;
my $same_count = 0;

unless ( @ARGV == 2 ) {
	die "Usage: remove_nearby_v3.0.pl input distance\n$!";
}

chomp (my $distance = $ARGV[1]);

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( open INPUT, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

unless ( open OUTPUT, ">", "$ARGV[0].not_in_$distance" ) {
	die "Cannot open file $ARGV[0].not_in_$distance: $!";
}


while (<INPUT>) {

	chomp;
	my @line = split;
	# Here, $prev_site is the very first site.  Later, it serves as the last 
	# site to have been printed.
	$prev_chrom = $line[1];
		
	# Make the position the key, and make the count the value.
	$count_hash{$line[2]} = $line[0];
	last;
	
}


while (<INPUT>) {

	chomp;
	my @line = split;
	my $curr_chrom = $line[1];
	
	# Check if the current read is on the same chromosome as the previous read
	if ( $curr_chrom eq $prev_chrom ) {
		
		# This is the easy one - if you're on the same chromosome, just
		# put the new entry into the hash.
		$count_hash{$line[2]} = $line[0];

	} else {
		
		# If you're NOT on the same chromosome, it's time to look through
		# the count_hash, empty it, and fill it up again.
		
		# This produces a list of positions, in descending order
		# of the count (note that the list does NOT contain count info)
		my @pos_array_sort = sort {
			$count_hash{$b} <=> $count_hash{$a} or
			$a <=> $b
		} keys %count_hash;
		
		# Check each position in the hash, starting with the one with the most
		# fragments aligned.	
		foreach my $position (@pos_array_sort) {
			
			my $lower_bound = $position - $distance;
			my $upper_bound = $position + $distance;
			
			# Look at all positions in the output hash
			foreach ( keys %output_hash ) {
				
				if ( ($_ >= $lower_bound) && ($_ <= $upper_bound) ) {
					
					# If there's something else within 5 bp of the current 
					# integration, set the flag to note that the current 
					# integration should NOT be printed.
					$flag = 1;
					
					if ($count_hash{$_} == $count_hash{$position}) {
						
						# If the aforementioned "something else" has the same
						# integration count as the current integration, keep
						# track of how often that comes up.
						$same_count++;
						
					}
					
				}
				
			}
			
			if ($flag) {
				$flag = 0;
				next;
			}
			
			# If you've looked through the whole output hash and nothing
			# else is nearby, you can add that entry to the output hash.
			$output_hash{$position} = 1;
			
			# Need to print count, chromosome, and position.
			print OUTPUT "$count_hash{$position}\t$prev_chrom\t$position\n";
			
		}
		
		# Once you've printed all the good entries for that chromosome,
		# empty all the lists and hashes and continue to the next.
		
		%count_hash = ();
		%output_hash = ();
		@pos_array_sort = ();
		
		# Now, put the current entry into the count hash, and reset the
		# chromosome.
		$count_hash{$line[2]} = $line[0];
		
		$prev_chrom = $curr_chrom;
		
	}
	
}

# Once the file has been read through, the last chromosome still has all
# its entries in the hash, not yet printed.  Print them as above.

my @pos_array_sort = sort {
	$count_hash{$b} <=> $count_hash{$a} or
	$a <=> $b
} keys %count_hash;

foreach my $position (@pos_array_sort) {
	
	my $lower_bound = $position - $distance;
	my $upper_bound = $position + $distance;
	
	# Look at all positions in the output hash
	foreach ( keys %output_hash ) {
		
		if ( ($_ >= $lower_bound) && ($_ <= $upper_bound) ) {
			$flag = 1;
			if ($count_hash{$_} == $count_hash{$position}) {
				$same_count++;
			}
					
		}
		
	}
	
	if ($flag) {
		$flag = 0;
		next;
	}

	
	# If you've looked through the whole output hash and nothing
	# else is nearby, you can add that entry to the output hash.
	$output_hash{$position} = 1;
	
	# Need to print count, chromosome, and position.
	print OUTPUT "$count_hash{$position}\t$prev_chrom\t$position\n";
	
}

close OUTPUT;

close INPUT;

print "There were $same_count reads removed with the same count as their neighbor.\n";

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;