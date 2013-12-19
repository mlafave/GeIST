#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 9/27/13 (modified from remove_nearby_v3.0.pl)

# Input is a list of integration sites & the number of fragments detected at
# that site.
# This script takes that list of entries sorted by the 2nd column (chromosome, 
# presumably), and removes entries within X bp of entries with a higher
# value in the first column (count, presumably).  X is given on the 
# command line; in the main shell script, X=5.

# The way it goes about doing so is somewhat more complicated than the other
# scripts.  In essence, it collects all the entries on a chromosome, sorts them
# by integration count, and starts moving sites to a hash of things to be 
# printed as output.  It creates a "position blacklist" in a hash, such that 
# the base that serves to represent the integration, and each position within 
# 5 bp in both directions, are added to the hash before an entry is printed.  If
# another entry falls within those bases, it won't be printed.
# In this way, preference is given to sites with
# more fragments associated with them.  Ties are broken by position on the 
# chromosome, with smaller-numbered positions winning out.  There's no 
# biological reason for breaking ties that way - rather, it's there so, if bias
# is introduced, we at least know what kind of bias there would be.


################################################################################

my $prev_chrom;
my %count_hash;
my %blacklist_hash;
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
			
			if ( exists $blacklist_hash{$position}) {
				
				if ($blacklist_hash{$position} == $count_hash{$position}) {
					
					# If the other site has the same
					# integration count as the current integration, keep
					# track of how often that comes up.
					$same_count++;
										
				}
				
				# That position is disallowed because there was already an 
				# integration near it or in it.  Try the next position.
				
				next;

			}
			
			# If the current position is not in the blacklist hash, good.  It 
			# can be printed.  Before that, though, add it (and the nearby 
			# bases) to the blacklist hash.
			
			for ( my $i = ($position - $distance); $i <= ($position + $distance) ; $i++) {
				
				# If it's not already in the hash, add it with the value equal
				# to the fragment count.
				
				if ( ! exists $blacklist_hash{$i}) {
					$blacklist_hash{$i} = $count_hash{$position};
				}
				
			}
			
			# Need to print count, chromosome, and position.
			print OUTPUT "$count_hash{$position}\t$prev_chrom\t$position\n";
			
		}
		
		# Once you've printed all the good entries for that chromosome,
		# empty all the lists and hashes and continue to the next.
		
		%count_hash = ();
		%blacklist_hash = ();
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
	
	if ( exists $blacklist_hash{$position}) {
		
		if ($blacklist_hash{$position} == $count_hash{$position}) {
			
			# If the other site has the same
			# integration count as the current integration, keep
			# track of how often that comes up.
			$same_count++;
			
		}
		
		# That position is disallowed because there was already an 
		# integration near it or in it.  Try the next position.
		
		next;

	}
	
	# If the current position is not in the blacklist hash, good.  It 
	# can be printed.  Before that, though, add it (and the nearby 
	# bases) to the blacklist hash.
	
	for ( my $i = ($position - $distance); $i <= ($position + $distance) ; $i++) {
		
		# If it's not already in the hash, add it with the value equal
		# to the fragment count.
		
		if ( ! exists $blacklist_hash{$i}) {
			$blacklist_hash{$i} = $count_hash{$position};
		}
		
	}
	
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
