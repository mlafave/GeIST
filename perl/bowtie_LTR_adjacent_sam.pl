#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 5/2/12


# Manipulates a bowtie output file, such that the positions reported correspond
# to the leftmost base of the sequence in which the virus integrated.  For MLV,
# this is the leftmost base of the 4 bp used for (and duplicated during) 
# integration.

################################################################################


unless ( @ARGV == 4 ) {
	die "Usage: bowtie_LTR_adjacent_v4.0.pl LTR_F_trimmed LTR_R_trimmed bowtie repeat_length\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

unless ( -e $ARGV[2] ) {
	die "$ARGV[2] doesn't exist: $!";
}

# On the command line, define the length of the integration site duplicated
# sequence.  For example, the value is '4' for MLV; the value would be '5' for
# HIV.
my $repeat = $ARGV[3];

my %LTR_F_hash;
my %LTR_R_hash;
my %position_hash;
my %strand_hash;


unless ( open LTRF, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

# Make a hash that lists LTR_F names; do the same for LTR_R.  Use files that
# are post-trim, so there won't be any duplicates between the two files.

# First, get the names of reads with forward-facing LTR.
while (<LTRF>) {
	
	chomp;
	
	# if it looks like a name, put it in the hash and skip to the next place a
	# name should be.  Otherwise, check the next line.
	if (m%^(@(.*)/(1|2)(\|p)?)$%) {
		
		$LTR_F_hash{$_} = 1;
		<LTRF>;
		<LTRF>;
		<LTRF>;
		
	} else {
		
		next;
	
	}
	
}

close LTRF;

unless ( open LTRR, "<", "$ARGV[1]" ) {
	die "Cannot open file $ARGV[1]: $!";
}

# Do the same for reads with reverse LTR.
while (<LTRR>) {
	
	chomp;
	
	# if it looks like a name, put it in the hash and skip to the next place a
	# name should be.  Otherwise, check the next line.
	if (m%^(@(.*)/(1|2)(\|p)?)$%) {
		
		$LTR_R_hash{$_} = 1;
		<LTRR>;
		<LTRR>;
		<LTRR>;
		
	} else {
		
		next;
	
	}
	
}

close LTRR;

# With two hashes of names, go through the bowtie file, using orientation and
# strand to decide how to manipulate the alignment value.

unless ( open BOWTIE, "<", "$ARGV[2]" ) {
	die "Cannot open file $ARGV[2]: $!";
}

while (<BOWTIE>) {
	
	chomp;
	my @bowtie_line = split;
	
	my $name = '@' . $bowtie_line[0];
	my $strand = $bowtie_line[1];
	# $bowtie_line[3] is the position, and $bowtie_line[9] is the read sequence.
	
	# Figure out if the read is LTR_F or LTR_R
	my $orientation;
	if ( exists $LTR_F_hash{$name}) {
		
		$orientation = 'F';
		delete $LTR_F_hash{$name};
		
	} elsif ( exists $LTR_R_hash{$name}) {
		
		$orientation = 'R';
		delete $LTR_R_hash{$name};
		
	} else {
		
		# If the name isn't in either hash, which shouldn't happen, it's not
		# something I want.
		next;
		
	}
	
	
	# Determine what side the adjacent base is on, and make the change, if 
	# needed.
	if ( $orientation eq 'F') {
		
		if ( $strand == 0 ) {
			
			# Base you want is 5' base, which is the leftmost.
			# No change is needed.
			
			# Integration was on the + strand.
			$strand_hash{$bowtie_line[0]} = "+";
			
		} elsif ( $strand == 16 ) {
			
			# Want the 5' base, which here is the rightmost.
			$bowtie_line[3] += length ($bowtie_line[9]) - $repeat;
			
			# Integration was on the - strand.
			$strand_hash{$bowtie_line[0]} = "-";
			
		} else {
			
			next;
			
		}
		
	} elsif ( $orientation eq 'R') {
		
		if ( $strand == 0 ) {
			
			# Want the 3' base, which here is the rightmost.
			$bowtie_line[3] += length ($bowtie_line[9]) - $repeat;
			
			# Integration was on the - strand.
			$strand_hash{$bowtie_line[0]} = "-";
			
		} elsif ( $strand == 16 ) {
			
			# Want the 3' base, which here is the leftmost.
			# No change is needed.
			
			# Integration was on the + strand.
			$strand_hash{$bowtie_line[0]} = "+";
			
		} else {
			
			next;
			
		}
		
	} else {
		
		# This should never happen (having neither + nor - orientation), but 
		# just in case:
		next;
	
	}
	
	# Now that the position is accurate, add it to the position hash.
	# The key is the name (without an '@'), and the value is the new position.
	
	$position_hash{$bowtie_line[0]} = $bowtie_line[3];
		
}

close BOWTIE;

%LTR_F_hash = ();
%LTR_R_hash = ();

# Use the bowtie file as a list of names to remove fragments with inconsistent
# positions between reads.
unless ( open BOWTIE, "<", "$ARGV[2]" ) {
	die "Cannot open file $ARGV[2]: $!";
}

while (<BOWTIE>) {
	
	chomp;
	my @bowtie_line = split;
	
	my $name = $bowtie_line[0];
	my $pair_name = $name;
	
	# Generate the name of the paired read
	if ( substr($pair_name, -1) eq '1' ) {
		substr($pair_name, -1) = '2';
	} elsif ( substr($pair_name, -1) eq '2' ) {
		substr($pair_name, -1) = '1';
	} else {
		next;
	}
	
	if ( exists $position_hash{$name} && exists $position_hash{$pair_name} ) {
		
		# If both halves are still in the hash, make sure they have the same
		# position.
		if ( $position_hash{$name} == $position_hash{$pair_name} ) {
			
			# If they are the same, great!  Picking the base directly adjacent
			# to the LTR should have the same position in both the + and the - 
			# reads.  Carry on.
			next;
			
		} else {
			
			# If the positions differ, remove both.  You've already determined
			# that both exist, so no extra 'if' test is needed.
			delete $position_hash{$name};
			delete $position_hash{$pair_name};
			next;
			
		}
		
	} else {
		
		# If one or both of the fragments are not in the hash, it means either
		# 1) it was a long fragment, or 2) the positions have already been found
		# to be inconsistent.  Either way, there's nothing more to delete, so 
		# move on.
		next; 
		
	}
	
}

close BOWTIE;




unless ( open BOWTIE, "<", "$ARGV[2]" ) {
	die "Cannot open file $ARGV[2]: $!";
}

unless ( open OUTPUT, ">", "$ARGV[2]_LTRadjacent" ) {
	die "Cannot open file $ARGV[2]_LTRadjacent: $!";
}	

while (<BOWTIE>) {
	
	chomp;
	my @bowtie_line = split;
	
	my $name = $bowtie_line[0];
	
	if ( exists $position_hash{$name}) {
		
		# If the entry is still in the hash at this point, it's safe to print.
		# Any bad reads would have already been removed.
		
		print OUTPUT "$bowtie_line[0]\t$bowtie_line[1]\t$bowtie_line[2]\t$position_hash{$name}\t$bowtie_line[4]\t$bowtie_line[5]\t$bowtie_line[6]\t$bowtie_line[7]\t$bowtie_line[8]\t$bowtie_line[9]\t$bowtie_line[10]\t$bowtie_line[11]\t$bowtie_line[12]\t$bowtie_line[13]\t$strand_hash{$name}\n";
				
	} else {
		
		# If the name isn't in the hash, that's because the read positions
		# were not consistent.  Carry on.
		next;
		
	}
	
}


close OUTPUT;

close BOWTIE;


################################################################################


my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
