#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 2/28/12

# Splits a FASTQ file into short reads (those with LTR detected on both reads)
# and long reads (everything else).

################################################################################

unless (@ARGV == 3) {
	die "Usage: split_LTR_to_long_and_short_v1.0.pl LTR_R LTR_F LTR_cat\n$!";
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

unless ( open LTR_R, '<', $ARGV[0]) {
	die "Can't open $ARGV[0]: $!";
}

my %LTR_R_name_hash;
my %short_hash;
my %long_hash;
my $count;
my $LTRcat_line;


# First, go through the file with reverse LTR, and put all the names in a hash.
while (<LTR_R>) {

	chomp;
	
	# If the line looks like a name, put it in the hash.  Otherwise, skip to the
	# next line.
	if (m%^(@(.*)/(1|2))$%) {
		
		$LTR_R_name_hash{$_} = 1;

		# Skip to what should be the next name line
		<LTR_R>;
		<LTR_R>;
		<LTR_R>;
		
	} else {
		
		next;
		
	}

}

close LTR_R;

unless ( open LTR_F, '<', $ARGV[1]) {
	die "Can't open $ARGV[1]: $!";
}

# Next, go through the forward LTR file.  If the line looks like a name, check 
# to see if the PAIR of that name is in the reverse LTR file.  If it is, the 
# fragment is probably short; put it in the short hash.  Otherwise, it's 
# probably long, so put it in the long hash.
while (<LTR_F>) {
	
	chomp;
	
	if (m%^(@(.*)/(1|2))$%) {
	
		my $pair_name = $_;
	
		if ( substr($pair_name, -1) == 1) {
			substr($pair_name, -1) = 2;
		} elsif (substr($pair_name, -1) == 2) {
			substr($pair_name, -1) = 1;
		} else {
			# If the name doesn't end with a 1 or a 2, it's weird.  
			# Go to the next entry.
			next;
		}
		
		# If the pair is in the reverse hash, put the original name in the 
		# appropriate hash.
		if ( exists ($LTR_R_name_hash{$pair_name})) {
		
			$short_hash{$_} = 1;
	
		} else {
			
			$long_hash{$_} = 1;
		
		}
		
		delete $LTR_R_name_hash{$pair_name};

	}

}

close LTR_F;

unless ( open LTR_CAT, '<', $ARGV[2]) {
	die "Can't open $ARGV[2]: $!";
}

unless ( open SHORTFRAG, '>', "$ARGV[1]_short.fastq") {
	die "Can't open $ARGV[1]_short.fastq: $!";
}

unless ( open LONGFRAG, '>', "$ARGV[1]_long.fastq") {
	die "Can't open $ARGV[1]_long.fastq: $!";
}



# Finally, go through the LTR_cat file (which is a concatenation of LTR F and R 
# files, but with duplicate entries removed).  If a name matches the short hash,
# print it to the appropriate file, etc. Many entries won't have forward LTR at 
# all; these entries should be skipped, as we already have the LTR_R file to
# deal with those.
while ( defined ($LTRcat_line = <LTR_CAT>)) {
	
	chomp $LTRcat_line;
	
	if ( exists ($short_hash{$LTRcat_line})) {
		
		chomp (my $seq = <LTR_CAT>);
		chomp (my $plus = <LTR_CAT>);
		chomp (my $qual = <LTR_CAT>);
		print SHORTFRAG "$LTRcat_line\n$seq\n$plus\n$qual\n";
		next;
		
	} elsif ( exists ($long_hash{$LTRcat_line})) {
		
		chomp (my $seq = <LTR_CAT>);
		chomp (my $plus = <LTR_CAT>);
		chomp (my $qual = <LTR_CAT>);
		print LONGFRAG "$LTRcat_line\n$seq\n$plus\n$qual\n";
		next;
		
	} else {
		
		next;
		
	}

}

close LONGFRAG;

close SHORTFRAG;

close LTR_CAT;


my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
