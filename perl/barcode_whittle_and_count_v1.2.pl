#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 4/5/12
# Last revised: 7/5/12

# This script adds barcode information to a file with bowtie-output-style 
# formatting.  The "whittle" refers to the fact that entries without associated
# barcodes won't make it through this script (though there shouldn't be any at
# this step, it doesn't hurt to be cautious), and the "count" refers to counting
# the number of fragments that have each barcode.

# The way it goes about its job is a little roundabout, but the extra steps
# help slim down memory usage.

# Note that this script requires that reads with inconsistent barcodes be 
# removed in order to give an accurate count.

################################################################################

# The bowtie file is intended to be the 'island' or 'island_sort' file - the one
# near the very end, when things have already been whittled down to only those
# inserts that are not with x bases of another.

unless ( @ARGV == 3 ) {
	die "Usage: barcode_whittle_and_count_v1.2.pl bowtie_file barcode_cat barcode_ref\n$!";
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

# Declare global variables
my %name_hash;
my %read_hash;
my %fragment_hash;
my %count_hash;
my ($name, $seq);



unless ( open BOWTIE, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

# First, make a hash of the names of all the barcodes you need to retrieve.
while (<BOWTIE>) {
	
	chomp;
	my @line = split;
	
	my $name = '@'.$line[0];
	my $pair_name = $name;
	
	if ( substr($pair_name, -1) eq '1' ) {
		substr($pair_name, -1) = '2';
	} elsif ( substr($pair_name, -1) eq '2' ) {
		substr($pair_name, -1) = '1';
	} else {
		next;
	}
	
	$name_hash{$name} = 1;
	$name_hash{$pair_name} = 1;
	
}

close BOWTIE;




unless ( open BARCODEIN, "<", "$ARGV[1]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

while (<BARCODEIN>) {
	
	chomp;
	
	# Check to see if you're currently on a 'name' line.  If so, read in the
	# rest of the entry; if not, move on until you reach a name line.
	if (m%^(@(.*)/(1|2))$%) {
		
		# Verify that the barcode you're about to grab is one that will be used
		# in the bowtie file.
		if ( exists $name_hash{$_}) {
			
			$name = $_;		
			chomp ($seq = <BARCODEIN>);
			<BARCODEIN>; # Skip the plus
			<BARCODEIN>; # Skip the quality line
			
			$read_hash{$name} = $seq;
			delete $name_hash{$_};
			
		} else {
			
			<BARCODEIN>; # Skip the sequence
			<BARCODEIN>; # Skip the plus
			<BARCODEIN>; # Skip the quality line
			
		}
	
	} else {
	
		next;
		
	} # if

} # while

close BARCODEIN;

# Empty the now-useless name hash.
%name_hash = ();


unless ( open BARCODEREF, "<", "$ARGV[2]" ) {
	die "Cannot open file $ARGV[2]: $!";
}

while (<BARCODEREF>) {
	
	chomp;
	$count_hash{$_} = 0;
	
}

close BARCODEREF;



unless ( open BOWTIE, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

unless ( open BARCODEOUT, ">", "$ARGV[0]_barcodes" ) {
	die "Cannot open file $ARGV[0]_barcodes.fastq: $!";
}

while (<BOWTIE>) {
		
	chomp;
	my @line = split;
		
	my $name = '@'.$line[0];
	
	my $pair_name = $name;
	
	if ( substr($pair_name, -1) eq '1' ) {
		substr($pair_name, -1) = '2';
	} elsif ( substr($pair_name, -1) eq '2' ) {
		substr($pair_name, -1) = '1';
	}
	
	
	if ( exists $read_hash{$name}) {
		
		# if the read_hash contains the barcode for that read, print that
		# to the output file.
		
		print BARCODEOUT "$_\t$read_hash{$name}\t(x)\n";
		
		# Next, keep a running count of how many different fragments have 
		# been detected.
		# Capture the name of the FRAGMENT in $shortname.  The fragment name is
		# the portion of the name shared by both reads (basically, everything
		# that comes before /1 or /2).		
		my $shortname = substr($name, 0, (length($name) - 2));
		if (! exists $fragment_hash{$shortname}) {
			
			# If the basename ISN'T in the hash, it hasn't been seen yet.
			# Increment the count hash and add the name to the hash.
			$count_hash{$read_hash{$name}}++;
			$fragment_hash{$shortname} = 1;
			
		} else {
			
			# If the basename IS in the hash, it has been seen before.
			# Rather than count the same fragment twice, go to the next 
			# iteration.
			next;
			
		} # if
		
	} elsif ( exists $read_hash{$pair_name}) {
		
		# If the read_hash didn't have the barcode for the read, but DID have
		# the barcode for the pair of the read (such as would be the case for
		# long fragments, or short fragments with quality problems on one read),
		# then print out the barcode with a note that it came from the paired
		# read.
		
		print BARCODEOUT "$_\t$read_hash{$pair_name}\tBarcode_from_pair\n";
				
		my $shortname = substr($pair_name, 0, (length($pair_name) - 2));
		if (! exists $fragment_hash{$shortname}) {
			
			# If the basename ISN'T in the hash, it hasn't been seen yet.
			# Increment the count hash and add the name to the hash.
			$count_hash{$read_hash{$pair_name}}++;
			$fragment_hash{$shortname} = 1;
			
		} else {
			
			# If the basename IS in the hash, it has been seen before.
			# Go to the next iteration.
			next;
			
		} # if
	} # if
	
	# If, for some reason, a fragment did not have a barcode in the hash, don't
	# print either read from it at all.
	
}

close BARCODEOUT;

close BARCODEIN;

# Empty the fragment hash
%fragment_hash = ();


unless ( open BCCOUNT, ">", "$ARGV[0]_barcodesum" ) {
	die "Cannot open file $ARGV[0]_barcodesum: $!";
}

while ( my ($key, $value) = each %count_hash ) {
	
	# Print a summary of how many times each barcode was detected
	print BCCOUNT "$key\t$value\n";

}

close BCCOUNT;


my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
