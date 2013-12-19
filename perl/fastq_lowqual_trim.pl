#!/usr/bin/env perl

use 5.8.8;
use strict;
use warnings;

# Author: Matt LaFave
# Written: 9/27/13

# This script takes a fastq file from a pipe and trims off the rightmost reads
# that have low quality (indicated by #).  It trims from the point after which
# there are only # quality calls.



while (<>) {
	
	chomp (my $name = $_);
	chomp (my $seq = <>);
	chomp (my $plus = <>);
	chomp (my $qual = <>);
	
	if ($qual =~ m/^(.*?)#*$/) {
		
		my $length = length($1);
		$seq = substr($seq, 0, $length);
		$qual = substr($qual, 0, $length);
		
	}
	
	print "${name}\n${seq}\n${plus}\n${qual}\n";
	
}



exit;
