#!/usr/bin/perl
use 5.8.8;

while (<>) {
	
	chomp;
	$name = $_;
	chomp ($seq = <>);
	chomp ($plus = <>);
	chomp ($qual = <>);
	
	$seq = reverse $seq;
	$seq =~ tr/ACGTNacgtn/TGCANtgcan/;
	
	$qual = reverse $qual;
	
	print $name."\n".$seq."\n".$plus."\n".$qual."\n";
	
}

exit;
