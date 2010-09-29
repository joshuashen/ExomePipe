#!/usr/bin/perl
#written by Rajan Banerjee
#reb2143
#A script to remove lines in a VCF file that the GenomicAnnotator cannot handle.

if(@ARGV != 2)
	{ die "please enter at least one VCF input file and one VCF output file";} #forces users to use proper format


open (IN, "$ARGV[0]") or die "Couldn't properly clean the file";
open OUT, ">", "$ARGV[1]" or die "Couldn't properly clean the files";



while(<IN>){
my $line = $_;
if ($_=~/^[0-9]{1,2}\t/){
$line = "chr".$_;
}
if($_=~/^(X|Y|N|M)/){
$line = "";
}
if($_=~/^\S*\t\S*\t\S*\t\S*\t\S{3,}\t/){
$line = "";
}
if($_=~/^#/){
$line = $_
}
print OUT $line;
}


close(IN);
close(OUT);

