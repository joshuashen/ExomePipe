#!/usr/bin/perl -w
use strict;
use Carp;

## From  Ryan Kim.

my $usage = q(
qseq2fastq.pl  - a script to convert all qseq files in a directory into a single fastq file with sanger-style ASCII q-score encoding
USAGE: qseq2fastq.pl <qseq.txt file> <output file>
);

if (scalar(@ARGV) != 2) {
    print $usage;
    exit;
}

my $in_file = $ARGV[0];
my $output_fastq_file = $ARGV[1];

my $qfilter = "";
open(OUTFASTAQFILE, "> $output_fastq_file");

open INFILE, "< $in_file" || die "Error: Couldn't open $in_file\n";
while(<INFILE>)
{
    chomp;
    my @this_line = split/\t/, $_;
    croak("Error: invalid column number in $in_file\n") unless(scalar(@this_line) == 11);
    #if($this_line[10] == 1) {
    #$qfilter = "Y";
    #} else {
    #	$qfilter = "N";
    #}
    # Convert quality scores
    my @quality_array = split(//, $this_line[9]);
    my $phred_quality_string = "";
    # convert each char to Phred quality score
    foreach my $this_char (@quality_array){
	my $phred_quality = ord($this_char) - 64; # convert illumina scaled phred char to phred quality score
	my $phred_char = chr($phred_quality + 33); # convert phred quality score into phred char (sanger style)
	$phred_quality_string = $phred_quality_string . $phred_char;
    }
    # replace "." gaps with N
    $this_line[8] =~ s/\./N/g;
    # output line
    print OUTFASTAQFILE "@" . $this_line[0] . ":" . $this_line[2] . ":" . $this_line[3] . ":" . $this_line[4] . ":" . $this_line[5] . "#" . $this_line[6] . "/" . $this_line[7] . "\n" .  #header line
	$this_line[8] . "\n" . # output sequence
	"+" . $this_line[0] . ":" . $this_line[2] . ":" . $this_line[3] . ":" . $this_line[4] . ":" . $this_line[5] . "#" . $this_line[6] . "/" . $this_line[7] . "\n" . # header line
	$phred_quality_string . "\n"; # output quality string
}
close(INFILE);
close(OUTFASTAQFILE);
exit;

