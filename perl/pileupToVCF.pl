#!/usr/bin/perl
#written by Rajan Banerjee
#reb2143



$ArgsLen = @ARGV;  #variable for keeping track of the number of files.
$output = $ArgsLen-1;
$tempfiles = 0;


$samConvert = "perl /ifs/home/c2b2/ip_lab/yshen/usr/bin/sam2vcf.pl";  #location of converter
$vcfUpdate = "perl /ifs/home/c2b2/ip_lab/reb2143/VCFtools/perl/vcf-convert";		#location of a version updater for VCF files
$vcfMerge = "perl /ifs/home/c2b2/ip_lab/reb2143/VCFtools/perl/merge-vcf";			 #location of a VCF file merger
$vcfValidate = "perl /ifs/home/c2b2/ip_lab/reb2143/VCFtools/perl/vcf-validator";   #location of a VCF file validator

if($ArgsLen < 2)
	{ die "please enter at least one input pileup file and exactly one VCF output file.\n
		Use as: perl ~./pileupToVCF snps1.samtools snps2.samtools out.vcf\n";}
		 #forces users to use proper format

if($ARGV[$output]!~/\.vcf$/){
	die "the last file entered should be your final .vcf format destination file\n
		Use as: perl ~./pileupToVCF snps1.samtools snps2.samtools out.vcf\n";}

system("$samConvert < $ARGV[0] > /ifs/home/c2b2/ip_lab/reb2143/temp/pileupToVCFTEMPFILEDONOTDUPLICATE.vcf");  #converts the file

system("cat /ifs/home/c2b2/ip_lab/reb2143/temp/pileupToVCFTEMPFILEDONOTDUPLICATE.vcf | $vcfUpdate -r ref.fa > $ARGV[$output]"); #updates the file to newest VCF version

system("rm /ifs/home/c2b2/ip_lab/reb2143/temp/pileupToVCFTEMPFILEDONOTDUPLICATE.vcf"); #removes temporary files created.


=head
if($ArgsLen>2) {
	for($i=1; $i<$output; $i++){
	system("$samConvert < $ARGV[$i] > pileupToVCFTEMPFILEDONOTDUPLICATE.vcf");
	
	system("$vcfMerge 	
	}
}
=cut
