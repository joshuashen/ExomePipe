## atomic mapping wrapper for Illumina/Solexa reads
#  paired or single-end reads
#  batch size:  single run, <100M reads.

### author: Yufeng Shen, c2b2, Columbia


require 'getoptlong'

def main

  opts = GetoptLong.new(
     ["--reference", "-r", GetoptLong::OPTIONAL_ARGUMENT],
     ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
     ["--pair", "-p", GetoptLong::OPTIONAL_ARGUMENT],
     ["--help", "-h", GetoptLong::NO_ARGUMENT],
     ["--bwa", "-b", GetoptLong::OPTIONAL_ARGUMENT]
  )
  
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  
  if optHash.key?("--help")
    $stderr.puts "Usage: ruby __.rb -i foo.1.fastq [-p foo.2.fastq] [-r reference.fasta]"
    $stderr.puts "Note: mapping paired reads requires -p argument."
    exit
  end

  if optHash.key?("--reference") 
    ref = optHash["--reference"]
  else
    ref = "/ifs/data/c2b2/ip_lab/shares/DATA/Sequencing/resources/human_g1k_v37.fasta"
  end

  f1 = optHash["--input"]
  if optHash.key?("--pair")
    f2 = optHash["--pair"]
  else
    f2 = nil
  end
 
  if optHash.key?("--bwa")
    bwa = File.expand_path(optHash["--bwa"])
  else
    bwa = "/ifs/data/c2b2/ip_lab/shares/SOFTWARE/bwa_titan/bwa"
  end

## todo: optionally provide samtools location through command line.
  samtools = "/ifs/home/c2b2/ip_lab/yshen/usr/bin/samtools"

## todo: provide more options of bwa for fine-tuning. 
  threads = 2 ## number of threads used by bwa
  options = "" # other bwa options


  ######### align step
  cmd="#{bwa} aln -t #{threads} #{options}  #{ref} #{f1} > #{f1}.sai" 
  
  system(cmd)

  if f2 != nil # paired 
    cmd="#{bwa} aln -t #{threads} #{options}  #{ref} #{f2} > #{f2}.sai" 
    system(cmd)

    #### map
    output = File.dirname(File.expand_path(f1)) + "/" +   getPrefix(File.basename(f1)) + ".aligned.bam"
    cmd = "#{bwa} sampe #{ref} #{f1}.sai #{f2}.sai #{f1} #{f2} | #{samtools} view -bS -o #{output}  - "
    system(cmd)
    ## sort and index the BAM file
    cmd= "#{samtools} sort #{output} #{output}.sorted"
    system(cmd)
    cmd = "#{samtools} index #{output}.sorted"
    system(cmd)
  else # single reads
    output = f1 + ".bam"
    cmd = "#{bwa} samse #{ref} #{f1}.sai #{f1} | #{samtools} view -bS -o #{output} - "
    system(cmd)
# We need to make a function for this.
    cmd= "#{samtools} sort #{output} #{output}.sorted.bam"
    system(cmd)
    cmd = "#{samtools} index #{output}.sorted.bam"
    system(cmd)
  end
end

def sortAndIndex(bam)
#

end

def getPrefix(fname)
   return fname
   # $stderr.puts fname
   # return fname.match(/^(.*)[\.|\-]1[\.|\-](.*)$/)[1]
end

main()

