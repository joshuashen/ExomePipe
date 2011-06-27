#!/usr/bin/env ruby
require 'getoptlong'

## extract a set of samples


def main
  optHash = getopt()

  samples = readSamples(optHash["--samples"])

  extract(optHash["--vcf"], samples) 

end

def readSamples(pheno)
  samples = {}
  if File.file?(pheno) # a file lists subjectNames
    File.new(pheno, 'r').each do |line|
      cols = line.chomp.split(/\s+/)
      sid = cols[0]
      samples[sid] = 1
    end
  else  # the string is the subjectName
    samples[pheno.chomp] = 1
  end
  
  return samples 
end

def extract(vcf, samples)
  sid = []
  sindex = []

  File.new(vcf, 'r').each do |line|
    puts line if line.match("^##")
    
    cols = line.chomp.split(/\t/) 

    if line.match("^#CHROM")  # header of VCF
      cols[9..-1].each do |cc|
        cc = cc.sub(" ","_")
        
        sid << cc
        #        sid=cols[9..-1]
      end
      
      0.upto(sid.size - 1) do |i|
        if samples.key?(sid[i]) ## sample of interest
          sindex << i
        end
      end    
      puts "#{cols[1..8].join}"
    else
    end  

  end
end


def getopt

  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--samples", "-s", GetoptLong::REQUIRED_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") or !optHash.key?("--samples")
    $stderr.puts "Usage: ruby __.rb -v foo.vcf -s samplelist"
    exit
  end
  return optHash
  
end

main()
