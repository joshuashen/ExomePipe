#!/usr/bin/env ruby

require 'getoptlong'


### filter VCF based on:
# 1. calling rate (AN / 2 / N_samples)
# 2. ratio of zero-mapping-qual reads ( MQ0 / DP) or RMS of MQ 
# 3. allele frequency (AF or AC / 2 / N_samples)
#  etc

#default values
def main 
  settings = {}
  settings["--missing"]=0.8
  settings["--mq"]=20
  settings["--zmq"]=0.4  # ((MQ0 / (1.0 * DP)) > $zmq 
  settings["--freq"] = -1  # default no filter
  settings["--ac"] = 0  
  settings["--maxfreq"]=1 # default no filter
  # ab=0.95  # allele balance
  settings["--qual"]=50.0 # min qual
#  settings["--clusterWinSize"]=10 #   --clusterWindowSize = 10
  settings["--HRun"]=5 # homopolymer
  settings["--qd"]=5.0 # qual over depth cutoff 
  settings["--sb"]=-0.1 # strand bias
  settings["--minDP"]=5 # average DP per sample, need to multiply by nsample
  
  optHash = getopt()
  vcf = optHash["--vcf"]
  
  settings.keys.sort.each do |s|
    if optHash.key?(s)
      settings[s] = optHash[s].to_f
    end
  end
  
#  if optHash.key?("--indelMask")
#    settings["--indelMask"] = optHash["--indelMask"]
#  end
  
  nsample=countSamples(vcf)
  
  filterVCF(vcf,settings,nsample)  # gt: gene -> pos -> sample -> genotype, 

end

def filterVCF(vcf, settings, nsample)
  o = File.new(vcf + ".filtered.vcf", 'w')

  File.new(vcf, 'r').each do |line|
    if line.match("^#") 
      o.puts line
    else
      cols=line.chomp.split(/\s+/)
      qual, info = cols[5].to_f, cols[7].split(';') 
      dp = 0     
      flag = 0 
      if qual < settings["--qual"]
        flag = 1
      else
        info.each do |item|
          k,v=item.split("=")[0..1]
#          $stderr.puts "#{k}\t#{v}"
          if k == "DP" 
            dp = v.to_f
            if dp < settings["--minDP"] * nsample
              flag = 1
#              $stderr.puts "#{cols[0]}\t#{cols[1]}\t : dp is too small"
            end
          elsif k == "QD" and v.to_f < settings["--qd"]
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t :  QD too small"
          elsif k == "HRun" and v.to_f > settings["--HRun"]
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t : Hrun "
          elsif k == "SB" and  v.to_f > settings["--sb"]
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t :  SB "
          elsif k == "AN" and v.to_f < nsample * 2.0 * ( 1 - settings["--missing"] )
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t :  AN "
          elsif k == "AF" and ( v.to_f < settings["--freq"] or v.to_f > settings["--maxfreq"])
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t : AF " 
          elsif k == "AC" and v.to_f < settings["--ac"]
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t :  AC  "
          elsif k == "MQ" and v.to_f < settings["--mq"]
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t :  MQ "
          elsif k == "MQ0" and v.to_f > settings["--zmq"] * dp
            flag = 1
#            $stderr.puts "#{cols[0]}\t#{cols[1]}\t : MQ0 "
          end
        end
      end

      if flag == 0
        o.puts line
      else
        o.puts "#{cols[0..5].join("\t")}\tStandardFilters\t#{cols[7..-1].join("\t")}"
      end
    end
  end
  o.close
end

def countSamples(vcf)
  n = 0
  File.new(vcf, 'r').each do |line|
    n += 1
    if line.match("^#CHROM") 
      cols = line.chomp.split(/\s+/)
      return cols.size - 9 
    elsif n > 1000
      return 0
    end
  end
end

def getopt
  
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        # ["--indelMask", "-i", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--missing", "-m", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--ac", "-c", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--mq", "-q", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--zmq", "-z", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--freq", "-f", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--maxfreq", "-F", GetoptLong::OPTIONAL_ARGUMENT],
#                        ["--clusterWinSize", "-w", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--qd", "-Q", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") 
    $stderr.puts "Usage: ruby __.rb -v VCF [options]"
    $stderr.puts "     options: "
    exit
  end
  return optHash
  
end


main()

exit
