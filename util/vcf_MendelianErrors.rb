## check Mendelian errors 

# input: one vcf file with family subjects

require 'getoptlong'

def main
  optHash = getopt()
  samples = readFam(optHash["--fam"])

  if optHash.key?("--snp")
    compare(optHash["--vcf"], samples, 1)
  else
    compare(optHash["--vcf"], samples, 0)
  end
end

def compare(vcf, samples, snpOnly)
  pidx, midx, oidx = nil, nil, nil


  novelErr = 0
  knownErr = 0
  synonErr = 0
  missenseErr = 0
  nonsenseErr = 0
  readthroughErr = 0
  missing = 0 
  totalCodinghomo = 0 
  lastpos = -1
 
  File.new(vcf, 'r').each do |line|
    next if line.match(/^\##/)
    cols=line.chomp.split(/\s+/)



    if line.match(/^#CHROM/)  # header
      i = 0
      cols[9..-1].each do |cc|
        cc = cc.sub(" ","_")
        if samples[:p] == cc
          pidx = i
        elsif samples[:m] == cc
          midx = i
        elsif samples[:o] == cc
          oidx = i
        end
        i += 1
      end
    else  # var lines
      pos,name,ref,alt, passflag, info, gt = cols[1].to_i,cols[2],cols[3],cols[4], cols[6], cols[7], cols[9..-1]
      next if ref.size != alt.size 

# #      next if ref.size != alt.size ## indel   || pos == lastpos ## indels or same var
      next if passflag != "PASS"
      next if pos == lastpos 
      
      synonFlag, missenseFlag, nonsenseFlag, knownFlag, novelFlag, readthroughFlag = 0,0,0, 0, 0, 0,0, 0
      err = 0 
      coding = 0
    
      info.split(';').each do |l|  
        k,v = l.split('=')
        
        if  k == "refseq.functionalClass" and v == "silent"  ## syn
          synonFlag = 1
          coding = 1
        elsif  k == "refseq.functionalClass" and v == "missense"  # missense
          missenseFlag = 1
          coding = 1
        elsif k == "refseq.functionalClass" and v == "nonsense"
          nonsenseFlag = 1
          coding = 1
        elsif k == "refseq.functionalClass" and v == "readthrough"
          coding = 1
          readthroughFlag = 1
        end
      end
      
      if name =~ /^rs/  # known
        knownFlag = 1
      else
        novelFlag = 1
      end

      pgt = gt[pidx]
#      mgt = gt[pidx] 
      ogt = gt[oidx]

      pgenotype = pgt.split(":")[0]
      ogenotype = ogt.split(":")[0]
#      next if pgenotype == '0/1'  or ogenotype == "0/1"  #
      
      if pgenotype != ogenotype 
        if ogenotype == './.'  or pgenotype == './.' # 
          missing += 1
        elsif pgenotype == '0/1'  or ogenotype == "0/1"
          1
        elsif ogenotype == "1/1"   # only interested in the proband
          err  = 1
          $stderr.puts line  ## need to investigate the depth coverage of these samples
        end
      end

      
      novelErr += novelFlag * err
      knownErr +=  knownFlag * err
      synonErr += synonFlag * err
      missenseErr +=  missenseFlag * err
      nonsenseErr += nonsenseFlag * err
      readthroughErr += readthroughFlag * err

      if ogenotype == "1/1"
        totalCodinghomo += coding 
      end
      lastpos = pos
    end
    
  end
  puts "Comparison b.t.w samples[:p] (P) and samples[:o] (O)"
  puts "Total Mendelian errors:\t#{knownErr + novelErr}"
  puts "Total Mendelian errors in coding:\t#{synonErr + missenseErr + nonsenseErr + readthroughErr}"
  puts "Non-call in one subject:\t#{missing}"
  puts "Mendelian errors in SNPs from dbSNP129:\t#{knownErr}"
  puts "Mendelian errors in novel SNVs:\t#{novelErr}"
  puts "Mendelian errors in synonymous SNVs:\t#{synonErr}"
  puts "Mendelian errors in missense SNVs:\t#{missenseErr}"
  puts "Mendelian errors in nonsense SNVs:\t#{nonsenseErr}"
  puts "Mendelian errors in readthrough SNVs:\t#{readthroughErr}"
  puts "Number of coding and homozygous alternative in proband:\t#{totalCodinghomo}"
end

def readFam(fam)
  # for one-parent <-> offspring pair
  samples = {:p => nil , :m => nil , :o => nil}
  File.new(fam, "r").each do |line|
    cols = line.chomp.split(/\s+/)
    if cols[1] == "0" # maternal
      samples[:m] = cols[0]
    elsif cols[1] == "1"  # paternal
      samples[:p] = cols[0]
    elsif cols[1] == "2"  # offspring
      samples[:o] = cols[0]
    end
  end
  return samples
end

def getopt
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--fam", "-f", GetoptLong::REQUIRED_ARGUMENT],
                        ["--snp", "-s", GetoptLong::NO_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") or !optHash.key?("--fam")
    $stderr.puts "Usage: ruby __.rb -v VCF_or_list -f fam [-s]"
    $stderr.puts " if -s is provided, only consider SNV sites"
    exit
  end
  return optHash
  
end

main()
