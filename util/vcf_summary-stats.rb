## check Mendelian errors 

# input: one vcf file with family subjects

require 'getoptlong'

def main
  optHash = getopt()
  samples = readFam(optHash["--fam"])

  if optHash.key?("--switch")  # codingswitch: 1 mean only consider coding var 
    compare(optHash["--vcf"], samples, 1)
  else
    compare(optHash["--vcf"], samples, -1)
  end
end

def compare(vcf, samples, codingSwitch)
  pidx, midx, oidx = nil, nil, nil


  novelErr = 0
  knownErr = 0
  synonErr = 0
  missenseErr = 0
  nonsenseErr = 0
  readthroughErr = 0
  missing = 0 
  totalCodinghomo = 0
  
  
  novelti = {}
  noveltv = {}
  knownti = {}
  knowntv = {}
  synon = {}
  missense = {}
  nonsense = {}
  readthrough = {}
  homo = {}
  het  = {}
  sid = [] # sample array

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
        sid << cc
        novelti[cc] = 0
        noveltv[cc] = 0
        knownti[cc] = 0 
        knowntv[cc] = 0
        synon[cc] = 0
        missense[cc] = 0
        nonsense[cc] = 0
        readthrough[cc] = 0
        homo[cc] = 0
        het[cc] = 0
      end
    else  # var lines
      pos,name,ref,alt, passflag, info, gt = cols[1].to_i,cols[2],cols[3],cols[4], cols[6], cols[7], cols[9..-1]
      next if ref.size != alt.size 
      next if passflag != "PASS"
      next if pos == lastpos 
      
      synonFlag, missenseFlag, nonsenseFlag, knownFlag, novelFlag, readthroughFlag = 0,0,0, 0, 0, 0,0, 0
      
      knowntiFlag, knowntvFlag ,noveltiFlag, noveltvFlag = 0,0,0, 0

      err = 0 
      coding = -1
    
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
      
      next unless coding * codingSwitch > 0  # only take coding var
     
      
      if name =~ /^rs/  # known
        knownFlag = 1
      else
        novelFlag = 1
      end

      
      ti = judge(ref,alt)
      if name =~ /^rs/  # known
        if ti == 1
          knowntiFlag = 1
        else
          knowntvFlag = 1
        end
      else  # novel
        if ti == 1
          noveltiFlag = 1
        else
          noveltvFlag = 1
        end
      end
      
      
      
      #      $stderr.puts gt
      i = 0 
      gt.each do |genotypeinfo|
        subject = sid[i]
        genotype = genotypeinfo.split(":")[0]
        if genotype == '0/1' and coding > 0  ## het
          het[subject] += 1
          
        elsif genotype == '1/1' and coding > 0 ## homo
          homo[subject] += 1
        end
        
        if genotype == '0/1' or genotype == '1/1'
          novelti[subject] += noveltiFlag
          noveltv[subject] += noveltvFlag
          knownti[subject] += knowntiFlag
          knowntv[subject] += knowntvFlag
          synon[subject] += synonFlag
          missense[subject] += missenseFlag
          nonsense[subject] += nonsenseFlag
          readthrough[subject] += readthroughFlag
        end
        i += 1
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

      if ogenotype == "1/1" and coding > 0
        totalCodinghomo += 1 
      end
      lastpos = pos
    end
    
  end
  puts "Comparison b.t.w #{samples[:p]} (P) and #{samples[:o]} (O)"
  puts "Total Mendelian errors:\t#{knownErr + novelErr}"
  puts "Total Mendelian errors in coding:\t#{synonErr + missenseErr + nonsenseErr + readthroughErr}"
  puts "Non-call in one subject:\t#{missing}"
  puts "Mendelian errors in SNPs from dbSNP129:\t#{knownErr}"
  puts "Mendelian errors in novel SNVs:\t#{novelErr}"
  puts "Mendelian errors in synonymous SNVs:\t#{synonErr}"
  puts "Mendelian errors in missense SNVs:\t#{missenseErr}"
  puts "Mendelian errors in nonsense SNVs:\t#{nonsenseErr}"
  puts "Mendelian errors in readthrough SNVs:\t#{readthroughErr}"
  puts "Number of coding/homozygous-non_ref SNVs in proband:\t#{totalCodinghomo}"
  
  puts "\n\n\#items\t#{sid.join("\t")}"
  print "known:ti/tv"
  sid.each {|s| print "\t#{knownti[s]}/#{knowntv[s]}"}
  print "\n"

  print "known:ti/tv-ratio"
  sid.each {|s| print "\t#{(knownti[s]/knowntv[s].to_f).round(3)}"}
  print "\n"


  print "novel:ti/tv"
  sid.each {|s| print "\t#{novelti[s]}/#{noveltv[s]}"}
  print "\n"

  print "novel:ti/tv-ratio"
  sid.each {|s| print "\t#{(novelti[s]/noveltv[s].to_f).round(3)}"}
  print "\n"

  print "silent"
  sid.each {|s| print "\t#{synon[s]}" }
  print "\n"
  
  print "missense"
  sid.each {|s| print "\t#{missense[s]}" }
  print "\n"

  print "nonsense"
  sid.each {|s| print "\t#{nonsense[s]}" }
  print "\n"

  print "readthrough"
  sid.each {|s| print "\t#{readthrough[s]}" }
  print "\n"

  print "homozygous"
  sid.each {|s| print "\t#{homo[s]}" }
  print "\n"

  print "heterozygous"
  sid.each {|s| print "\t#{het[s]}" }
  print "\n"
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

def judge(a1, a2)
  ti = 0
  if a1 == 'A'
    if a2 == 'G'  # ti
      ti = 1
    end
  elsif a1 == 'C'
    if a2 == 'T'
      ti = 1
    end
  elsif a1 == 'T'
    if a2 == 'C'
      ti = 1
    end
  elsif a1 == 'G'
    if a2 == 'A'
      ti = 1
    end
  end
  return ti
  
end

def getopt
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--fam", "-f", GetoptLong::REQUIRED_ARGUMENT],
                        ["--switch", "-s", GetoptLong::NO_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") or !optHash.key?("--fam")
    $stderr.puts "Usage: ruby __.rb -v VCF_or_list -f fam [-s]"
    $stderr.puts " if -s is provided, only consider coding sites"
    exit
  end
  return optHash
  
end

main()

