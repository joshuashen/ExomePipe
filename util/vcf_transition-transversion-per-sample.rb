## calcualte transition/transversion ratio from VCF file
## per sample
## also count the number of synonymous, missense, silent var per sample

def main
  vcf = ARGV[0]
  codingOnly = ARGV[1]

  codingSwitch = 1
  if codingOnly != nil and codingOnly.to_i <= 0 ## 
    codingSwitch = -1
  end

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
#  while line=ARGF.gets do 
  File.new(vcf, 'r').each do |line|
    next if line.match(/^\##/)
    cols=line.chomp.split(/\s+/)
    if line.match(/^#CHROM/)  # header
      cols[9..-1].each do |cc|
        cc = cc.sub(" ","_")
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
    else
      pos,name,ref,alt, passflag, info, gt = cols[1].to_i,cols[2],cols[3],cols[4], cols[6], cols[7], cols[9..-1]
      next if ref.size != alt.size ## indel   || pos == lastpos ## indels or same var
      next if passflag != "PASS"
      ## SYN AND NONSYN

      synonFlag, missenseFlag, nonsenseFlag, knowntiFlag, knowntvFlag ,noveltiFlag, noveltvFlag, readthroughFlag = 0,0,0, 0, 0, 0,0, 0
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
      
      next if pos == lastpos       

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
        if genotype == '0/1'  ## het
          het[subject] += 1
          
        elsif genotype == '1/1' ## homo
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
      lastpos = pos
    end
#    puts "#knownti/tv: #{knownti}/#{knowntv} = #{knownti/knowntv.to_f}"
#    puts "#novelti/tv: #{novelti}/#{noveltv} = #{novelti/noveltv.to_f}"
  end

  puts "\#items\t#{sid.join("\t")}"
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

main()

