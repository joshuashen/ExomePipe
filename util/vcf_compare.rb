#!/usr/bin/env ruby
require 'getoptlong'


def main
  optHash = getopt()
  vcf1, vcf2 = optHash["--vcf1"], optHash["--vcf2"]
  
  pflag = 0
  samples = {}
  if optHash.key?("--position")
    pflag = 1 ## do not care about homo/het status for non-ref alleles
  end

  if optHash.key?("--samples") ## 
    samples = readSamples(optHash["--samples"])
  end

  gt1 = readVCF(vcf1, samples)  # gt:  pos => snv  ;  snv: sample => genotype, 
  $stderr.puts "# of variants: #{gt1.keys.size}"
  gt2 = readVCF(vcf2, samples)
#   printout(gt, outprefix)
  nsample = samples.keys.size
  compare(gt1, gt2, nsample)
end

def compare(gt1, gt2, nsamples)
  a1 = gt1.size  # all var
  a2 = gt2.size  # 
  
  n1 = 0    # novel var
  n2 = 0
  mismatch = {}
  
  v1, v2, v0 = 0,0,0 ## Venn diagram: v1: specific to gt1, v2: specific to gt2; v0: overlap
  
  gtdiscorarray = Array.new(nsamples+1,0)
  posdiscorarray = Array.new(nsamples+1,0)
  
  gt1.each do |pos, snv|
    if snv[:dbsnp] == 0  ## novel
      n1 += 1
    else  
      next  ## only do novel SNV for now
    end 
    
    gtconcor, gtdiscor, posconcor,posdiscor = 0,0,0, 0
    
    nonref1 = (snv[:gt].values.select {|gtc| gtc > 0}).size  # number of
    
    if gt2.key?(pos)  # also called in gt2
      snv2 = gt2[pos]
      
      snv[:gt].each do |sample, gtc|
        gtc2 = snv2[:gt][sample]
        if gtc > 0  or gtc2 > 0 # non-ref
          v0 += 1
          if gtc2 == gtc 
            gtconcor += 1
            posconcor += 1
          elsif gtc > 0 and gtc2 > 0 ## 
            posconcor += 1
            gtdiscor += 1
          else
            gtdiscor += 1
            posdiscor += 1
          end                    
        
        end
      end
      gt2.delete(pos)  # remove pos => snv2 from gt2
    elsif nonref1 > 0  ## not called in gt2
      v1 += 1
      nonref = (snv[:gt].values.select {|gtc| gtc > 0}).size  # number of non-ref gt
      gtdiscor += nonref
      posdiscor += nonref
    end
    gtdiscorarray[0] += gtconcor
    gtdiscorarray[gtdiscor] += 1
    posdiscorarray[0] += posconcor
    posdiscorarray[posdiscor] += 1

    gt1.delete(pos)
  end
  
  gt2.each do |pos, snv| ## all v2
    
    nonref = (snv[:gt].values.select {|gtc| gtc > 0}).size  # number of non-ref gt
    if nonref > 0
      v2 += 1
    end 
 #   gtdiscor += nonref
 #   posdiscor += nonref
  end
  
  puts "Venn on positions of noval variants: #{v1} | #{v0} | #{v2}\n\n"
  puts "Histogram of discordant on the existance of non-reference alleles:"
  
  0.upto(posdiscorarray.size-1) do |i|
    puts "#{i}\t#{posdiscorarray[i]}"
  end

  puts "\n\nHistogram of discordant on genotypes:"
  
  0.upto(gtdiscorarray.size-1) do |i|
    puts "#{i}\t#{gtdiscorarray[i]}"
  end

#  puts "Histogram of discordant on the existance of non-reference alleles:"
  
#  posdiscorarray.each do |posdiscor|
#    puts "#{posdiscor+1}\t#{posdiscorarray[posdiscor]}"
#  end
  
  
end

def printout(gt, prefix)
  snpo = File.new(prefix + ".csnp", "w")
  gto = File.new(prefix + ".cgeno", "w")
  summary = File.new(prefix + ".table", "w")
  gt.keys.sort.each do |gene|
    num = 0
    gt[gene].keys.sort.each do |pos|
      num += 1
      snv = gt[gene][pos]
      snpo.puts "#{gene}\t#{num}\t0.5\t#{snv[:chr]}\t#{pos}\t#{snv[:fclass]}\t#{snv[:id]}"  # polyphen score is set to 0.5

      majorAllele = 0
      alleles = {0=> 0, 1=>0, 2=> 0, -9 => 0}
      gt[gene][pos][:controls].values.each do |a|
        alleles[a] += 1
      end
      
      if alleles[2] > alleles[0]
        majorAllele = 2
      end
      
      c1, c0 = 0,0
      n1 = snv[:cases].keys.size
      n0 = snv[:controls].keys.size
      snv[:cases].keys.sort.each do |sample|
        if  snv[:cases][sample] !=  majorAllele # carrier of a minor allele
          allele = (snv[:cases][sample] - majorAllele).to_i.abs
          gto.puts "#{gene}\t#{sample}\t#{num}\t#{allele}"
          c1 += 1
        end
      end
      snv[:controls].keys.sort.each do |sample|
        if  snv[:controls][sample] !=  majorAllele # carrier of a minor allele
          allele = (snv[:controls][sample] - majorAllele).to_i.abs
          gto.puts "#{gene}\t#{sample}\t#{num}\t#{allele}"
          c0 += 1
        end

      end
      
      summary.puts "#{gene}\t#{snv[:id]}\t#{c1}\t#{c0}\t#{n1}\t#{n0}"
    end
  end
  snpo.close
  gto.close
  summary.close
end

def readVCF(vcf, samples)
  gt = {}
  sid = []
  sindex = []

  File.new(vcf, 'r').each do |line|
    next if line.match("^##")
    
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
    else
      chr,pos,id,qual,filter,info, gtdetails = cols[0], cols[1].to_i, cols[2], cols[5].to_f, cols[6], cols[7], cols[9..-1]
      fclass = 0  # 0: synonymous; >0: non-syn
      dbsnp = 1
      
      info.split(';').each do  |item|
        k,v = item.split('=')[0..1]
        if k =~ /refseq.changesAA/
          if v == "true"  # non-syn
            fclass = 1
          end
        elsif k=~ /refseq.functionalClass/
          if v == "nonsense" or v == "readthrough" # set it as damaging by default
            fclass = 2
          end
        end
      end
      
      snv = {}

      if !id.match("^rs")
        id = "var_#{chr}_#{pos}"
        dbsnp = 0 # novel SNP
      else  ## dbSNP known position
        next ## for now only deal with novel SNPs
      end
      #snv[:pos] = pos
      snv[:fclass] = fclass
      snv[:dbsnp] = dbsnp
      snv[:gt] = {}
      gt[id]= snv
      
      sindex.each do  |i|  
     # 0.upto(gtdetails.size-1) do |i|
        gtc = gtdetails[i].split(':')
        sname = sid[i]
        snv[:gt][sname] = encodeGT(gtc[0])
      end
    end
  end
  
  return gt
end

def encodeGT(vcfcode)
  if vcfcode == '0/0'
    return 0
  elsif vcfcode == "0/1"
    return 1
  elsif vcfcode == "1/1"
    return 2
  else
    return -9
  end
end

def readSamples(pheno)
  samples = {}
  File.new(pheno, 'r').each do |line|
    cols = line.chomp.split(/\s+/)
    sid = cols[0]
    samples[sid] = 1
  end
  return samples 
end


def getopt

  opts = GetoptLong.new(
                        ["--vcf1", "-1", GetoptLong::REQUIRED_ARGUMENT],
                        ["--vcf2", "-2", GetoptLong::REQUIRED_ARGUMENT],
                        ["--position", "-p", GetoptLong::NO_ARGUMENT],
                        ["--samples", "-s", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf1") or !optHash.key?("--vcf2")
    $stderr.puts "Usage: ruby __.rb -1 1.vcf -2 2.vcf [-p] [-s samplelist]"
    exit
  end
  return optHash
  
end

main()
