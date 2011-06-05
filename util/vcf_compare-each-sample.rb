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
  compare(gt1, gt2, samples)
end

def compare(gt1, gt2, samples)
  
  samples.each_key do |sample|
    ### :n1, :n2, :n0 -> all variants
    ### :k1, :k2, :k0 -> known variants
    venn = {:n1=>0, :n2=>0, :n0=>0, :k1=>0,:k2=>0, :k0=>0, :diff => 0}
    samples[sample] = venn
    
  end
  
  gt1.each do |pos, snv|
    
    nonref1 = (snv[:gt].values.select {|gtc| gtc > 0}).size  # number of
    
    if nonref1 > 0 
      if gt2.key?(pos)  # also called in gt2
        snv2 = gt2[pos]
        
        snv[:gt].each do |sample, gtc|
          gtc2 = snv2[:gt][sample]
          if gtc > 0  and gtc2 > 0 
            samples[sample][:n0] += 1
            if snv[:dbsnp] == 1  # known SNPs
              samples[sample][:k0] += 1
            end 
            if gtc2 != gtc    
              ## het/hom status discordant
              samples[sample][:diff] += 1                    
            end 
          elsif gtc > 0 and gtc2 <= 0 # 
            samples[sample][:n1] += 1
            if snv[:dbsnp] == 1  # known SNPs
              samples[sample][:k1] += 1
            end 
          elsif gtc <= 0 and gtc2 > 0 
            samples[sample][:n2] += 1
            if snv[:dbsnp] == 1  # known SNPs
              samples[sample][:k2] += 1
            end 
          end
        end
        gt2.delete(pos)  # remove pos => snv2 from gt2
      else   ## not called in gt2
        snv[:gt].each do |sample, gtc|
          if gtc > 0 
            samples[sample][:n1] += 1
            if snv[:dbsnp] == 1  # known SNPs
              samples[sample][:k1] += 1
            end 
          end
        end
      end
    end
    gt1.delete(pos)
  end
  
  gt2.each do |pos, snv| ## all v2
    
    nonref = (snv[:gt].values.select {|gtc| gtc > 0}).size  # number of non-ref gt
    if nonref > 0
      snv[:gt].each do |sample, gtc|
        if gtc > 0 
          samples[sample][:n2] += 1
            if snv[:dbsnp] == 1  # known SNPs
              samples[sample][:k2] += 1
            end 
        end
      end
    end
    
    gt2.delete(pos)
  end
  
  # puts "#FID\tVenn1\tVenn0\tVenn2\Venn1k"
  samples.each do |sample, venn|
    puts "#{sample}\t#{venn[:n1]}\t#{venn[:n0]}\t#{venn[:n2]}\t#{venn[:k1]}\t#{venn[:k0]}\t#{venn[:k2]}\t#{venn[:diff]}"
  end  
  
end


def readVCF(vcf, samples)
  gt = {}  ## gt: sample -> pos -> snv
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
      next if filter!= "PASS"

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
   #   else  ## dbSNP known position
   #     next ## for now only deal with novel SNPs
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
