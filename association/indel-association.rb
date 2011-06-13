## looking for genes with more than one non-synonymous changes for each sample


def main
  vcf=ARGV[0]
  pheno = ARGV[1]

  puts ["chr", "pos", "gene", "qual", "c1", "c0", "n1", "n0", "g1", "g0","hweChisq1","hweChisq0", "function", "rare"].join("\t")
  
  freqCutoff = 0.1 # variants with freq higher are deemed as common var.
  distCutoff = 10  # ignore variants within 10bp window
  
  pheno = readpheno(pheno)
  
  n1, n0 = 0, 0
  pheno.each do |s,p|
    if p == "case"
      n1 += 1
    elsif p == "control"
      n0 += 1
    end
  end
#  n1 = pheno.values.map {|v| v == "case"}.size
#  n0 = pheno.values.map {|v| v == "control"}.size

  sid=[]
  gt = []
  gene = {}
  sample = {}

  lastpos = -100


  File.new(vcf, 'r').each do |line|
    cols = line.chomp.split(/\t/)
   
    if line.match("^#CHROM") # header
      cols[9..-1].each do |i|
    
        i = i.sub(" ","")
        sample[i] = {}
        sid << i
      end

    elsif line.match("^#")
      next
    else  # normal line
      chr,pos,id,a1, a2, qual,filter,info, gt = cols[0], cols[1].to_i, cols[2], cols[3], cols[4], cols[5].to_f, cols[6], cols[7], cols[9..-1]
      flag = 0 
      rare = 0
      functionclass = "unknown"
      geneName = "unknown"
      
      # discard insertions and SNPs
      next if a1.size <= a2.size


      info.split(';').each do  |item|
        k,v = item.split('=')[0..1]
        
        if k == "AF" # allele freqency
          if v.to_f < freqCutoff  ## 
            rare = 1
          else
            rare = 0
          end
        elsif k =~ /refseq.changesAA/
          if v == "true"  # non-syn
            flag = 1
          end
        elsif k =~ /refseq.functionalClass/
          functionclass = v

        elsif k =~ /refseq.name2/
          geneName = v
        end
      end
#      $stderr.puts functionclass
      #      if rare == 1 and flag == 1 # non-syn and rare
      alleles = {"case" => 0, "control" => 0}
      genotypes = {"case" => [0,0,0], "control" => [0,0,0]}

      0.upto(gt.size-1) do |i|
        gtc = gt[i].split(':')
        sname = sid[i]
        next unless pheno.key?(sname)
#        if gtc[0] == '0/1'  ##  only consider hets or gtc[0] == '1/1' # carrier
        if gtc[0] == '0/1'
          
          alleles[pheno[sname]] += 1
          genotypes[pheno[sname]][1] += 1
        elsif gtc[0] == "1/1"
          alleles[pheno[sname]] += 2
          genotypes[pheno[sname]][2] += 1
        else 
          genotypes[pheno[sname]][0] += 1
        end
      end
      
      chisq1 = chisq(genotypes["case"])
      chisq0 = chisq(genotypes["control"])
      puts "#{chr}\t#{pos}\t#{geneName}\t#{qual}\t#{alleles["case"]}\t#{alleles["control"]}\t#{n1}\t#{n0}\t#{genotypes["case"].join('/')}\t#{genotypes["control"].join('/')}\t#{chisq1}\t#{chisq0}\t#{functionclass}\t#{rare}"
    end
  end



end

def chisq(arr)
  total = (arr[0] + arr[1] + arr[2]) 
  p = (arr[1] + arr[2] * 2 ) / 2.0 /  total
  x = 0
  pvector = [ (1 - p) * (1 - p), 2 * p * (1 - p), p * p ]
  0.upto(2)  do |i|
    obs = arr[i]
    expected = total * pvector[i]
    x += (obs - expected) ** 2 / expected
  end
  return x
end


def readpheno(fam)
  s = {}
#  s[:cases] = {}
#  s[:controls] = {}
  File.new(fam,'r').each do |line|
    cols = line.chomp.split(/\s+/)
    if cols[-1] == "1" # control
      s[cols[0]] = "control"
#      s[:controls][cols[0]] = 1
    elsif cols[-1] == "2" # case
#      s[:cases][cols[0]]= 1
      s[cols[0]] = "case"
    end
  end
  return s
end


main()

