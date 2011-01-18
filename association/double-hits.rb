## looking for genes with more than one non-synonymous changes for each sample


def main
  vcf=ARGV[0]
  
  sid=[]
  gt = []
  gene = {}
  sample = {}

  File.new(vcf, 'r').each do |line|
    cols = line.chomp.split(/\s+/)
    if line.match("^#CHROM") # header
      sid=cols[9..-1]
      sid.each do |i|
        sample[i] = {}
      end

    elsif line.match("^#")
      next
    else  # normal line
      chr,pos,id,qual,filter,info, gt = cols[0], cols[1].to_i, cols[2], cols[5].to_f, cols[6], cols[7], cols[9..-1]
      flag = 0 
      functionclass = ""
      geneName = ""
      info.split(';').each do  |item|
        k,v = item.split('=')[0..1]
        if k =~ /refseq.changesAA/
          if v == "true"  # non-syn
            flag = 1
          end
        elsif k =~ /refseq.functionalClass/
          functionclass = v

        elsif k =~ /refseq.name2/
          geneName = v
          if !gene.key?(geneName) 
            gene[geneName] = {}
          end
          gene[geneName][pos]={}
          gene[geneName][pos]["class"] = functionclass
        end
      end
      $stderr.puts functionclass
      if flag == 1 # non-syn
        0.upto(gt.size-1) do |i|
          gtc = gt[i].split(':')
          sname = sid[i]
          if gtc[0] == '0/1'  ##  only consider hets or gtc[0] == '1/1' # carrier
            if !sample[sname].key?(geneName) 
              sample[sname][geneName] = []
            end
            sample[sid[i]][geneName] << pos
#            if gtc[0] == '1/1' 
#              sample[sid[i]][geneName] << pos
#            end
          end
        end
      end
    end
  end


  # look for double hits
  gene.keys.sort.each do |geneName|
    n = 0
    str = ""
    sid.each do |i|
      if sample[i].key?(geneName) 
        if sample[i][geneName].size > 1 
          n += 1
          str += "\t#{i}:"
          sample[i][geneName].each do |pos|
            fclass = gene[geneName][pos]["class"]
            str = str + "#{pos}/#{fclass};"
          end
        end
      end
    end
    
    puts "#{geneName}\t#{n}#{str}"
  end
end


main()

