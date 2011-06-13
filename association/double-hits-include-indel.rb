## looking for genes with more than one non-synonymous changes for each sample


def main
  vcf=ARGV[0]
  
  freqCutoff = 0.1 # variants with freq higher are deemed as common var.
  distCutoff = 10  # ignore variants within 10bp window

  sid=[]
  gt = []
  gene = {}
  sample = {}

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
      functionclass = ""
      geneName = ""
      
      # discard insertions
      next if a1.size < a2.size


      info.split(';').each do  |item|
        k,v = item.split('=')[0..1]
        
        if k == "AF" # allele freqency
          if v.to_f < freqCutoff  ## 
            rare = 1
          else
            break
          end
        elsif k =~ /refseq.changesAA/
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
#      $stderr.puts functionclass
      if rare == 1 and flag == 1 # non-syn and rare
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
        flag = 0
        lastpos = -2 * distCutoff 
        
        numVar = 0
        sample[i][geneName].sort.each do |pos|
          if pos - lastpos >= distCutoff ## not too close
            numVar += 1
          end
          lastpos = pos
        end

        if numVar > 1 
          sample[i][geneName].sort.each do |pos|
            fclass = gene[geneName][pos]["class"]
            #   str = str + "#{pos}/#{fclass};"
            if fclass == "nonsense" or fclass == "readthrough" or fclass == "frameshift" # 
              flag = 1
             end 
          end
          
          if flag == 1
            n += 1
            str += "\t#{i}:"
            sample[i][geneName].each do |pos|
              fclass = gene[geneName][pos]["class"]
              str = str + "#{pos}/#{fclass};"
            end
          end
        end
      end
    end
    
    puts "#{geneName}\t#{n}#{str}"
  end
end


main()

