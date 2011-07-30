## extract gene and functionClass info from vcf

def main
  while line=ARGF.gets do
    cols = line.chomp.split(/\t/)

    snp, info = cols[2], cols[7]
    gene, func = "", {}
    info.split(';').each do |f|
      k,v = f.split("=")
      if k.match("refseq.functionalClass")
        func[v] = 1
      elsif k.match("refseq.name2")
        gene = v
      end

    end
    puts "#{snp}\t#{gene}\t#{func.keys.join(";")}"
  end
end

main()
