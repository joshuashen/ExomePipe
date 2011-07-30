## extract indels from vcf file


def main
  vcf = ARGV[0]
  snvo = File.new(vcf + ".snv", 'w')
  indelo = File.new(vcf + ".indel", 'w')
  
  File.new(vcf,'r').each do |line|
    if line.match(/^#/) # comments
      snvo.puts line
      indelo.puts line
    else
      cols = line.split(/\s+/)
      ref,alt = cols[3], cols[4]
      if ref.size != alt.size ## indels
        indelo.puts line
      else
        snvo.puts line
      end
    end
  end
  
  snvo.close
  indelo.close
end

main()
