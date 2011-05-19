## case control counting

def main
  dhits = ARGV[0]
  pheno = ARGV[1]

  samples = readpheno(pheno)
  
  parse(dhits, samples)
end

def parse(dhits, samples)
  tc1, tc0 = 0, 0
  n1 = samples[:cases].size
  n0 = samples[:controls].size
  
  File.new(dhits, 'r').each do |line|
    cols = line.chomp.split(/\s+/)
    ns = cols[1].to_i
    next if ns == 0
    
    c1, c0 = 0, 0
    cols[2..-1].each do |info|
      sname = info.split(':')[0]
      if samples[:cases].key?(sname)  # a case
        c1 += 1
      else
        c0 += 1
      end
    end

    puts "#{cols[0]}\t#{c1}\t#{c0}\t#{n1}\t#{n0}"
    
    tc1 += c1
    tc0 += c0
    
  end

  $stderr.puts "#{tc1}\t#{tc0}\t#{n1}\t#{n0}"
end
def readpheno(fam)
  s = {}
  s[:cases] = {}
  s[:controls] = {}
  File.new(fam,'r').each do |line|
    cols = line.chomp.split(/\s+/)
    if cols[-1] == "1" # control
      s[:controls][cols[0]] = 1
    elsif cols[-1] == "2" # case
      s[:cases][cols[0]]= 1
    end
  end
  return s
end

main()
