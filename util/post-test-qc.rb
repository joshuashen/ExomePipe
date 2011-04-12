## 

snplist = ARGV[0]
hardy = ARGV[1]
missing = ARGV[2]

snpa = []
snph = {}
gt = {}
hwe  = {}
miss = {}
File.new(snplist, 'r').each do |line|
  # example: 
  # ACACB 1 0 12 109577344 1 12_109577344
  # ACACB 2 0 12 109617728 1 rs16940029

  cols=line.chomp.split(/\s+/)
  if cols[-1] =~ /^rs/
    s = cols[-1]
  else
    s = "var_" + cols[-1]
  end
  snpa << s
  snph[s] = line.chomp
  gt[s] = []
  hwe[s] = []
  miss[s] = ""
end

File.new(hardy, 'r').each do |line|
  # example:   
#  12               rs2300455      ALL    A    G           21/131/327   0.2735   0.2959       0.1206
#  12               rs2300455      AFF    A    G              6/28/73   0.2617    0.304       0.1954
#  12               rs2300455    UNAFF    A    G           15/103/254   0.2769   0.2936        0.288
  cols = line.strip.split(/\s+/)
  s = cols[1]
  next unless snph.key?(s)
  
  if cols[2] == "AFF" or cols[2] == "UNAFF"
    gt[s] << cols[5]
    hwe[s] << cols[8]
  end
end

File.new(missing, 'r').each do |line|
  # example:    1              rs75062661        0.271       0.2258       0.3659

  cols = line.strip.split(/\s+/)
  s = cols[1]
  if snph.key?(s)
    miss[s] = cols[2..-1].join("\t")
  end
end

# header:
    
snpa.each do |s|
  puts "#{snph[s]}\t#{gt[s].join("\t")}\t#{miss[s]}\t#{hwe[s].join("\t")}"
end
