## get polyphen score for VT test

snpl = ARGV[0]
pscore = ARGV[1]

snps = {}
File.new(pscore, 'r').each do |line|
# example: 
#  #o_snp_id                     pph2_prob
#  chr10:100013464.CT.uc001kpa.1     0.994
#  chr10:100015474.GA.uc001kpa.1     0.022
#  chr10:100017453.TG.uc001kpa.1         1

  next if line.match("^#")
  cols=line.chomp.split(/\s+/)
  info=cols[0].split('.')[0]
  snps[info] = cols[1] 
  
end


File.new(snpl, 'r').each do |line|
  
  # example :
  #       "#{gene}\t#{num}\t0.5\t#{snv[:chr]}\t#{pos}\t#{snv[:fclass]}"  # polyphen score is set to 0.5                 
  cols = line.chomp.split(/\s+/)
  
  info = "chr#{cols[3]}:#{cols[4]}"
  if snps.key?(info)
    p = snps[info]
  else
    p = 0.5
  end
  if cols[5].to_i > 1 ## nonsense or readthrough
    p = 1
  end

  puts "#{cols[0]}\t#{cols[1]}\t#{p}\t#{cols[3..-1].join("\t")}"
end

  
 
