sname = ARGV[0]
rfile = ARGV[1]

flag = 0
parray = []
File.new(rfile,'r').each do |line|
  if line=~ /^p-values/
    flag = 1
    
  elsif flag == 1 # start
    cols = line.chomp.split(/\s+/)
    if cols[0] == "1"  # p-value line
      cols[1..-1].each do |p|
        parray << p
      end
    end
  end

end

if parray.size == 8
  puts "#{sname}\t#{parray.join("\t")}"
end
    
## header: 
#geneName  CMC_1%  CMC_1%+polyphen2 CMC_5%  CMC_5%+polyphen2  Madson-Browning  Madson-Browning+polyphen2  VT-test VT-test+polyphen2
