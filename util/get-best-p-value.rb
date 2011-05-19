## 

while line=ARGF.gets do 
  cols = line.split(/\t/)
  gene=cols[0]
  pa= []
  cols[1..-1].each do |p|
    pa << p.to_f unless p.match("NA")
  end
  bestp = pa.sort[0]
#  puts pa.join("\t")
  puts "#{gene}\t#{bestp}"
end

