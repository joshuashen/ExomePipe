order=ARGV[0]
ori = ARGV[1]
cnum = ARGV[2]

if order== nil
   $stderr.puts "Usage: ruby __.rb order old_file"
end

if cnum == nil
  cnum = 1 # second column
else
  cnum = cnum.to_i
end

output = "#{ori}_re-ordered"

chr = []
chrf = {}
chrh = {}
File.new(order, 'r').each do |line|
  c = line.chomp
  chr << c
  chrf[c] = "temp_#{c}.temp"
  chrh[c] = File.new(chrf[c],'w')
end

File.new(ori, 'r').each do |line|
  cols = line.split(/\s+/)
  if chrh.key?(cols[cnum])   # 
    chrh[cols[cnum]].puts line
  end
end

chr.each do |c|
  chrh[c].close
  `cat #{chrf[c]} >> #{output}`
  `rm -f #{chrf[c]}`
end


    
