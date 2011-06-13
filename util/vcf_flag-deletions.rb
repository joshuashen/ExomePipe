## flag deletion or insertion

while line = ARGF.gets do
  cols = line.chomp.split(/\t/)
  if cols[1].size > cols[2].size  # deletion
    puts "DEL\t#{line}"
  elsif cols[1].size < cols[2].size
    puts "INS\t#{line}"
  else
    puts "SNV\t#{line}"
  end
end
