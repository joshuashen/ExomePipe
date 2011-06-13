## remove duplicates

last = ""

while line=ARGF.gets do 
  if line.match("^#")
    puts line
  else
    cols = line.split(/\s+/)
    chr,pos,a1,a2=cols[0], cols[1], cols[3], cols[4]
    current = [chr, pos, a1, a2].join(":")

    if current != last
      puts line
    end
    
    last = current
  end
end


