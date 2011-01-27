## remove duplicates

last = ""

while line=ARGF.gets do 
  if line.match("^#")
    puts line
  else
    cols = line.split(/\s+/)
    chr,pos=cols[0], cols[1]
    current = chr + ":" + pos

    if current != last
      puts line
    end
    
    last = current
  end
end


