## fgrep -f INPUT foo.txt | ruby __.rb INPUT

def main
  list = ARGV.shift

  keyHash = {}
  keyArray = []
  File.new(list,'r').each do |line|
    k = line.chomp
    #keyHash[k] = ""
    keyArray << k
  end

  while line=ARGF.gets do 
    cols = line.chomp.split(/\s+/)
    k,v = cols[0],cols[1]
    
    
    keyHash[k] = v
    
  end

  keyArray.each do |k|
    if keyHash.key?(k)
      puts "#{k}\t#{keyHash[k]}"
    else
      puts "#{k}\tNA"
    end
  end
end

main()
