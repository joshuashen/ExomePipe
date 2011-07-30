## randomly take N reads from fastq file

$VERBOSE = nil

def main
  fq = ARGV[0]
  num = ARGV[1]
  
  chunk = 2000000  # 2 million reads
  
  if fq == nil or num == nil
    puts "Usage: ruby _.rb fastq NUM > random_NUM.fastq"
    exit
  end

  num = num.to_i
  total = `wc -l #{fq} | cut -f1 -d ' '`
  totalN = total.to_i / 4
  
  io = File.new(fq, 'r')
  array = []
  while !io.eof?
    s = io.take(4)
    array << s

    if array.size > chunk
      randomTake(array, num, totalN)
      array = []
    end
  end
  randomTake(array, num, totalN)
  io.close
  
end

def randomTake(array, num, t)
  n = num * array.size / t
  return if  n < 1
  puts (array.sort_by {rand})[0..n-1]
  

end

main()
