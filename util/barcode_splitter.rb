## split fastq files based on barcode,
## barcode matching using hamming distance ( <=2 to real and >2 to the other)
##  alternatively, we can generate a set of possible codes with hd <=2 and then use hash

#require "narray"

#class String   ## hamming distance and xor function from Kirk Haines
#  def xor(other)
#    if other.empty?
#      self
#    else
    #  left = self
#      right = other

#      if left.length < right.length
#        n,r = right.length.divmod(left.length)
#        left = left * n + left[0, r]
#      elsif right.length < left.length
#        n,r = left.length.divmod(right.length)
#        right = right * n + right[0, r]
#      end

#    (NArray.to_na(self, "byte") ^ NArray.to_na(other, "byte")).to_s
#    end
#  end
  
#  def hamming_distance(other)
#    self.xor(other).tr("\x00",'').length
#  end
# end

require 'zlib'

def main
  barcode = ARGV[0]
  barcodefastq = ARGV[1]
  targetfastq = ARGV[2]

  barcodesize = 6

  coding = readBar(barcode)
  multiplex = {}
  outputio = {}
  coding.each do |str,sampleID|
    mutate1(str).each do |strmut|   # only allow hamming distance of 1
      multiplex[strmut] = sampleID
    end
    outputio[sampleID] = File.new(targetfastq + "_" + sampleID + ".fastq",'w')
  end
  $stderr.puts "allowed multiplexing codes: #{multiplex.keys.size}"
  $stderr.puts "multiplex mapping: #{multiplex}"

  outputio["discarded"] = File.new(targetfastq + "_discarded.fastq", "w" )
  assignment = decode(barcodefastq, targetfastq, multiplex, outputio, barcodesize)

  outputio.each  do |sampleID, fio|
    fio.close
  end
end

def decode(barcodefq, targetfq, multiplex, outio, barcodesize)
  if barcodefq.match(/.gz$/) # gzipped
    bio = Zlib::GzipReader.new(File.open(barcodefq))
  else
    bio = File.new(barcodefq, 'r')
  end
  
  if targetfq.match(/.gz$/)
    tio = Zlib::GzipReader.new(File.open(targetfq))
  else
    tio = File.new(targetfq, 'r')
  end
#   Zlib::GzipReader.new
 
  chunk = 1000000  # 1 million lines each time
  while !bio.eof? # not end
    i = 0
    j = 0
    decoded = []

    bio.lines.take(chunk).each do |str|
      bcode = str[0..barcodesize-1]
      i += 1
      j = (i - 2) / 4
      if i % 2 == 0 and i % 4 != 0 ## barcode line
        if multiplex.key?(bcode)  # good code
          decoded[j] = multiplex[bcode]
        else 
          decoded[j] = "discarded"
        end
      end
    end
#    puts decoded.size
#    puts decoded.join("\t")
    k = 0
###    i = 0
    # bigstring = tio.lines.take(chunk)
    for k in 0..(decoded.size - 1)
      str = tio.lines.take(4) 
     # puts str
      outio[decoded[k]].puts str
    end
  end
  bio.close
  tio.close
end

def mutate2(string)
  ## mutation the string within hammingDist of 2
  hd1 = mutate1(string)
  hd2=[]
  hd1.each do |s|
    hd2.concat(mutate1(s))
  end
  return hd1.concat(hd2)
end

def mutate1(string)
  l = string.length
  hd1 = []
  0.upto(l-1) do |i|
    cp = String.new(string)
    ["A", "T", "G", "C", "."].each do |nt|
      cp[i] = nt
      hd1 << String.new(cp)
    end
  end
  return hd1
end


def readBar(b)
  # example:
  # # ID   barcode
  # CDH088  CGATGT
  # CDH089  TGACCA
  coding = {}
  File.new(b, 'r').each do |line|
    next if line.match(/^#/) # header line
    cols = line.chomp.split(/\s+/)
    sampleID, code = cols[0], cols[1]
    coding[code] = sampleID
  end
  return coding
end

main()


## NOte: optimal file reading:
# File.open(ARGV.first, 'r') do |fi|
#  fsize = fi.size
#  fi.read(fsize).lines { |l| 
#  }
# end
