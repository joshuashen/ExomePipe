#!/usr/bin/env ruby

def main
  input = ARGV[0]
  barcodeSize = ARGV[1]

  fqout = File.new(input + "_Forward.fastq", 'w')
  bout = File.new(input + "_barcode.fastq", 'w')

  if barcodeSize != nil
    bsize = barcodeSize.to_i
  else
    bsize = 6 ## default size
  end


  io = File.new(input, 'r')
  while !io.eof? 
    unit = io.take(4)
    bc = unit[1][0..bsize-1]
    seq = unit[1][bsize..-1]
    bcq = unit[3][0..bsize-1]
    qual = unit[3][bsize..-1]

    fqout.puts unit[0] + seq + unit[2] + qual
    
    bout.puts unit[0] + bc + "\n" + unit[2] + bcq 
  end

  io.close
  fqout.close
  bout.close

end


main()
