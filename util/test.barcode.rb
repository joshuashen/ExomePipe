require "narray"

class String   ## hamming distance and xor function from Kirk Haines                                                   \
                                                                                                                        
  def xor(other)
    if other.empty?
      self
    else
      left = self
      right = other

      if left.length < right.length
	n,r = right.length.divmod(left.length)
	left = left * n + left[0, r]
      elsif right.length < left.length
	n,r = left.length.divmod(right.length)
	right = right * n + right[0, r]
      end

      left_na = NArray.to_na(left, "byte")
      right_na = NArray.to_na(right, "byte")

      (left_na ^ right_na).to_s
    end
  end

  def hamming_distance(other)
    self.xor(other).tr("\x00",'').length
  end
end


barcodes = ARGV[0]

r1 = "TGACCA"
r2 = "CGATGT"

File.new(barcodes, 'r').each do |line|
  c = line.strip.split(/\s+/)
  n, code = c[0].to_i, c[1]
  hd1 =	r1.hamming_distance(code[0..5])
  hd2 = r2.hamming_distance(code[0..5])
  puts "#{code}\t#{n}\t#{hd1}\t#{hd2}"
end

