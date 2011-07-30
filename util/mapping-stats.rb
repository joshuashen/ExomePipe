#!/usr/bin/env ruby

### use samtools to read the mapping stats, then calculate the number of mapped reads

bam = ARGV[0]

a = `samtools idxstats #{bam}`.split(/\n/)

mapped = 0
unmapped = 0

a.each do |line|
  # example: 2   243199373   13531738   170957

  b = line.chomp.split(/\s+/)
  mapped += b[2].to_i
  unmapped += b[3].to_i
end

ratio = mapped.to_f/(mapped+unmapped)

puts "#{bam}\t#{mapped}\t#{unmapped}\t#{ratio.round(2)}\t#{mapped + unmapped}"



