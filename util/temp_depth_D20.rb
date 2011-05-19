## 
depthStat = ARGV[0]

samples = {}
d15 = {}
d10 = {}
d20 = {}
File.new(depthStat, 'r').each do |line|
  cols = line.chomp.split(/\s+/)
  sample, depth = cols[0], cols[2]
  samples[sample] = depth
  d15[sample] = cols[6].to_f
end

samples.keys.each do |sample|
  file = "#{sample}.cleaned.bam_pipe/all.recalibrated.bam.coverage.sample_interval_statistics"
  x = `cut -f12,17,22 #{file} | tail -1`
  cols = x.chomp.split(/\s+/)
  d10[sample] = cols[0].to_f / cols[1].to_f * d15[sample]
  d20[sample] = cols[2].to_f / cols[1].to_f * d15[sample]
  puts "#{sample}\t#{samples[sample]}\t#{d10[sample]}\t#{d15[sample]}\t#{d20[sample]}"
end

  
