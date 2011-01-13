## screen vcf files by only taking positions within a targeted list

## assume both vcf and target files are sorted

vcf=ARGV[0]
targets=ARGV[1]

intervals = {}

File.new(targets,'r').each do |line|
  cols=line.chomp.split(/\s+/)
  chr,s,e=cols[0],cols[1].to_i, cols[2].to_i
  if !intervals.key?(chr)
    intervals[chr] = []
  end
  intervals[chr] << [s,e]
end

# $stderr.puts '@@ intervals loaded'

File.new(vcf, "r").each do |line|
  next if line.match("^#")
  if line=~/^(\S+)\s+(\d+)/
    chr,p=$1, $2.to_i
    

    while (intervals[chr].size > 0 &&  p > intervals[chr][0][1]) do
      intervals[chr].shift
    end
    
    if intervals[chr].size < 1 or  p < intervals[chr][0][0] # smaller than top -> out of range, discard
      $stderr.puts line
    else
      puts line
    end
  end
end



