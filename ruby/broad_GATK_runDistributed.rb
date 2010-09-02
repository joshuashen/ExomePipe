#!/util/bin/ruby


# some global variables

# replace with your GATK command-line here
$gatkCommandFront = "java -jar dist/GenomeAnalysisTK.jar -T GATKPaperGenotyper -R /broad/1KG/reference/human_b36_both.fasta -nt 4 -I /broad/hptmp/aaron/bams/22sorted.bam"
$chromSize = 247249719

# the command we use to parallelize processing
$bsub = "bsub -q gsa"


# validate the input
if (ARGV[0].to_i < 1 || ARGV[0].to_i > $chromSize)
  puts "#{ARGV[0].to_i} is not a valid split size"
  exit(1)
end

if (ARGV.size != 2)
  puts "two required arguments: the split size, and the output directory location (minus the trailing slash please)"
  exit(1)
end

$splitSize = ARGV[0].to_i
outputDir = ARGV[1]

if (!File.directory?(outputDir))
  puts "#{outputDir} appears not to be a directory"
  exit(1)
end

# some variables - start, expected output, and files to merge
start = 1
expectedFiles = []
mergeFiles = []

# start spliting the file
while (start < $chromSize) 
  stop = start + $splitSize
  if (stop > $chromSize)
    stop = $chromSize
  end
  # put together the job command
  command = "#{$bsub} -o #{outputDir}/#{$splitSize}.#{start}.out #{$gatkCommandFront} -cl #{outputDir}/#{$splitSize}.#{start}.calls.out -L 22:#{start}-#{stop}"
  # add it to the expected file list
  expectedFiles.push("#{outputDir}/#{$splitSize}.#{start}.out")
  mergeFiles.push("#{outputDir}/#{$splitSize}.#{start}.calls.out")
  # dispatch the job
  # puts "#{command}"
  output = `#{command}`
  puts output
  start = start + $splitSize
end

runtimes = []

while (expectedFiles.size > 0)
  expectedFiles.each { |ex|
    # puts "#{ex}"
    if (File.file?(ex))
      expectedFiles.delete(ex)
      fl = File.open(ex,"r")
      found = false
      fl.each { |line|
        # puts "#{line}"
        if (line =~ /\s+CPU\s+time\s+\:\s+(\d+\.\d+)/)
          runtimes.push($1.to_f)
          found = true
        end
      }
      if (!found)
        puts "couldn't find CPU time in output for file #{ex}"
        exit(1)
      end
      File.delete(ex);
    end
  }
  sleep(2)
end

# cat the files together, given a list of files in order
def catFiles(fileList)
  files = "cat"
  fileList.each { |file|
    if (File.exists?(file))
      files = files + " #{file}"
    end
  }
  files = files + " > #{$splitSize}results.from.cat"
  t1 = Time.now
  output = `#{files}`
  t2 = Time.now
  puts "output result was : #{output}"
  fileList.each { |file|
    if (File.exists?(file))
      File.delete(file)
    end
  }
  return (t2-t1)
end

# cat the files, and get the time taken
catTiming = catFiles(mergeFiles)

min = 10000000000
max = 0
average = 0
runtimes.each { |runtime|
  if (runtime > max)
    max = runtime
  end
  if (runtime < min)
    min = runtime
  end
  average = average + runtime
  puts runtime
}
puts "#{min},#{max},#{average/runtimes.size}"
puts "cat time = #{catTiming}"
