# Author: Yufeng "Joshua" Shen
# map Solexa (Illumina) short reads using bwa

require 'getoptlong'

$bwa = ''
def main 
  # logic of the program
  bwa = tryBWA()
  if (bwa == false)
    $stderr.puts "Error: couldn't find bwa. "
    exit
  else
    $maq = File.expand_path(maq)
  end

  optHash = getOptions() ## "--inputDir"
    
  if optHash.key?("--help") 
    help()
    exit
  end

# maxMismatch, maxDist
  if optHash.key?("--maxMismatch")
    maxMismatch = optHash["--maxMismatch"].to_i
  else
    maxMismatch = 3
  end

  if optHash.key?("--maxDist")
    maxDist = optHash["--maxDist"].to_i
  else
    maxDist = 300
  end

  if optHash.key?("--batchSize")
    batchSize = optHash["--batchSize"].to_i * 4
  else
    batchSize = 2000000  # default 500,000 reads per batch
  end

  absInputPath = File.expand_path(optHash["--inputDir"])
  absOutputPath = File.dirname(absInputPath) + "/Mapping_of_" + File.basename(absInputPath)
  if !File.exist?(absOutputPath)
    system("mkdir #{absOutputPath}")
  end
  $stderr.puts "convert Solexa export files into fastq files \nand split into smaller batches..."


  if optHash.key?("--do") 
    if optHash["--do"] == "1" # write qsub scripts
      flag = 1
    elsif optHash["--do"] == "2"  # write a bsub script
      flag = 2
    else
      flag = 0  # do in a loop locally
    end


  else
    ## if no do, just do maq
    flag = 0
  end

  if !optHash.key?("--nosplit")
    export2fastqAndDivide(absInputPath, absOutputPath, batchSize)
  end
  
  ref = File.expand_path(optHash["--ref"])
  if !optHash["--ref"].match("bfa$") 
    system("#{$maq} fasta2bfa #{ref} #{ref}.bfa")
    ref = ref + '.bfa'
  end

  log = 'maq_log.txt'
  pairedEndMapping(absOutputPath, ref, flag, log, maxMismatch, maxDist)
  
end

def getOptions 
  opts = GetoptLong.new(
      ["--inputDir", "-i", GetoptLong::REQUIRED_ARGUMENT],
      ["--ref", "-r", GetoptLong::REQUIRED_ARGUMENT],
      ["--maxMismatch", "-n", GetoptLong::OPTIONAL_ARGUMENT],
      ["--maxDist", "-a", GetoptLong::OPTIONAL_ARGUMENT],
      ["--batchSize", "-b", GetoptLong::OPTIONAL_ARGUMENT],
      ["--do", "-d", GetoptLong::OPTIONAL_ARGUMENT],
      ["--nosplit", "-z", GetoptLong::NO_ARGUMENT],                        
      ["--help", "-h", GetoptLong::NO_ARGUMENT]
  )

  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  return optHash
end

def help 
  $stderr.puts "Usage: ruby __.rb -i InputDir -r Reference.bfa [-b batchSize] [ --do 0] [-n maxMismatch] [-a maxDist] [-z]"
  $stderr.puts " --do    default 0: "
  $stderr.puts "         0  execute Maq on local machine one by one in a loop"
  $stderr.puts "         1  write qsub shell scripts for mapping through Maq"
  $stderr.puts "         2  write a bsub shell script for mapping through Maq"
  $stderr.puts "\n -b     batchSize, default 1000000 reads per batch"
  $stderr.puts "\n -z     just do maq mapping, assuming the split step has been completed"
end
# see if maq is avaible
def tryMaq
  trymaq = `which maq`
  if trymaq == '' # maq not found
    $stderr.puts "Failed to find Maq. Please make sure Maq is installed and you have put the path to Maq in $PATH of your shell environment"
    return false
  else
    return trymaq.chomp
  end
end

def export2fastqAndDivide(input, output, batchSize)
# get the file names
  files = Dir.entries(input).sort.select {|file| file.match("export.txt$") }


  files.each do |file|
    $stderr.puts "work on #{file}.."
    absFile = input + "/" + file
    absOutPrefix = output + '/' + file + "_"
    cmd = %q[awk '{print "@"$1":"$2":"$3":"$4":"$5"/"$6"\n"$7"\n+\n"$8}'] +  " #{absFile} " + %q[ | split -l ] + " #{batchSize} - #{absOutPrefix}" 
    system(cmd)
    batches = Dir.entries(output).sort.select {|batch| batch.match(file) and !batch.match("bfq") and !batch.match("fastq")}
    batches.each do |batch|
      batchabs = output + '/' + batch
      fastq = batchabs + ".fastq"
      bfq = batchabs + ".bfq"
      solexa2MaqFastq(batchabs, fastq)
      fastq2bfq(fastq, bfq)
      File.delete(batchabs)
    end
  end
end


def solexa2MaqFastq(sf, mf)
  system("#{$maq} sol2sanger #{sf} #{mf}")
end

def fastq2bfq(f,b)
  system("#{$maq} fastq2bfq #{f} #{b} 2> /dev/null")
end

def pairedEndMapping(dir, ref, flag, log, maxMismatch, maxDist)
  # find pairs of bfq files, do maq
  bfqsFreads = Dir.entries(dir).sort.select {|file| file.match("bfq$") and file.match('_1_export')}
  
  bfqsFreads.each do |fReads|
    rReads = fReads.sub("_1_export", "_2_export")
    r1 = dir + '/' + fReads
    r2 = dir + '/' + rReads
    errorLog = r1 + '.stderr_log'
    stdoutlog = r1 + '.stdout_log'
    mapresult = r1 + '.map'
    maqstring = ''
    if File.exist?(r2)
      maqstring = "#{$maq} map  -n #{maxMismatch} -a #{maxDist}  #{mapresult}  #{ref} #{r1} #{r2} 2 >> #{log}"
      $stderr.puts "do maq: #{maqstring}"
    else
      $stderr.puts "Warning: cannot find the mate pair of #{fReads} :  #{rReads}"
    end
    
    if flag == 0 ## do locally
      $stderr.puts maqstring
      system(maqstring)
    elsif flag == 1 ## qsub
      shstring = '#!/bin/sh' + "\n" + '#$ -cwd' + "\n" + maqstring

      shf = File.new(r1+".sh", "w")
      shf.puts shstring
      shf.close
      system("qsub #{r1}.sh -e #{errorLog} -o #{stdoutlog}")
      
    elsif flag == 2 # bsub
      system("bsub -o #{stdoutlog} -e #{errorLog} #{maqstring}")
    end
  end
  return
end

main()    # do
