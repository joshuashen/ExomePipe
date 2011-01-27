#!/usr/bin/env ruby
require 'getoptlong'

def main
  perm1 = 10000
  perm2 = 5000000
  pcut = 0.05
  optHash = getopt()
  if optHash.key?("--perm") 
    perm2 = optHash["--perm"].to_i
  end
  #  sname = ARGV[0]
  # rfile = ARGV[1]
  
  parray=runRscript(optHash["--dirbase"], optHash["--pheno"], optHash["--snp"], optHash["--geno"], perm1, optHash["--out"], optHash["--gene"])

#   parray=collectResult(optHash["--out"])
  flag = 0
  parray.each do |p|
    if p.to_f < pcut
      flag = 1
      break
    end
  end

  if flag == 1 # do further permuation
    runRscript(optHash["--dirbase"], optHash["--pheno"], optHash["--snp"], optHash["--geno"], perm2, optHash["--out"], optHash["--gene"])
  end
end

def runRscript(base, pheno, snp, geno, perm, out, gene)
  
  cmd="Rscript --slave #{base}/rareVariantTests.R -p #{perm}  -a #{pheno} -b #{snp} -c #{geno} --multicore  > #{out}"
  $stderr.puts cmd
  system(cmd)
  return collectResult(out, gene)
end

def collectResult(rfile,gene)
  flag=0
  parray = []
  File.new(rfile,'r').each do |line|
    if line=~ /^p-values/
      flag = 1
      
    elsif flag == 1 # start
      cols = line.chomp.split(/\s+/)
      if cols[0] == "1"  # p-value line
        cols[1..-1].each do |p|
          parray << p
        end
      end
    end
    
  end
  
  
  if parray.size == 8
    info=File.new(rfile + ".stats", 'w')
    info.puts "#{gene}\t#{parray.join("\t")}"
    info.close
  end
  return parray
end


def getopt
  
  opts = GetoptLong.new(
                        ["--pheno", "-a", GetoptLong::REQUIRED_ARGUMENT],
                        ["--snp", "-b", GetoptLong::REQUIRED_ARGUMENT],
                        ["--geno", "-c", GetoptLong::REQUIRED_ARGUMENT],
                        ["--dirbase", "-d", GetoptLong::REQUIRED_ARGUMENT],
                        ["--gene", "-g", GetoptLong::REQUIRED_ARGUMENT],
                        ["--out", "-o", GetoptLong::REQUIRED_ARGUMENT],
                        ["--perm", "-p", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") 
    $stderr.puts "Usage: ruby __.rb -a pheno -b snp -c geno -o ouput -g geneName -d basedir"
    exit
  end
  return optHash
  
end

main()
    
## header: 
#geneName  CMC_1%  CMC_1%+polyphen2 CMC_5%  CMC_5%+polyphen2  Madson-Browning  Madson-Browning+polyphen2  VT-test VT-test+polyphen2
