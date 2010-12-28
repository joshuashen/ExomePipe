bamlist = ARGV[0]
settingf = ARGV[1]
vm = ARGV[2]

bamfiles = []
setting = {}
outprefix = File.basename(bamlist) + "_SNV_joint_"
reftype = "hg"  # default hg18 or hg19
chrprefix = 'chr'
heap = 8000 # default heap size in M

if vm != nil
  heap = vm.to_i
end

File.new(bamlist, 'r').each do |line|
  bamfiles << line.chomp
end

File.new(settingf, 'r').each do |line|
  if line =~ /REF\=(\S+)/
    setting["ref"] = $1.gsub('"', '')
  elsif line =~ /DBSNP\=(\S+)/
    setting["dbsnp"] = $1.gsub('"', '')
  elsif line =~ /GATKJAR\=(\S+)/
    setting["gatk"] = $1.gsub('"', '')
  elsif line =~ /REFTYPE\=(\S+)/
    reftype = $1.gsub('"', '')    
  end
end

if reftype != 'hg'
  chrprefix = ''
end

1.upto(24) do |i|
  if i == 23 
    i = "X"
  elsif i == 24
    i = 'Y'
  end

  system("mkdir -p temp/chr#{i}")
  out=File.new("#{outprefix}chr#{i}.sh", 'w')
  out.puts "\#!/bin/bash \n\#\$ -cwd"
  
  cmd = "java -Xmx#{heap}m -Djava.io.tmpdir=temp/chr#{i}/  -jar #{setting["gatk"]} -T UnifiedGenotyper  -R #{setting["ref"]}  -D #{setting["dbsnp"]}  -nt 2 -o #{outprefix}chr#{i}.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 50 -L #{chrprefix}#{i}"

  bamfiles.each do |bam|
    cmd = cmd + " -I #{bam}"
  end
  out.puts cmd
  out.close
end


