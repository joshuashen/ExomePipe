bamlist = ARGV[0]

bamfiles = []

File.new(bamlist, 'r').each do |line|
  bamfiles << line.chomp
end


1.upto(24) do |i|
  if i == 23 
    i = "X"
  elsif i == 24
    i = 'Y'
  end

  system("mkdir -p temp/chr#{i}")
  out=File.new("joint_call.chr#{i}.sh", 'w')
  out.puts "\#!/bin/bash \n\#\$ -cwd"
  
  cmd = "java -Xmx8000m -Djava.io.tmpdir=temp/chr#{i}/  -jar /ifs/home/c2b2/ip_lab/yshen/data/shares/SOFTWARE/Sting/dist/GenomeAnalysisTK.jar  -T UnifiedGenotyper  -R /ifs/home/c2b2/ip_lab/yshen/data/shares/DATA/Sequencing/resources/bcm_hg18.fasta -D /ifs/home/c2b2/ip_lab/yshen/data/shares/DATA/Sequencing/resources/dbsnp_130_hg18_bcm.rod -nt 2 -o snps.joint.chr#{i}.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 50 -L chr#{i}"

  bamfiles.each do |bam|
    cmd = cmd + " -I #{bam}"
  end
  out.puts cmd
  out.close
end


