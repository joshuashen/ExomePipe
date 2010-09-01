# step 0: bwa index of reference genome

# command: bwa index -a bwtsw human_g1k_v37.fasta 

ref = ARGV[0]

cmd = "bwa index -a bwtsw " + ref
system(cmd)
 
