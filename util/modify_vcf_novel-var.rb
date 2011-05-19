## add names to novel variants -- output from older version of GATK

# example: 
#1       865665  .       G       A       104.85  PASS    AC=2;AF=0.0029;AN=680;DP=2362;Dels=0.00;HRun=0;HaplotypeScore=0.5121;MQ=59.32;MQ0=0;QD=17.41;SB=-29.01;refseq.changesA...
# change to:
#1       865665  var_1_865665       G       A       104.85  PASS    AC=2;AF=0.0029;AN=680;DP=2362;Dels=0.00;HRun=0;HaplotypeScore=0.5121;MQ=59.32;MQ0=0;QD=17.41;SB=-29.01;refseq.changesA...


while line=ARGF.gets do 
  cols=line.split(/\t+/)
  if cols[2] == "."
    cols[2] = "var_" + cols[0]+"_" + cols[1]
  end
  puts cols.join("\t")
end


