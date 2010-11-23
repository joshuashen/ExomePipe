# jobs on the cluster may be terminated due to various reasons
# including longer time than allocated, java runtime errors etc

# this script looks up the submitted jobs and jobs that are running,
# determine if completed jobs really did the work, and resubmit otherwise

## check jobs status : qstat -r 

## For now, only check truncated jobs (due to time overrun)




def main
  checklist = ARGV[0]
  global = ARGV[1]
  refix = ARGV[2]

  File.new(checklist, 'r').each do |line|
    line.chomp
    if line=~ /^(\S+)\.(\w*\d+)$/  # match pattern: "TCGA-25-1323-01A-01W-0494-09_IlluminaGA-DNASeq_exome.bam.19"
      bam, chr = $1, $2
      
      dir= bam+"_pipe"
      
#       Dir.chdir(dir); Dir.pwd

      fixed = dir + "/" + chr + ".fixed.bam"
      cleaned = dir + "/" + chr + ".cleaned.bam"
      
      jobstdout = dir + "/" + chr + ".realign.o"
      jobstderr = dir + "/" + chr + ".realign.e"
      
#      puts "#{jobstdout}\t#{jobstderr}"
      cmd="grep \'net.sf.picard.sam.FixMateInformation done\' #{jobstderr}"
      complete = `#{cmd}`
#      puts cmd
      if complete == ""  # not complete in fixmate
        # then check if realign is complete
        runtimec  = `grep \'Total runtime\' #{jobstdout}`  # 
#        puts "#{bam}\t#{runtimec}"
        
        if runtimec.split(/\n/).size == 2  # realign is finished
          $stderr.puts "#{bam}\t#{runtimec.split(/\n/).size}"
          # therefore we only need to re-do fixmate
          cmd="qsub -N #{bam}.#{chr}.fixmate -l mem=3.5G,time=28:: -o #{jobstdout}.new -e #{jobstderr}.new #{refix} -I #{bam} -L #{chr} -g #{global}"
          puts cmd
          system(cmd)
        else
          $stderr.puts "need to start from scratch for #{bam}"
        end
      end
    end
  end
end



main()
