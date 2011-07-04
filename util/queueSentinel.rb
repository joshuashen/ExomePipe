# get the job ID or running qjobs

$VERBOSE=nil

ENV['SGE_ROOT'] = '/opt/gridengine/titan'

def main 
  check = ARGV[0]
  
  if check == nil ## do resub
    doit = 1
  else
    doit = 0  # just check, don't resub
  end

  run = {}   ## { jobID => nodeInfo }
  list = `/opt/gridengine/titan/bin/lx24-amd64/qstat`
  list.each_line do |line|
#    $stderr.puts line
    cols = line.chomp.split(/\s+/)
    next if cols[2] == "QRLOGIN"  # 
    if cols[4] == 'r'  ## running
      run[cols[0]] = cols[7] # job ID
    end
  end
  
  pwd = Dir.pwd()
  
  run.each do |qid, node|
    bdir, o, e = nil, nil, nil
    cwd = nil
    hour,min,sec = 0, 0, 0
    
    f = `/opt/gridengine/titan/bin/lx24-amd64/qstat -j #{qid}`
    f.each_line do |line|
      cols = line.chomp.split(/\s+/)
      if line.match("sge_o_workdir") 
        bdir = cols[1]
      elsif line.match("cwd:")
        cwd = cols[1]
      elsif line.match("stderr_path_list")
        e = cols[1].gsub("NONE:", "")
      elsif line.match("stdout_path_list")
        o = cols[1].gsub("NONE:", "")
      elsif line =~  /cpu\=(\d+)\:(\d+)\:(\d+)\,/
        hour,min,sec = $1.to_i, $2.to_i, $3.to_i
        
      end
    end

    if o.match(/^\//) #  abs path
      jobout = o
    else
      jobout = "#{bdir}/#{o}"
    end

    if e.match(/^\//) 
      joberr = e
    else
      joberr = "#{bdir}/#{e}"
    end

    if !File.exist?(jobout) and !File.exist?(joberr) and  hour == 0 and min == 0 and sec == 0
      $stderr.puts "#{qid} is stuck on #{node}"
 
      if doit == 1
        ## need to goto the working dir
        Dir.chdir(cwd)
        system("/opt/gridengine/titan/bin/lx24-amd64/qresub #{qid}")
        system("/opt/gridengine/titan/bin/lx24-amd64/qdel #{qid}")
      end
    end
    
  end
  
  Dir.chdir(pwd)
  
# 3102451 0.01530 realign.1. saec         r     07/02/2011 21:12:07 pow18.q@b14b.c8.titan              1        
end

main()
