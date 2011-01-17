## calcualte transition/transversion ratio from VCF file


def main
  novelti = 0
  noveltv = 0
  knownti = 0
  knowntv = 0
  
  while line=ARGF.gets do 
    cols=line.chomp.split(/\s+/)
    name,ref,alt = cols[2],cols[3],cols[4]
    ti = judge(ref,alt)
    if name =~ /^rs/  # known
      if ti == 1
        knownti += 1
      else
        knowntv += 1
      end
    else  # novel
      if ti == 1
        novelti += 1
      else
        noveltv += 1
      end
    end
  end
  puts "#knownti/tv: #{knownti}/#{knowntv} = #{knownti/knowntv.to_f}"
  puts "#novelti/tv: #{novelti}/#{noveltv} = #{novelti/noveltv.to_f}"
end

def judge(a1, a2)
  ti = 0
  if a1 == 'A'
    if a2 == 'G'  # ti
      ti = 1
    end
  elsif a1 == 'C'
    if a2 == 'T'
      ti = 1
    end
  elsif a1 == 'T'
    if a2 == 'C'
      ti = 1
    end
  elsif a1 == 'G'
    if a2 == 'A'
      ti = 1
    end
  end
  return ti
  
end

main()

