
if (ARGV.size() == 0)  then 
  print "Error: no intensity value\n"; 
  exit
end

$colorTable = Array.new(6,0);

$colorTable[0] = [0,0,255];
$colorTable[1] = [255,0,255];
$colorTable[2] = [255,127,0];
$colorTable[3] = [255,255,0];
$colorTable[4] = [127,255,0];
$colorTable[5] = [0,255,0];

def findColor(normalizedDC)
  intensity = 1.0;
  tableindex = normalizedDC*($colorTable.size()-1)  
  firstcolor = tableindex.floor();

  if (firstcolor==$colorTable.size()-1)
    intensity = 1.0;
    firstcolor-=1;
  else
    intensity = tableindex - firstcolor;
  end
  
  secondcolor = firstcolor+1;
  
  color = [(1.0-intensity)*$colorTable[firstcolor][0] + intensity*$colorTable[secondcolor][0],
           (1.0-intensity)*$colorTable[firstcolor][1] + intensity*$colorTable[secondcolor][1],
           (1.0-intensity)*$colorTable[firstcolor][2] + intensity*$colorTable[secondcolor][2]];
           
  
  print "#{color[0].floor()} #{color[1].floor()} #{color[2].floor()}\n"
end


findColor(ARGV[0].to_f)

