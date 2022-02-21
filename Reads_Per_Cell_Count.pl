# Reads per cell from a FABA file



$filename = $ARGV[0];
open($file, '<', $filename);

$outputfilenameCG      = substr($filename,0,-5)."_CountsPerCellFaba.txt";
open($outputCG, '>', $outputfilenameCG);


for ($j=0;$j<=95;$j++)
    {
	$ReadsPerCell[$j] = 0;

    }
	
while ($line = <$file>)
{
    @x = split("\t",$line);
    #$ReadsPerCell[$x[0]-1]++;
	$ReadsPerCell[$x[0]]++;
	
}

for ($j=0;$j<=95;$j++)
    {
print $outputCG ("$ReadsPerCell[$j]\n");
}
