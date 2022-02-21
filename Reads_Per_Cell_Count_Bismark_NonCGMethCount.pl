# Reads per cell from a FABA file



$filename = $ARGV[0];
open($file, '<', $filename);

$outputfilenameCG      = substr($filename,0,-4)."_MethylatedCountsPerCell.txt";
open($outputCG, '>', $outputfilenameCG);

$outputfilenamezCG      = substr($filename,0,-4)."_UnmethylatedCountsPerCell.txt";
open($outputzCG, '>', $outputfilenamezCG);


for ($j=0;$j<=95;$j++)
    {
	$MethReadsPerCell[$j] = 0;
	$UnMethReadsPerCell[$j] = 0;
    }
	
while ($line = <$file>)
{
    @x = split(" ",$line);
    #$ReadsPerCell[$x[0]-1]++;
    $Cell = $x[0]-1;
	if ($x[1] <= 23){
	if (($x[5] eq H) | ($x[5] eq X)){$MethReadsPerCell[$Cell]++;}
	if (($x[5] eq h) | ($x[5] eq x)){$UnMethReadsPerCell[$Cell]++;}
	}
}

for ($j=0;$j<=95;$j++)
    {
print $outputCG ("$MethReadsPerCell[$j]\n");
print $outputzCG ("$UnMethReadsPerCell[$j]\n");

}
