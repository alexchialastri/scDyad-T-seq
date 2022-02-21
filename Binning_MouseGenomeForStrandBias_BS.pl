# Binning gDNA reads into bins

use warnings;


#@chrlen = (248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415);      #Human, chr23 is chrX, chr24 is Y
@chrlen = (195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698,16299);      #chr22 is chrM

$l = scalar(@chrlen); #length of vector

#opendir("C:/Users/alexc/Documents/Research/SisterChromatidExchange/scAba-Seq-OldDatabySid");
#open($input, "C:/Users/alexc/Documents/Research/SisterChromatidExchange/scAba-Seq-OldDatabySid/E14HYDRNA_AbaReads_Full5hmC_Final_Rmdup_sortFirstTwoCells.txt");
open($input, $ARGV[0]);
$inName = substr($ARGV[0],0,-4);
#open($input, "E14HYDRNA_AbaReads_Full5hmC_Final_Rmdup_sort.txt");


#Creates bins and then numbers them in order starting with chromosome 1.
$chrstart[0] = 0;
$numbCells=96;
$binsize = $ARGV[1];
$ReadsNeed=50;
#$cutoff=-1;
#$numbChrom=21;

if ($ARGV[1] > 248956422){
$ARGV[1]="Full_Chr";
$out1 = ">".$inName."_MethPercent_BS_".$ARGV[1].".txt";
$out2 = ">".$inName."_Meth_".$ARGV[1].".txt";
$out3 = ">".$inName."_Unmeth_".$ARGV[1].".txt";
}
else{
$out1 = ">".$inName."_MethPercent_BS_".$ARGV[1]."bp.txt";
$out2 = ">".$inName."_Meth_".$ARGV[1]."bp.txt";
$out3 = ">".$inName."_Unmeth_".$ARGV[1]."bp.txt";
}

open($output, $out1);
open($output2, $out2);
open($output3, $out3);

for ($i=0; $i<$l; $i++)
{
    $chrbin[$i] = int($chrlen[$i]/$binsize); #number of bins in each chromosome
    if ($i > 0)
    {
	$chrstart[$i] = $chrstart[$i-1] + $chrbin[$i-1] + 1; #running list of all bins that tells if you are at a certain bin, what chromosome are you in.
    }
    for ($j=$chrstart[$i];$j<=($chrstart[$i]+$chrbin[$i]);$j++) #for inside the chromosome
    {
	$bincoord[$j]  = ($j - $chrstart[$i])*$binsize + 1;  #gives all the bins in that chromosome
	$chrcolumn[$j] = $i + 1; #cooresponding chromosome
	
    }
    # print " $chrlen[$i] $chrbin[$i] $chrstart[$i]/n";
}


#initializes empty matrix for each bin
for ($i=0;$i<=($chrstart[$l-1]+$chrbin[$l-1]);$i++) #Gives a row for every bin
{
    for ($j=0;$j<$numbCells;$j++)
    {
	$ForwardBinData[$i][$j] = 0; #Reused this variable from old script and it will be our Methylated Counter per cell per bin
	$ReverseBinData[$i][$j] = 0; #Same but for unmethylated
	$BinData[$i][$j] = 0;
    }
}

while ($in = <$input>)
{
    @R= split(" ",$in);
	

    if ($R[1] == 23) {$R[1] = 20;}
    if ($R[1] == 24) {$R[1] = 21;}
    if ($R[1] == 25) {$R[1] = 22;}
	

	
	
	$cell= $R[0] - 1; #because perl starts counting at zero
	$chr = $R[1]-1;
	
	$bin  = $chrstart[$chr] + int($R[2]/$binsize);
	if ($R[5] eq "Z")
	{
	$ForwardBinData[$bin][$cell] = $ForwardBinData[$bin][$cell] + 1;
	}
	
	if ($R[5] eq "z")
	{
	$ReverseBinData[$bin][$cell] = $ReverseBinData[$bin][$cell] + 1;
	
	}

}


for ($i=0;$i<=($chrstart[$l-1]+$chrbin[$l-1]);$i++)
{
    print $output ("$chrcolumn[$i] $bincoord[$i] ");
    print $output2 ("$chrcolumn[$i] $bincoord[$i] ");
    print $output3 ("$chrcolumn[$i] $bincoord[$i] ");

 for ($j=0;$j<$numbCells;$j++)
    {
	print $output2 ("$ForwardBinData[$i][$j] ");
	print $output3 ("$ReverseBinData[$i][$j] ");
		
	if ($ForwardBinData[$i][$j]+$ReverseBinData[$i][$j]>$ReadsNeed) #set some limits on how many reads are needed in a bin to count
	{
	
	$BinData[$i][$j]=$ForwardBinData[$i][$j]/($ForwardBinData[$i][$j]+$ReverseBinData[$i][$j]);
	  
        print $output ("$BinData[$i][$j] ");
	    
		
		}
	else {
		print $output ("NAN "); #These are poorly Read areas
		
	}
    }   
    print $output ("\n");
	print $output2 ("\n");
	print $output3 ("\n");
}


 close($input);

 
 
 
 
 
