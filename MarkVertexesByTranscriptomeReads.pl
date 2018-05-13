use strict;
use Bio::SeqIO;


my $graphBCalmFile = $ARGV[0];
my $mappedBGREATFile = $ARGV[1];
my $maxVertexSeqLen = $ARGV[2];
my $minCovNumber = $ARGV[3];
my $outFile = $ARGV[4];
my $outCoveredVertexies = $ARGV[5];


my %graphEdges = ();
my %verSeqLen = ();
my %verSeq = ();
my %coveredVertexies = ();



open(FTR,"<$graphBCalmFile");


open(FTW,">$outFile") or die;
my $outFasta = Bio::SeqIO->new(-file=>">$outCoveredVertexies",-format=>'fasta');


while (my $input = <FTR>) {
    #print "$input";
    chomp($input);
    my $val = substr($input,0,1);
    
    if (substr($input,0,1) eq ">") {
        $input = substr($input,1);
        
        my @arrInp = split(/\s+/,$input);
        my $verNam = $arrInp[0];
        die if(exists $graphEdges{"$verNam"});
        my @arrNew = ();
        $graphEdges{"$verNam"} = \@arrNew;
        for(my $i = 0;$i<=$#arrInp;$i++)
        {
            my $currVal = $arrInp[$i];
            
            if(substr($currVal,0,2) ne "LN")
            {

                if(substr($currVal,0,1) eq "L")
                {
                    my @arrTmp = split(/\:/,$currVal);
                    my $conVertex = $arrTmp[2];
                    my $ptrArrNew = $graphEdges{"$verNam"};
                    push @$ptrArrNew,$conVertex;
                    
                }
            }
            
        }
        
        
    }
    

}

close(FTR);


my $fasta=Bio::SeqIO->new(-file=>"$graphBCalmFile",-format=>'fasta');
while (my $seq=$fasta->next_seq()) {
	my $seqId = $seq->id;
	my $seqTxt = $seq->seq;
    
    $verSeqLen{"$seqId"} = length($seqTxt);
    $verSeq{"$seqId"} = $seq;
    
}
$fasta->close();

open(FTR,"<$mappedBGREATFile") or die;

while(my $input = <FTR>)
{
    chomp($input);
    $input = substr($input,0,length($input)-1);
    if(substr($input,0,1) eq "-")
    {
        $input = substr($input,1);
    }
    
    die if(not exists $verSeqLen{"$input"});
    
    next if($verSeqLen{"$input"} > $maxVertexSeqLen);
    
    if(exists $coveredVertexies{"$input"})
    {
        $coveredVertexies{"$input"} = $coveredVertexies{"$input"} + 1;
    }
    else
    {
        $coveredVertexies{"$input"} = 1;
    }
    
}


close(FTR);

my @arrCovered = keys %coveredVertexies;

for(my $i = 0;$i <= $#arrCovered;$i++)
{
    my $verCov = $arrCovered[$i];
    my $verCovCount = $coveredVertexies{"$verCov"};
    
    if($verCovCount > $minCovNumber)
    {
        print FTW "$verCov\n";
        die if(not exists $verSeq{"$verCov"});
        my $currVerSeq = $verSeq{"$verCov"};
        $outFasta->write_seq($currVerSeq);
        next;
    }
    
    die if(not exists $graphEdges{"$verCov"});
    
    my $ptrArr = $graphEdges{"$verCov"};
    my $isFoundCovered = 0;
    for(my $i = 0;$i<=$#$ptrArr;$i++)
    {
        my $currVertex = $ptrArr->[$i];
        if(exists $coveredVertexies{"$currVertex"})
        {
            $isFoundCovered = 1;
        }
    }
    
    if($isFoundCovered == 1)
    {
        print FTW "$verCov\n";
        die if(not exists $verSeq{"$verCov"});
        my $currVerSeq = $verSeq{"$verCov"};
        $outFasta->write_seq($currVerSeq);
    }
    
}

close(FTW);
$outFasta->close();