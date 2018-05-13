use strict;


my $assembledGenomeFile = $ARGV[0];
my $dwgsimPath = $ARGV[1];
my $sgaPath = $ARGV[2;
my $outPrefix = $ARGV[3];



system("$dwgsimPath -C 20 -e 0 -E 0 -d 400 -1 127 -2 127 -r 0 -y 0 $assembledGenomeFile $outPrefix");
system("gunzip $outPrefix.bwa.read1.fastq.gz");
system("gunzip $outPrefix.bwa.read2.fastq.gz");

system("$sgaPath preprocess --pe-mode=1 --out=$outPrefix .$outPrefix.bwa.read1.fastq $outPrefix.bwa.read2.fastq");
system("$sgaPath index -a ropebwt $outPrefix");
system("$sgaPath correct $outPrefix");
system("$sgaPath filter $outPrefix");
system("$sgaPath overlap -t 20 -m 50 $outPrefix");
system("$sgaPath assemble $outPrefix.asqg.gz");

