use strict;

my $inputvcf='/home/karine_user/Documents/Spodo/geneflow_divmigrate/outgroup.vcf.gz';
my $distF="/home/karine_user/Documents/Spodo/VCF_11052020_phylo/phylo_092020/dist.tbl";
my $resultD="/home/karine_user/Documents/Spodo/VCF_11052020_phylo/phylo_092020/boots";

my $B=1000;
my @names;
my $nsamples;

open(IN, "zcat $inputvcf |");
L:while (<IN>)
{
	$_=~s/\n//g;
	my @single=split("\t",$_);

	if($_=~/CHROM/)
	{
		for(my $i=0;$i<9;$i++) {shift @single}
		foreach my $s (@single)
		{
			my $ns=10-(length $s);
			for(my $i=0;$i<$ns;$i++) {$s.=" "}
			push @names,$s;
		}
		$nsamples=($#single+1);		
	}
	if($_!~/#/) {last L}
}
close (IN);

open my $fd,$distF;
my @ddata=<$fd>;
close $fd;

my @K;
my @dists;
foreach my $line (@ddata)
{
	$line=~s/\n//;
	my @single=split("\t",$line);
	my $key=shift @single;
	push @K,$key;

	for(my $i=0;$i<=$#single;$i++)	{$dists[$i].="$single[$i]\t"}
}

foreach my $line (@dists) {$line=~s/\t$//}

#my $Fres;
for(my $b=0;$b<$B;$b++)
{
	print "$b\n";
	my %DIST;
	for(my $i=0;$i<=$#dists;$i++)
	{
		my @s=split("\t",$dists[int rand $#dists]);
		
		for(my $j=0;$j<=$#s;$j++) {$DIST{$K[$j]}+=$s[$j]}
	}

	my $iutput="   $nsamples\n";
	for(my $i=0;$i<$nsamples;$i++)
	{
		my $res=$names[$i];
		for(my $j=0;$j<$nsamples;$j++)
		{
			my $key="$i\_$j";
			my $d=int((sqrt $DIST{$key})*10000)/10000;
			$res.="$d ";
		}	
		$res=~s/ $/\n/;
		$iutput.=$res;
	}

#	$Fres.="$iutput\n\n";

#open my $fd,">$resultF";
#print $fd $Fres;
#close $fd;

	my $outF="$resultD/$b.dist";

	print "ls $outF\n";
	open my $fd,">$outF";
	print $fd $iutput;
	close $fd;
}

#open my $fd,">$resultF";
#print $fd $output;
#close $fd;

