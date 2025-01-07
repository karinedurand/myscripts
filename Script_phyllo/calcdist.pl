use strict;

my $inputvcf='/home/karine_user/Documents/Spodo/geneflow_divmigrate/outgroup.vcf.gz';
my $inputgenotype='/home/karine_user/Documents/Spodo/VCF_11052020_phylo/phylo_092020/reducedvcf.txt.gz';
my $resultF="/home/karine_user/Documents/Spodo/VCF_11052020_phylo/phylo_092020/VCF_11052020_phylo/phylo_092020/dist.tbl";

#my @names;
my %dist;
my $n=0;
my $nsamples;
my %results;

open(IN, "zcat $inputvcf |");
L:while (<IN>)
{
	$_=~s/\n//g;
	my @single=split("\t",$_);

	if($_=~/CHROM/)
	{
		for(my $i=0;$i<9;$i++) {shift @single}
		$nsamples=($#single+1);		
	}
	if($_!~/#/) {last L}
}
close (IN);


open(IN, "zcat $inputgenotype | ");
while(<IN>)
{
#		foreach my $s (@single)
#		{
#			my $ns=10-(length $s);
#			for(my $i=0;$i<$ns;$i++) {$s.=" "}
#			push @names,$s;
#		}
#	}
#	else
#	{
#		if($_=~/##/) {next}
#
#		shift @single;
#		shift @single;
#		shift @single;
#		my $ref=shift @single;
#		my $alt=shift @single;
#		if($alt=~/,/) {next}
#		my $tstv=2;
#		if($ref eq 'A' and $alt eq 'G') {$tstv=1}
#		if($ref eq 'G' and $alt eq 'A') {$tstv=1}
#		if($ref eq 'T' and $alt eq 'C') {$tstv=1}
#		if($ref eq 'C' and $alt eq 'T') {$tstv=1}
#		shift @single;
#		shift @single;
#		shift @single;
#		shift @single;
#
#		my @genotypes;
#		for(my $i=0;$i<=$#single;$i++)
#		{
#			$single[$i]=~/^(\d)[\/|\|](\d)/;
#			if($1 eq '') {next L}
#			my $geno=($1+$2)*$tstv;
#			push @genotypes,$geno;
#		}

	$_=~s/\n//;
	my @genotypes=split('',$_);
	for(my $i=0;$i<=$#genotypes;$i++)
	{
		for(my $j=0;$j<=$#genotypes;$j++)
		{
			my $key="$i\_$j";
			my $diff=($genotypes[$i]-$genotypes[$j]);
			$dist{$key}+=($diff*$diff);
		}
	}

	$n++;

	if($n%1000==0)
	{
		savedata();
		print "$n\n";
		undef %dist;
	}
}
close(IN);

my @idx=sort {$a cmp $b} keys %results;
my $res;
foreach my $id (@idx)
{
	$results{$id}=~s/\t$/\n/;
	$res.="$id\t$results{$id}";
}

open my $fd,">$resultF";
print $fd $res;
close $fd;

###########################################

sub savedata
{
	for(my $i=0;$i<$nsamples;$i++)
	{
		for(my $j=0;$j<$nsamples;$j++)
		{
			my $key="$i\_$j";
			$results{$key}.="$dist{$key}\t";
		}
	}

	my @idx=sort {$a cmp $b} keys %results;

	my $res;
	foreach my $id (@idx)
	{
        	$res.="$id\t$results{$id}\n";
	}

	$res=~s/\tt\n/\n/g;

	open my $fd,">$resultF";
	print $fd $res;
	close $fd;
}
		
#	my $(utput="   $nsamples\n";
#	for(my $i=0;$i<$nsamples;$i++)
#{
#	my $res=$names[$i];
#
#	for(my $j=0;$j<$nsamples;$j++)
#	{
#		my $key="$i $j";
#		my $d=int((sqrt $dist{$key})*10000)/10000;
#		$res.="$d ";
#	}	
#	$res=~s/ $/\n/;
#	$output.=$res;
#}

#open my $fd,">$resultF";
#print $fd $output;
#close $fd;

