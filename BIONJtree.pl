use strict;

for(my $i=0;$i<1000;$i++)
{
	print "$i\n";
	`fastme -i /home/karine_user/Documents/Spodo/VCF_11052020_phylo/phylo_092020/boots/$i.dist -o /home/karine_user/Documents/Spodo/VCF_11052020_phylo/phylo_092020/NJtree/boottree/$i`;
}
 



