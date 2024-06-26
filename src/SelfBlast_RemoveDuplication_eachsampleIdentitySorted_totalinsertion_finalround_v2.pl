#!/usr/bin/perl
use strict;
use warnings;
#use Data::Dump qw(dump);
#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### function: This script is to extract the potential eliminated insertions and added the removed insertion coverage into previous insertions
### The criterion of the script is as follows,
###   1. we considered the short insertions (short 50bp, 150 bp in total) have a high strandard for final deduplications as they have less sequence errors, 95% identity, 95% Coverage, less than 6 gapsize and less than 6 mismatches as previously showed (no need to set up as before).

###  2 for the longer insertions that happened very few, we set up less stringent due to high sequence errors. The main setup is the matches/ read length >0.95, identity >80%, wherease we ignore the gapsize, mismatches. These parameter is adjustable if we try to have a high standard.




my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","o:s","i:s","b:s","h:s","t:s","g:s","c:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} ||!defined $opts{o} ||!defined $opts{i}||!defined $opts{b}|| defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -f Insfasta -i Read with evaluation -b SelfBlast -o Output of Quality Control Reads
	
	Request Parameters:
	-b Blast Results (Self Blast results))
	-i The evaluation of each inserted reads(reads count, quality, identity)
	-f Inserted fasta (The final fasta with insertion events)
	-o The final results strings of files, including the final forward/reverse deduplicated reads, final clustering file with statistical resutls and their representive read
	
	Optional Parameters for long insertion because of the low quality within sequence length
	-r Cut off whole long insertions because of the read length(default 300)
	-t Identity of two reads (default 80)
	-g Shift size (default 6)
	-c Coverage (Matched size/Full length, default 95%)
	-h Help
************************************************************************\n";
}



########################################################################################
#### The first step to generate a unique cluster results from the self blast results

#### Criterion : Identity more than 90%, less than , two reads coverage more than 95% for short insertions as previous 

########################################################################################

my $blast=$opts{b};

my $ident=(defined $opts{t})?$opts{t}:80;
my $cover=(defined $opts{c})?$opts{c}:0.95;
#my $mismatches = (defined $opts{m})?$opts{m}:10;
my $shiftsize=(defined $opts{g})?$opts{g}:6;
my $output=$opts{o};
my $Lengthcutoff=(defined $opts{r})?$opts{r}:300;


my %remove; my %con; my %contain; my %hash; my %name; my $n=0;

open BLAST,"$blast" or die $!;
while (<BLAST>){
	chomp;
	my ($id1,$id2,$qidentity,$qmatches,$qmismatches,$qgapsize,$qlenth,$rlength)=(split/\t/,$_)[0,1,2,3,4,5,8,11];
	next if ($id1 eq $id2);
	my $av=$qmatches/$qlenth;
	my $av2=$qmatches/$rlength;
	my $maxindel=0.15*$qlenth;
	
	### considering the begining of high quality at the upstream, here we restricted the shorter insertion events with high standard:
 
	### we ignore the dramatic sequence error but require the difference of overall alignments 
	if ($qlenth<$Lengthcutoff && $rlength<$Lengthcutoff ){
		
		next unless ($av>=0.95 && $av2>=0.95 && $qidentity >=90 && $qmismatches<=$maxindel && $qgapsize<=$maxindel);

	}else{
		next unless ($av>=$cover && $av2>=$cover && $qidentity >=$ident && abs($qlenth-$rlength) <=$shiftsize && $qmismatches<=$maxindel && $qgapsize<=$maxindel );
	}
	
	# next unless ($av>=$cover && $qidentity >=$ident && abs($qlenth-$rlength) <=$shiftsize);
	#
	# ### here we restricted the shorter insertion events with high standard:
	# next unless ($qlenth<140 && $av>=$cover && $qidentity >=$ident && abs($qlenth-$rlength) <=$shiftsize)
	
	if (!exists $name{$id1} && !exists $name{$id2}){
		$n++;
		#push @{$name{$id1}},($id1,$id2);
		my $string=join "\t",($id1,$id2);
		$name{$id1}=$n;
		$name{$id2}=$n;
		$hash{$n}=$string;
	}elsif (exists $name{$id1} && !exists $name{$id2} ){
		my $num=$name{$id1};
		#push @{$name{$id1}},$id2;
		$name{$id2}=$num;
		$hash{$num}.="\t$id2";
	}elsif (exists $name{$id2} && !exists $name{$id1} ){
		my $num=$name{$id2};
		#push @{$name{$id2}},$id1;
		$hash{$num}.="\t$id1";
		$name{$id1}=$num;
	}elsif (exists $name{$id2} && exists $name{$id1} && $name{$id2} != $name{$id1}){
		#print "$name{$id2}\t$name{$id1}\n";
		#identify the min value and, add the keys from max value to min value and minumium keys,  delete the max hash and value
		my ($max,$min)=($name{$id2} > $name{$id1})?($name{$id2},$name{$id1}):($name{$id1},$name{$id2});
		my @arrayR=split/\t/,$hash{$max};
		foreach my $q (@arrayR){
			$name{$q}=$min;
			$hash{$min}.="\t$q";
		}
		delete $hash{$max};
	}	
}


######################################################################################################################################
### make the index files: 
### 		asssembled fasta with the large insertion events
######################################################################################################################################

### read the assembled fasta
my $Assembled=$opts{f};
open FASTA,"$Assembled" or die $!;
#
my $Fid;  my %sequence; my %hashF; my %inf;
while (<FASTA>) {
	chomp;
    # print "$_" ;
	if (/^>(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
		$Fid=$1;
		$hashF{$Fid}++;
		$inf{$Fid}->{cov}=$2;
		$inf{$Fid}->{quality}=$3;

		### we define the identiy as 0 if there is no alignments
		$inf{$Fid}->{identity}=$4;		
	}else{
		$sequence{$Fid}.=$_;
		#$len{$Fid}=length $_;
	}

}
close FASTA;


#### Read the statistical file with the identify and quality 
my $rqual=$opts{i};  my $type;

### this file generated by the (yYY423-2_S13_mappable.txt)
open RQU,"$rqual" or die $!;
while (<RQU>){
	chomp;
	my ($id,$iden,$Rcount,$qual)=(split/\t/,$_)[1,7,8,9];
	next if ($id eq "NID");
	$type=(split/\t/,$_)[2];
	$inf{$id}->{inf}.=$_."\n";
}

close RQU;


######################################################################################################################################
## Here we generated the clustering results for further eliminating the duplicates and determine the represent read for each cluster
## Criterion: 1. Ranking the clustering reads with the identiy, quality and Read Count. We used the rank first read to represent each cluster
###
######################################################################################################################################


open CLS,">$output.cls" or die $!;

print CLS "AsemStatus\tInsertionID\tReadNumber\tRepRead\tRepQual\tRepIden\tClsReads\tClsIden\tClsQual\tRankNum\n";

my %uniq; my %repeat; my %frcounts;
foreach my $i (sort keys %hash){
	
	my @array=split /\t/,$hash{$i};
	my $clsN=@array;
	print CLS "Assembly\t$i\t";
	my @infQ; my @infI; my %qual=(); my %iden=();my %scoreI=();my %scoreQ=(); my %cov=();
	
	
	### here is to define the quality, identity of the clustering reads	
	foreach my $p (@array){
		
		$qual{$p}=$inf{$p}->{quality};
		$cov{$p}=$inf{$p}->{cov};
		### We set up the no-aligned identity as 0, 
		### which would be better for the identity comparison that we can put them into a pretty low priority.
		### If all the groups cannot have the identity, we then only consider the inserted quality.
		$iden{$p}=($inf{$p}->{identity} eq "NA")?0:$inf{$p}->{identity};
		
	}

	
	
#### Here we rank the reads first with quality , and then read count support. In this case, we obtained the representive result with the best alignment and high quality	

	my @keyF = sort {$cov{$b} <=> $cov{$a} || $qual{$b} <=> $qual{$a} } keys %cov;

	#### push the representive read to hash, and put the number of reads from each cluster into hash value
	$uniq{$keyF[0]}=$clsN; 
	
	# # we sorted the read quality first
# 	my @keyQ = sort { $qual{$b} <=> $qual{$a} or $a cmp $b } keys %qual;

	my ($prev,$prevC, $rankQ);
	for my $k (@keyF) {
	    $rankQ++ unless defined($prev) && $prev==$qual{$k} && $prevC == $cov{$k};
		$scoreQ{$k}=$rankQ;
	    $prev = $qual{$k};
		$prevC =$cov{$k};
	}
	
	
	my $num=0; my $Fstring; my $Fidentity; my $FQual; my $SortID;my $FCounts=0; my $FCov;
	
	#### 
	
	foreach my $m(@keyF){
		$num++;
		$Fstring.="$num;";
		$Fidentity.="$inf{$m}->{identity};";
		$FQual.="$inf{$m}->{quality};";
		$FCov.="$inf{$m}->{cov};";
		$repeat{$m}++;
		$SortID.="$m;";
		$FCounts +=$inf{$m}->{cov};	
	}
	
	$frcounts{$keyF[0]}=$FCounts;
	print CLS "$FCounts\t$keyF[0]\t$inf{$keyF[0]}->{quality}\t$inf{$keyF[0]}->{identity}\t$SortID\t$Fidentity\t$FQual\t$FCov\t$Fstring\n";
}
	

######################################################################################################################################
## Here we generated the unique fastq files for large insertion
######################################################################################################################################

open FOUT,">$output.FAssuniq.fasta" or die $!;
open OUT,">$output.finalinsertion.txt" or die $!;
#open ROUT,">$output.Funiq.reverse.fastq" or die $!;
my $uncls=0; my $un=0;
foreach my $f (keys %hashF){
	
	if (exists $uniq{$f}){
		print FOUT ">$f\n$sequence{$f}\n";
		
		if (exists $inf{$f}->{inf}){		
			$inf{$f}->{inf}=~s/\n$//;
			print OUT "$inf{$f}->{inf}\n";
			
		}else{
			$un++;
			my $seq=$sequence{$f};
			my $lengthun=length($seq) -90;
			my $insertedseq=substr $seq, 44,$lengthun;
			print OUT "Un$un\t$f\t$type\t$seq\t$insertedseq\t$lengthun\t1\tUnknown\t$inf{$f}->{cov}\t$inf{$f}->{quality}\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\n";
		}
	
		
	}elsif(!exists $uniq{$f} && !exists $repeat{$f} ){
		$uncls++;
		my $Fcoverage=$inf{$f}->{cov};
		print FOUT ">$f\n$sequence{$f}\n";
		print CLS "Assembly\tUndefined$uncls\t$Fcoverage\t$f\t$inf{$f}->{quality}\t$inf{$f}->{identity}\tNoClustered\t$inf{$f}->{identity}\t$inf{$f}->{quality}\tNoRank\n";
		
		if (exists $inf{$f}->{inf}){
			
			$inf{$f}->{inf}=~s/\n$//;
			print OUT "$inf{$f}->{inf}\n";
		}else{
			$un++;
			my $seq=$sequence{$f};
			my $lengthun=length($seq) -90;
			my $insertedseq = substr $seq, 44, $lengthun;
			print OUT "Un$un\t$f\t$type\t$seq\t$insertedseq\t$lengthun\t1\tUnknown\t$inf{$f}->{cov}\t$inf{$f}->{quality}\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\n";
		}
		
	}
	
}

close CLS;
close FOUT;
close OUT;



