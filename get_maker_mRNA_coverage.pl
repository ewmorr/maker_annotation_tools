#!/usr/bin/perl
#Eric Morrison
#062217
#Usage get_maker_mRNA_coverage.pl [gff] [coverage] [file of ids or 'all']
#This script uses a maker gff and a coverage file from 'samtools depth -a' to find the 
#the average coverage of mRNAs. A file of IDs on separate lines can be entered to acquire 
#only those ids or enter 'all' as the last argument to get coverage for all mRNAs. Returns 
#contig ID, mRNA ID, avg. coverage of gene, and gene avg. cov./(mean/median/mode) genome cov

use strict;
use warnings;

my $gffFile = $ARGV[0];
my $coverage = $ARGV[1];
my $idSelection = $ARGV[2];

open(GFF, "$gffFile") || die "Can't open GFF.\n";
open(COV, "$coverage") || die "Can't open coverage file.\n";

my %ids;
if($idSelection ne "all")
	{
	open(IDS, "$idSelection") || die "Can't open $idSelection.";
	chomp(my @ids = <IDS>);
	if($ids[0] =~ /\r\n|\r/)
		{
		print "Unknown line feeds. Please correct id file with '\n' as line feed.";
		exit;
		}
	foreach my $id (@ids)
		{
		$ids{$id} = 1;
		}
	}

#hash indexed by contig name with arrays of coverage indexed by position
my %cov;
my%genome;#mean
my@genome;#median
my%valCount;#mode
while(my $cov = <COV>)
	{
	my @cov = split("\t", $cov);
	$cov{$cov[0]}[ $cov[1]-1 ] = $cov[2]; #substract 1 to start from zero index
	$genome{"cov"}+=$cov[2];
	$genome{"count"}++;
	push(@genome, $cov[2]);
	$valCount{$cov[2]}++;
	#print $cov, "\n";
	#push(@{ $cov{$cov[0]}{"position"} }, $cov[1]);
	#push(@{ $cov{$cov[0]}{"cov"} }, $cov[2]);
	}
	
my $medianCov = median(\@genome);
print $medianCov, "\n";

my $modeCov = max(\%valCount);
print $modeCov;

#my %gff;
while(my $gff = <GFF>)
	{
	#skip header
	if($gff =~ /##/)
		{
		next;
		}
	my @gff = split("\t", $gff);
	
	#skip non mRNA
	if(defined($gff[2]) == 0 || $gff[2] ne "mRNA")
		{
		next;
		}

	#get ID
	#$gff[8] =~ /ID=([^;]+)/; #[^character] excludes character from string
	$gff[8] =~ /ID=(.+?);/; #or use non greedy match of characters to ';'
	my $id = $1;
	
	#if id list is provided skip ids that are undefined
	if($idSelection ne "all" && defined($ids{$id}) != 1)
		{
		next;
		}

	#input positions to coverage sub
	my $covSum;
	my $count;
	for(my $i = $gff[3]-1; $i<$gff[4]; $i++)#need to subtract 1 to index properly
		{
		$covSum += $cov{$gff[0]}[$i];
		$count++;
		}
	
	print $gff[0], "\t", $id, "\t", sprintf("%.2f", $covSum/$count), "\t", sprintf("%.2f", (($covSum/$count)/($genome{"cov"}/$genome{"count"})) ), 
		"\t", sprintf("%.2f", (($covSum/$count)/$medianCov) ), "\t", sprintf("%.2f", (($covSum/$count)/$modeCov)), "\n";
	}
		

sub median
{
	my $valsRef = $_[0];
	my @vals = @$valsRef;
    @vals = sort {$a <=> $b} @vals;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub max
{
my $hRef = $_[0];

my $max = 0;
my $ind;
foreach my $key (keys %$hRef)
	{
	if($$hRef{$key} > $max)
		{
		$max = $$hRef{$key};
		$ind = $key;
		}
	}
return $ind;
}



	