#!/usr/bin/perl
#Eric Morrison
#062217
#Usage get_maker_mRNA_coverage.pl [gff] [assembly.fasta]
#This script uses a maker gff and a coverage file from 'samtools depth -a' to find the 
#the average coverage of mRNAs. A file of IDs on separate lines can be entered to acquire 
#only those ids or enter 'all' as the last argument to get coverage for all mRNAs. Returns 
#contig ID, mRNA ID, avg. coverage of gene, and gene avg. cov./(mean/median/mode) genome cov

use strict;
use warnings;


sub process_fasta{
	my $fasta = $_[0];
	open(FAS, "$fasta") || die "Can't open fasta file.\n";
	my $id;
	my %fas;
	while(my $fas = <FAS>){
		chomp($fas);
		if($fas =~ /^>/){
			$id = $fas;
			$id =~ s/^>//;
			$fas{$id} = "";
			next;
		}
		$fas{$id} .= $fas;
	}
	return(\%fas);
}

sub process_gff{
	my $gffFile = $_[0];
	open(GFF, "$gffFile") || die "Can't open GFF.\n";
	my %gff;
	while(my $gff = <GFF>){
		chomp($gff);
		#skip header
		if($gff =~ /^#/){
			next;
		}
		my @gff = split("\t", $gff);
		$gff[8] =~ /ID=(.+?);/; #or use non greedy match of characters to ';'
		my $fid = $1; #feature ID
		$fid =~ s/\.//; #feature IDs in gene annotation files from JGI do not use '.' to indicate mutiple features extracted from a contig
		$gff{$fid} = [$gff[0], $gff[3], $gff[4], $gff[6]]; #scaffold/contig ID, start coord., end coord., strand
	}
	return(\%gff);
}

sub extract_genes{
	my($fasRef, $gffRef, $out) = @_;
	my %fas = %$fasRef;
	my %gff = %$gffRef;
	open(OUT, ">$out") || die "Can't open output\n";
	
	foreach my $fid (keys %gff){
		my @fid = @{ $gff{$fid} };
		#print $fid[0], "\n";
		my $seq = substr($fas{$fid[0]}, $fid[1]-1, $fid[2]-$fid[1]+1);
		if($fid[3] == -1){
			$seq = reverse_comp($seq);
		}
	print OUT ">$fid\n$seq\n";
	}
}

sub reverse_comp{
	my $seq = $_[0];
	$seq = reverse($seq);
	$seq =~ tr/AGCTagct/TCGAtcga/;
	return($seq);
}

#MAIN
{
	
	my $gffFile = $ARGV[0];
	my $fasta = $ARGV[1];
	my $out = $ARGV[2];

	my $fasHashRef = process_fasta($fasta);
	my %fas = %$fasHashRef;
	my $gffHashRef = process_gff($gffFile);
	extract_genes($fasHashRef, $gffHashRef, $out);
}
