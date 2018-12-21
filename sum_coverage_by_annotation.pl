#!/usr/bin/perl
#Eric Morrison
#121918

use strict;
use warnings;

sub process_coverage{
	my $coverage = $_[0];
	open(COV, "$coverage") || die "can't open coverage file\n";
	my %cov;
	#hash of arrays with seq ID index and [total bases mapped, tot. mapped/seq. length]
	while(my $cov = <COV>){
		chomp($cov);
		my @cov = split("\t", $cov);
		$cov{$cov[0]} = [$cov[2], $cov[3]];
	}
	return(\%cov);
}

sub add_cov_to_anns{
	my($ann, $covRef, $categoryInd) = @_;
	open(ANN, "$ann") || die "Can't open annotation file\n";
	my %cov = %$covRef;
	#hash of annotation categories (unique names such as KO number or taxonomy
	#add cov to hash val by calling coverage hash as reading annotation file
	my %annCov;
	while(my $ann = <ANN>){
		chomp($ann);
		my @ann = split("\t", $ann);
		
		if(defined($annCov{$ann[$categoryInd]}) == 0 && defined($cov{$ann[0]}) == 1){
			$annCov{$ann[$categoryInd]}{"total"} = ${ $cov{$ann[0]} }[0];
			$annCov{$ann[$categoryInd]}{"avg"} = ${ $cov{$ann[0]} }[1]
		}elsif(defined($cov{$ann[0]}) == 1){
			$annCov{$ann[$categoryInd]}{"total"} += ${ $cov{$ann[0]} }[0];
			$annCov{$ann[$categoryInd]}{"avg"} += ${ $cov{$ann[0]} }[1];
		}
	}
	return(\%annCov);
}

sub print_anns{
	my $annRef = $_[0];
	my %anns = %$annRef;
	
	foreach my $function (sort{$a cmp $b} keys %anns){
		print OUT $function, "\t", ${ $anns{$function} }[0], "\t", ${ $anns{$function} }[1], "\n";
	}
}
#MAIN
{
	my $coverage = $ARGV[0];#Coverage by gene sequence
	my $annotations = $ARGV[1];#gene annotations indexed by sequence ID
	my $categoryInd = $ARGV[2]; #row number of desired category to summarise, for example KO number or taxonomy, (starting from row 0)
	
	my $covHashRef = process_coverage($coverage);
	my $annHashRef = add_cov_to_anns($annotations, $covHashRef, $categoryInd);
	print_anns($annHashRef);
}
