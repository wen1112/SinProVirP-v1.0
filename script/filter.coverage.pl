#! /usr/bin/perl -w
#use strict;
die "perl $0 [genomecov_file][idxstats_file][idxstats_output_file][coverag_cutoff]\n" unless @ARGV == 4;
my ($in_f,$idx,$out_f,$cut) = @ARGV;
open IN, $in_f or die $!;
open IDX,$idx or die $!;
open OUT, ">$out_f" or die $!;
my %ctgid_coverage;
my $cutoff = $cut;
while(<IN>){
	chomp;
	my $line=$_;
	my @b=split /\t/,$line;
	if($b[1]!=0){
		if(exists $ctgid_coverage{$b[0]}){
			$ctgid_coverage{$b[0]}+=$b[4];
		#	print "$ctgid_coverage{$b[0]}\n";
		}else{
			$ctgid_coverage{$b[0]}=$b[4];
		}
	}else{next;}
}


close IN;

while(<IDX>){
        chomp;
        my $line=$_;
        my @b=split /\t/,$line;
	if(exists $ctgid_coverage{$b[0]}){
		if($ctgid_coverage{$b[0]}>=$cutoff){
			print OUT "$line\n";
		}else{
			print OUT "$b[0]\t$b[1]\t0\t$b[3]\n";}
	}else{
		print OUT "$b[0]\t$b[1]\t0\t$b[3]\n";
	}
}	

close IDX;
close OUT;




