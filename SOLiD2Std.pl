#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = qq(
This script is used for converting SOLiD/ABi FASTQ to Standard(Sanger) FASTQ
Usage :  perl $0 [Option]
Option:  -seq           seq file
         -qual          qual file
         -fastq         FASTQ file
         -o             output file
         -help          Help Information
Example: perl $0 -fastq solid.fastq -o std.fastq

Author : BENM <BinxiaoFeng\@gmail.com>
Version: 1.1
Date   : 2009-10-06
\n);
my ($Seq,$Qual,$Fastq,$Output,$Header,$Help);
my %opts;
GetOptions
(
	\%opts,
	"seq:s"=>\$Seq,
	"qual:s"=>\$Qual,
	"fastq:s"=>\$Fastq,
	"o:s"=>\$Output,
	"header:s"=>\$Header,
	"help"=>\$Help
);

die($usage) if (((!defined $Seq)||(!defined $Qual))&&(!defined $Fastq)&&(!defined $Output)||($Help));

# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64) {
  $conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}


# SOLiD color code
my @code = ([0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]);
my @bases = qw(A C G T);
my %decode = ();
foreach my $i(0..3)
{
	foreach my $j(0..3)
	{
		$decode{$code[$i]->[$j]} -> {$bases[$i]} = $bases[$j];
	}
}

open (OUT,">$Output") || die "Fail to create FASTA file:$Output for writing\n";
if (defined $Fastq)
{
	open (IN,$Fastq) || die "Fail to open FASTQ file:$Fastq for reading\n";
	while (<IN>)
	{
		if (/^@/)
		{
			print OUT $_;
			$_ = <IN>;
			s/\s+$//;
			my $seq = ($_=~/\d+/) ? col2base($_) : $_;
			print OUT "$seq\n";
			$_ = <IN>; print OUT $_; $_ = <IN>;
			s/\s+$//;
			my @t = split;
			my $qual = '';
			if ($_=~/^\d+\s*/)
			{
				my @t=split;
				$qual .= $conv_table[$_+64]; #intfq: SOLiD -> Standard
			}
			else
			{
				$qual .= $_;
			}
			substr($qual,0,1,'') if (length($qual)>length($seq));
			print OUT "$qual\n";
		}
	}
	close IN;
}
elsif ((defined $Seq)&&(defined $Qual))
{
	my ($Title,$Name,$seq,$qual)=("","","","");
	open (IN1,$Seq) || die "Fail to open Seq file:$Seq for reading\n";
	open (IN2,$Qual) || die "Fail to open Quality file:$Qual for reading\n";
	while(<IN2>)
	{
		if ($_=~/^\#\s+Title\:\s+(\S+)/)
		{
			$Title=$1;
			<IN1>;
		}
		else
		{
			if (/^[\@\>](\S+)/)
			{
				substr($qual,0,1,"") if (length($qual)>$seq);
				warn "The length of quality is not match with seqence!\nSeq: $seq\nQual: $qual\n" if (length($qual)!=$seq);
				print OUT "$Name\n$seq\n+\n$qual\n" if ($seq ne "")&&($qual ne "");
				($seq,$qual)=("","");
				$Name="@".$Title.$1;
				$_ = <IN1>;
				if (/^[\@\>]/)
				{
					$_ = <IN1>;
				}
				my $tmp_seq;
				while (($_ ne "")&&($_ !~ /^[\@\>]/)&&(!eof))
				{
					s/\s+$//;
					$tmp_seq.=$_;
					$_=<IN1>;
					chomp $_;
				}
				$tmp_seq .= $_ if ((eof)&&($_ !~ /^[\@\>]/));
				chomp $tmp_seq;
				$seq = (($tmp_seq=~/\d+/)) ? col2base($tmp_seq) : $tmp_seq;
				$_ = <IN2>;
				if ($_=~/^\d+\s*/)
				{
					s/^(\d+)\s*//;
					my @t=split;
					$qual .= $conv_table[$_+64]; #intfq: SOLiD -> Standard
				}
				else
				{
					$qual .= $_; #ASCII: SOLiD -> Standard
				}
			}
			else
			{
				if ($_=~/^\d+\s*/)
				{
					s/^(\d+)\s*//;
					my @t=split;
					$qual .= $conv_table[$_+64]; #intfq: SOLiD -> Standard
				}
				else
				{
					$qual .= $_; 
				}
			}
		}
	}
	substr($_,0,1,'') if (length($qual)>length($seq));
	print OUT "$Name\n$seq\n+\n$qual\n" if ($seq ne "")&&($qual ne "");
	close IN1;
	close IN2;
}
close OUT;

sub col2base
{
	my $col = shift;
	my @colors = split '',$col;
	my $string = $colors[0];
	if($string !~ /[acgt]/i){
		warn "$col has no header base\n";
		return 0;
	}
	my $last_base = $string;
	my $current_base = '';
	for(my $i=1;$i<@colors;$i++)
	{
		if (($last_base=~/N/i)&&($colors[$i]==5))
		{
			$current_base = $bases[int(rand(@bases))];
		}
		else
		{
			$current_base = (exists $decode{$colors[$i]}->{$last_base}) ? $decode{$colors[$i]}->{$last_base} : "N";
		}
		$string .= $current_base;
		$last_base = $current_base;
	}
	substr($string,0,1,"");
	return $string;
}

