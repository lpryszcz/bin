#!/usr/bin/perl
#http://www.biostars.org/p/61905/
use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($min_reads) = (10);
GetOptions (
        "min-reads:s" => \$min_reads,
);

while (<>) {
        my $line = $_;
        my @columns = split("\t",$line);
        my $chr = $columns[0];
        my $start = $columns[1];
        my $end = $start + 1;
        my $num_reads = $columns[3];
        my $calls = $columns[4];
        my $id = "mpileup_number_" . $.;
        if($num_reads < $min_reads){ # not enough coverage to have good confidence in the call
                next;
        }
        my $num_ref = 0;
        while ($calls =~ /[,.]/g) { $num_ref++ }
        my $num_var = $num_reads - $num_ref;
        my $varAlleleFreq = ($num_var/$num_reads)*100;
        print("$chr\t$start\t$end\t$id\t$varAlleleFreq\n");
}
