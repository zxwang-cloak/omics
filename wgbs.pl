#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script designed for wgbs upstream analysis (form .fastq to CpG methylation counts).
        Usage: perl $0 config.txt
	Author: zxwang (zxwanghelloworld\@foxmail.com)

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my %conf;
&load_conf($conf, \%conf);

# start to write makefile
my $all='all: ';
my $mk;

# configure
my $out = abs_path($conf{OUTDIR});
my $sample_f = abs_path($conf{SAMPLE});

#### write you things ###
open IN,"$sample_f" || die $!;
while(<IN>){
        chomp;
        next if(/^#/);
        my @t = split /\t/;
        my $id = shift @t;

        my ($r1,$r2) = (abs_path($t[0]), abs_path($t[1]));
        # fastp QC
        make_path "$out/$id/00datafilter";
        $mk .= "00filter.finished: $r1 $r2\n";
        $mk .= "\tfastp -i $r1 -I $r2 -o $out/$id/00datafilter/$id\_f1_trimmed.fq.gz -O $out/$id/00datafilter/$id\_r2_trimmed.fq.gz --detect_adapter_for_pe --correction --cut_right --cut_right_window_size=4 --cut_right_mean_quality=15 --length_required=36 --trim_front1=3 --trim_front2=3 1>$out/$id/00datafilter/$id.filter.log 2>$out/$id/00datafilter/$id.filter.err && touch 00filter.finished\n";
        $all .= "00filter.finished ";

	my $clear_f1 = "$out/$id/00datafilter/$id\_f1_trimmed.fq.gz";
        my $clear_r2 = "$out/$id/00datafilter/$id\_r2_trimmed.fq.gz";

	# Bismark alignment
	make_path "$out/$id/01alignment";
	$mk .= "01align.finished: 00filter.finished\n";
	$mk .= "\tbismark --genome $conf{GENOME} -1 $clear_f1 -2 $clear_r2 -o $out/$id/01alignment -X 700 --dovetail --parallel 10 1>$out/$id/01alignment/$id.bismark-alignment.log 2>$out/$id/01alignment/$id.bismark-alignment.err && touch 01align.finished\n";
        $all .= "01align.finished ";

	# Bismark deduplication
	make_path "$out/$id/02deduplication";
	$mk .= "02dedup.finished: 01align.finished\n";
	$mk .= "\tdeduplicate_bismark --bam --paired $out/$id/01alignment/$id\_f1_trimmed_bismark_bt2_pe.bam --output_dir $out/$id/02deduplication --parallel 10 1>$out/$id/02deduplication/$id.bismark-deduplication.log 2>$out/$id/02deduplication/$id.bismark-deduplication.err && touch 02dedup.finished\n";
	$all .= "02dedup.finished ";

	# Bismark methylation extractor
	make_path "$out/$id/03methylation";
	$mk .= "03CpG.finished: 02dedup.finished\n";
	$mk .= "\tbismark_methylation_extractor --paired-end --no_overlap --parallel 10 --comprehensive --bedGraph --counts --genome_folder $conf{GENOME} $out/$id/02deduplication/$id\_f1_trimmed_bismark_bt2_pe.deduplicated.bam -o $out/$id/03methylation 1>$out/$id/03methylation/$id.bismark-methylation-extractor.log 2>$out/$id/03methylation/$id.bismark-methylation-extractor.err && touch 03CpG.finished\n";
	$all .= "03CpG.finished ";

	# clear
	$mk .= "clear.finished: 03CpG.finished\n";
	$mk .= "\tmv fastp.html fastp.json $out/$id/00datafilter && rm $clear_f1 $clear_r2 $out/$id/01alignment/$id\_f1_trimmed_bismark_bt2_pe.bam $out/$id/02deduplication/$id\_f1_trimmed_bismark_bt2_pe.deduplicated.bam $out/$id/03methylation/CHH_context_$id\_f1_trimmed_bismark_bt2_pe.deduplicated.txt $out/$id/03methylation/CHG_context_$id\_f1_trimmed_bismark_bt2_pe.deduplicated.txt && gzip $out/$id/03methylation/CpG_context_$id\_f1_trimmed_bismark_bt2_pe.deduplicated.txt && touch clear.finished\n";
	$all .= "clear.finished ";

        make_path abs_path($conf{OUTDIR});
        open OUT, ">$out/$id/makefile";
        print OUT $all, "\n";
        print OUT $mk, "\n";
        close OUT;
        $all= 'all: ';
        $mk = "";
}
close IN;

sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}


