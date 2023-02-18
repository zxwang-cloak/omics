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

        This script designed for rna-seq data analysis.
        Usage: $0 config.txt

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
my $thread = $conf{THREAD};

#### write you things ###
open IN,"$sample_f" || die $!;
while(<IN>){
	chomp;
	next if(/^#/);
	my @t = split /\t/;
	my $id = shift @t; 
	
	my ($r1,$r2) = (abs_path($t[0]), abs_path($t[1]));  

	# data filtering
	make_path "$out/$id/00datafilter";
	$mk .= "00fastqc.finished: $r1 $r2\n";
	$mk .= "\tfastqc --threads $thread -f fastq -o $out/$id/00datafilter/ $r1 $r2 > /dev/null 2>/dev/null && touch 00fastqc.finished\n";
	$all .= "00fastqc.finished ";

	$mk .= "00trim.finished: $r1 $r2\n";
	$mk .= "\ttrim_galore -q 25 --phred33 --fastqc --length 36 --stringency 3 --paired -o $out/$id/00datafilter --gzip $r1 $r2 > $out/$id/00datafilter/$id.trim.log 2>$out/$id/00datafilter/$id.trim.err && touch 00trim.finished\n";
	$all .= "00trim.finished ";

	my $clear_r1 = "$out/$id/00datafilter/$id\_R1_val_1.fq.gz";
	my $clear_r2 = "$out/$id/00datafilter/$id\_R2_val_2.fq.gz";
	
	# alignment
	make_path "$out/$id/01alignment";
	$mk .= "01align.finished: 00trim.finished\n";
	$mk .= "\tSTAR --genomeDir $conf{INDEX} --runThreadN $thread --readFilesCommand zcat --readFilesIn $clear_r1 $clear_r2 --outFileNamePrefix $out/$id/01alignment/$id.  --outSAMtype BAM Unsorted --outSAMstrandField intronMotif 1>$out/$id/01alignment/$id.star.log 2>$out/$id/01alignment/$id.star.err && touch 01align.finished\n";
	$all .= "01align.finished ";

	$mk .= "01sort.finished: 01align.finished\n";
	$mk .= "\tsamtools sort $out/$id/01alignment/$id.Aligned.out.bam -o $out/$id/01alignment/$id.bam.sorted 1>$out/$id/01alignment/$id.sort.log 2>$out/$id/01alignment/$id.sort.err && touch 01sort.finished\n";
	$all .= "01sort.finished ";

	$mk .= "01linkbam.finished: 01sort.finished\n"; 
	$mk .= "\tln -s $out/$id/01alignment/$id.Aligned.out.bam $out/$id/01alignment/$id.bam && samtools flagstat $out/$id/01alignment/$id.bam > $out/$id/01alignment/$id.bam.flagstat && touch 01linkbam.finished\n";
	$all .= "01linkbam.finished ";

	#quantification
	make_path "$out/$id/02quantification";
	$mk .= "02quantification.finished: 01linkbam.finished\n";
	$mk .= "\tcufflinks -o $out/$id/02quantification -p $thread -g $conf{GTF} $out/$id/01alignment/$id.bam.sorted  >$out/$id/02quantification/$id.cufflinks.log 2>$out/$id/02quantification/$id.cufflinks.err && touch 02quantification.finished\n";
	$all .= "02quantification.finished ";

	$mk .= "03count.finished: 02quantification.finished\n";
    	$mk .="\thtseq-count -f bam -s no $out/$id/01alignment/$id.bam $conf{GTF} >$out/$id/02quantification/$id.count.txt && touch 03count.finished\n";
    	$all .= "03count.finished ";
        
	# visualization
	make_path "$out/$id/03visual";
	$mk .= "03visual.finished: 02quantification.finished\n";
	$mk .= "\tmakeTagDirectory $out/$id/03visual/homer $out/$id/01alignment/$id.bam -format sam -genome $conf{GENOME} -checkGC > $out/$id/03visual/makeTagDirectory.log 2>$out/$id/03visual/makeTagDirectory.err && makeUCSCfile $out/$id/03visual/homer -style rnaseq -o auto -strand both > $out/$id/03visual/makeUCSCfile.log 2>$out/$id/03visual/makeUCSCfile.err && gunzip -d -c $out/$id/03visual/homer/homer.ucsc.bedGraph.gz | grep -v \"^track\" > $out/$id/03visual/$id.bedGraph && $conf{BIN}/bedGraphToBigWig $out/$id/03visual/$id.bedGraph $conf{CHROMSIZE} $out/$id/03visual/$id.bw > $out/$id/03visual/$id.bw.log 2> $out/$id/03visual/$id.bw.err && touch 03visual.finished\n";
	$all .= "03visual.finished ";
     

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
