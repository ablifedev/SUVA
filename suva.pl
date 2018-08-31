#!/usr/bin/perl -w
#
# Copyright (c)   ABLife 2017
# Writer:         Chao Cheng <chaocheng@ablife.cc>
# Program Date:   2017.11.17
my $ver = "1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
$Bin = "$Bin/bin";
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# Designed for splice-site usage variation analysising（suva）


my %opts;
GetOptions( \%opts, "c=s","g=s","a=s","b=s","j=s","t=i","o=s","h" );
if ( !defined( $opts{c} ) ) {
	print <<"	Usage End.";
	Description:This program is used for splice-site usage variation analysising

		Version: $ver

	Usage:perl $0

		-c               config file                               must be given
		-g               genome annotation gff file                default is "NONE",for annotate AS to genes
		-a               gene annotation                           default is "NONE",after annotating AS to genes, annotate genes
		-b               bam file name                             default is "accepted_hits.uniq.bam",which is uniq mapped reads bam file
		-j               splice junction file name                 default is "SJ.out.tab",which is junction file from STAR
		-t               cpu number                                default is 8
		-o               outdir                                    default is current dir

	Usage End.
	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################
my $current_dir = `pwd`;
chomp($current_dir);

## read sample config
my $configfile = $opts{c};
$configfile = "$current_dir/$configfile" if ( $configfile !~ /^\/|\~/ );
my %config = ();  # ConfigOption
&Get_config_option($configfile,\%config);


## read arguments
my $outdir = $opts{o} || $current_dir;
`mkdir -p $outdir` if ( !-d $outdir );
$outdir = "$current_dir/$outdir" if ( $outdir !~ /^\/|\~/ );
chdir($outdir);
`mkdir -p $outdir/outdir`;
# `mkdir -p $outdir/report`;
`mkdir -p $outdir/outdir/NIR`;
`mkdir -p $outdir/outdir/IR`;
`mkdir -p $outdir/outdir/ALL`;

my $bamname = $opts{b} || "accepted_hits.uniq.bam";
my $sjname = $opts{j} || "SJ.out.tab";

my $thread = $opts{t} || 8;

my $gff = $opts{g} || "NONE";
my $geneanno = $opts{a} || "NONE";

if ($gff ne "NONE"){
	$gff = &AbsolutePath("file",$gff);
}
if ($geneanno ne "NONE"){
	$geneanno = &AbsolutePath("file",$geneanno);
}

### NIR analysis
## make cmds
chdir("$outdir/outdir/NIR");
my $sh = "$outdir/outdir/NIR/run.sh";
&Make_nir_cmd($sh);


### IR analysis
## make cmds
chdir("$outdir/outdir/IR");
my $sh2 = "$outdir/outdir/IR/run.sh";
&Make_ir_cmd($sh2);

### combine AS and differential splicing analysis
## make cmds
chdir("$outdir/outdir/ALL");
my $sh3 = "$outdir/outdir/ALL/run.sh";
&Make_all_cmd($sh3);

# run cmds
my $runsh = "$outdir/outdir/runsuva.sh";
open OUT,">$runsh";
print OUT "cd $outdir/outdir/NIR/ && sh run.sh &\n";
print OUT "cd $outdir/outdir/IR/ && sh run.sh &\n";
print OUT "wait\n";
print OUT "cd $outdir/outdir/ALL/ && sh run.sh";
close OUT;
`sh $runsh`;

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";

##################################

sub Get_config_option {
	my ($file,$hash_ref) = @_;
	#print $file,"\n";
	#sample	control	test	control	test	control_star_outdir	test_star_outdir
	# CSCC_F0209	CSCC_F0209_N_E	CSCC_F0209_T_E	/data12/dev/chengc/CS2017/basic/astarget/ALL/basic/result/mapping/CSCC_F0209_N_E_mapping	/data12/dev/chengc/CS2017/basic/astarget/ALL/basic/result/mapping/CSCC_F0209_T_E_L_mapping
	open (CONFIG,"<$file") or die "Can't open $file!\n";
	while (<CONFIG>) {
		chomp;
		next if(/^\#/);
		next if(/^\s*$/);
		s/\r//g;		#chomp the \r character
		my @line = split(/\t/);
		$hash_ref->{$line[0]}->{"tn"} = $line[1];  # test name
		$hash_ref->{$line[0]}->{"cn"} = $line[2];  # control name
		$hash_ref->{$line[0]}->{"td"} = $line[3];  # test STAR outdir
		$hash_ref->{$line[0]}->{"cd"} = $line[4];  # control STAR outdir
	}
	close(CONFIG);
}


sub Make_nir_cmd {
	my ($shfile) = @_;
	my $cmd = "";
	open OUT,">$shfile";

	# print OUT "cd $outdir/outdir/NIR\n";

	print OUT "## 1_combineSJ_and_cluster\n";
	print OUT "mkdir -p $outdir/outdir/NIR/1_combineSJ_and_cluster && cd $outdir/outdir/NIR/1_combineSJ_and_cluster\n";

	foreach my $sample (keys %config){
		$cmd = "ln -s -f ".$config{$sample}{"cd"}."/$sjname ".$config{$sample}{"cn"}."_sj.txt";
		print OUT $cmd,"\n";
		$cmd = "ln -s -f ".$config{$sample}{"td"}."/$sjname ".$config{$sample}{"tn"}."_sj.txt";
		print OUT $cmd,"\n";
	}
	close OUT;

	`cat $Bin/make_cluster.sh >> $shfile`;


	open OUT,">>$shfile";

	# $cmd = 'less -S sj_clu_all_anno | perl -ne \'chomp;@line=split(/\t/);$str="+";if($line[4] eq "2"){$str="-";}print "$line[1]\t$line[2]\t$line[3]\t$line[0]\t0\t$str\n";\' > sj_clu_all_anno.bed';
	# print OUT $cmd,"\n";
	# $cmd = "python3 /users/ablife/RNA-Seq/Pipeline/AS_Pipeline/suva/bed_to_gene_NIR.py -g /data0/Genome/human/gencode_v23/v26_test/genome/gencode.v26.basic.annotation.gff3 -b sj_clu_all_anno.bed -o sj_clu_all_anno_geneinfo.txt ";
	# print OUT $cmd,"\n";




	print OUT "\n\n\n## 2_quantification\n";
	print OUT "mkdir -p $outdir/outdir/NIR/2_quantification && cd $outdir/outdir/NIR/2_quantification\n";

	$cmd = 'ls ../1_combineSJ_and_cluster/*_sj.txt | cat | perl -ne \'BEGIN{open OUT,">cal_sj_ratio.sh";}chomp;$nn=$_;print OUT "perl '.$Bin.'/cal_sj_ratio.pl $nn &\n";$i++;if($i%10==0){print OUT "wait\n";}END{print OUT "wait\n";}\'';
	$cmd.= "\n";
	$cmd.= "sh cal_sj_ratio.sh";
	print OUT "\n\n",$cmd,"\n";

	$cmd = "perl $Bin/pipeline_stat_sj.pl -i ./ -p _sjclu -m _sjclu -c 12 -o all_sj1_ratio.txt";
	print OUT "\n\n",$cmd,"\n";

	$cmd = 'cat '.$configfile.' | perl -ne \'BEGIN{open OUT,">combine_sj.sh";}chomp;next if(/^#/);s/\r//g;@tt=split;print OUT "perl '.$Bin.'/combine_sj.pl $tt[0] $tt[1] $tt[2] &\n";$i++;if($i%10==0){print OUT "wait\n";}END{print OUT "wait\n";}\'';
	$cmd.= "\n";
	$cmd.= "sh combine_sj.sh";
	print OUT "\n\n",$cmd,"\n";

	print OUT "\n\n\n## 3_differential_splicing\n";
	print OUT "mkdir -p $outdir/outdir/NIR/3_differential_splicing && cd $outdir/outdir/NIR/3_differential_splicing\n";

	$cmd = 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;print "R --slave < '.$Bin.'/sj_diff_fisher.r --args ../2_quantification/$tt[0]\_combine_sj_ratio $tt[0]\_combine_sj_diff ./ \&\n";$i++;if($i%10==0){print OUT "wait\n";}END{print OUT "wait\n";}\' >> diff.sh';
	$cmd.= "\n";
	# $cmd.= 'echo "wait" >> diff.sh';
	# $cmd.= "\n";
	$cmd.= "sh diff.sh";
	print OUT "\n\n",$cmd,"\n";

	# $cmd = 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;$f=$_;open OUT,">$tt[0]\_combine_sj_diff_filter_chip";open IN,"$tt[0]\_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-1] eq "NA";if($line[-1]<=0.05 && ($line[-4]>=0.15||$line[-4]<=-0.15)){print OUT $_,"\n";}}close IN;\'';
	# $cmd.= "\n";
	# print OUT $cmd,"\n";

	# $cmd = 'cat *_diff_filter_chip | awk \'$NF<=0.05\' | cut -f 1 | sort | uniq > all_diff_clus && sed -i "1iCluster" all_diff_clus';
	# $cmd.= "\n";
	# $cmd.= "perl $Bin/exp_add_anno.pl -exp all_diff_clus -anno ../1_combineSJ_and_cluster/all_sj1_ratio.txt -o all_diff_clus_sjratio -t &";
	# $cmd.= "\n";
	# $cmd.= "perl $Bin/pipeline_stat_sj.pl -i ../1_combineSJ_and_cluster/ -p _combine_sj_ratio_all -m _combine_sj_ratio_all -c 18 -o all_sj1_diff.txt && perl $Bin/exp_add_anno.pl -exp all_diff_clus -anno all_sj1_diff.txt -o all_diff_clus_sjdiff -t &";
	# $cmd.= "\nwait\n";
	# print OUT $cmd,"\n";

}





sub Make_ir_cmd {
	my ($shfile) = @_;
	my $cmd = "";
	open OUT,">$shfile";
	my $juncs = "";

	print OUT "cd $outdir/outdir/IR\n";
	open OUT2,">$outdir/outdir/IR/starsj2junc.sh";
	foreach my $sample (keys %config){
		my $starsj = $config{$sample}{"cd"}."/$sjname";
		my $junc = $config{$sample}{"cd"}."/SJ.out.tab.junc";
		print OUT2 "sh $Bin/starsj2junc.sh $starsj $junc &\n";
		$juncs .= " $junc";
		$starsj = $config{$sample}{"td"}."/$sjname";
		$junc = $config{$sample}{"td"}."/SJ.out.tab.junc";
		print OUT2 "sh $Bin/starsj2junc.sh $starsj $junc &\n";
		$juncs .= " $junc";
	}
	print OUT2 "wait\n";
	close OUT2;
	print OUT "sh starsj2junc.sh\n";
	print OUT "cat $juncs ".'|awk \'$5>2\'| cut -f 1,2,3,6 | sort | uniq > totalsj '."\n";

	open OUT3,">$outdir/outdir/IR/run_cal_boundaryreads.sh";
	
	my $tn = 0;
	foreach my $sample (keys %config){
		my $dir = $config{$sample}{"cd"};
		# print OUT3 "cd $dir && python $Bin/count_boundary_reads_nopa_sample_star.py -b SJ.out.tab.junc -l $bamname -o SJ.out.tab.junc.br -j 0 -t $outdir/outdir/IR/totalsj &&\n";
		print OUT3 "cd $dir && python $Bin/count_boundary_reads_nopa_sample_star.py -b SJ.out.tab.junc -l $bamname -o SJ.out.tab.junc.br -j 0 -t $outdir/outdir/IR/totalsj &\n";
		$tn++;
		if ($tn % $thread == 0) {
			print OUT3 "wait\n";
		}
		$dir = $config{$sample}{"td"};
		# print OUT3 "cd $dir && python $Bin/count_boundary_reads_nopa_sample_star.py -b SJ.out.tab.junc -l $bamname -o SJ.out.tab.junc.br -j 0 -t $outdir/outdir/IR/totalsj &&\n";
		print OUT3 "cd $dir && python $Bin/count_boundary_reads_nopa_sample_star.py -b SJ.out.tab.junc -l $bamname -o SJ.out.tab.junc.br -j 0 -t $outdir/outdir/IR/totalsj &\n";
		$tn++;
		if ($tn % $thread == 0) {
			print OUT3 "wait\n";
		}
	}
	close OUT3;

	## TODO: for sge
	# print OUT "perl /public/bin/qsub-sge.pl --queue new.q ./run_cal_boundaryreads.sh\n\n\n";
	
	print OUT "sh run_cal_boundaryreads.sh\n\n\n";

	$cmd = "mkdir -p samples_br";
	print OUT $cmd,"\n";
	foreach my $sample (keys %config){
		$cmd = "ln -s -f ".$config{$sample}{"cd"}."/SJ.out.tab.junc.br samples_br/".$config{$sample}{"cn"}."_br.txt";
		print OUT $cmd,"\n";
		$cmd = "ln -s -f ".$config{$sample}{"td"}."/SJ.out.tab.junc.br samples_br/".$config{$sample}{"tn"}."_br.txt";
		print OUT $cmd,"\n";
	}

	print OUT "cat samples_br/*_br.txt | cut -f 1-4 | sort | uniq > br_uniq\n";

	print OUT "cd samples_br\n";

	print OUT 'cat *_br.txt | perl -ne \'chomp;@line=split(/\t/);$line[1]+=1;next if $line[5]<2;next if $line[6]<2;print join("\t",@line),"\n";\' | cut -f 1-4 | sort | uniq | perl -ne \'chomp;$i++;print "cluir$i\t",$_,"\t",$_,"\tir\n";\' > IR_total_cluster'."\n";

	$cmd = 'ls *_br.txt | cat | perl -ne \'BEGIN{open OUT,">cal_ir_ratio.sh";}chomp;$nn=$_;print OUT "perl '.$Bin.'/cal_ir_ratio.pl $nn &\n";END{print OUT "wait\n";}\'';
	$cmd.= "\n";
	$cmd.= "sh cal_ir_ratio.sh";
	print OUT "\n\n",$cmd,"\n";

	print OUT "perl $Bin/pipeline_stat_sj.pl -i ./ -p _sjclu -m _sjclu -c 12 -o all_ir_ratio.txt\n";

	$cmd = 'cat '.$configfile.' | perl -ne \'BEGIN{open OUT,">combine_ir.sh";}chomp;next if(/^#/);s/\r//g;@tt=split;print OUT "perl '.$Bin.'/combine_ir.pl $tt[0] $tt[1] $tt[2] &\n";$i++;if($i%5==0){print OUT "wait\n";}END{print OUT "wait\n";}\'';
	$cmd.= "\n";
	$cmd.= "sh combine_ir.sh";
	print OUT "\n\n",$cmd,"\n";

	#filters
	# print OUT "sh $Bin/filter_ir.sh\n";

	$cmd = "rm -rf diff.sh \n";
	$cmd .= 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;print "touch $tt[0]\_combine_sj_diff && R --slave < '.$Bin.'/sj_diff_fisher.r --args $tt[0]\_combine_sj_ratio $tt[0]\_combine_sj_diff ./ \&\n";$i++;if($i%5==0){print "wait\n";}END{print "wait\n";}\' >> diff.sh';
	$cmd.= "\n";

	$cmd.= "sh diff.sh";
	print OUT "\n\n",$cmd,"\n";

	# diff
	# $cmd = 'cat '.$configfile.' | perl -ne \'chomp;s/\r//g;@tt=split;$f=$_;open OUT,">$tt[0]\_combine_sj_diff_filter_chip";open IN,"$tt[0]\_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-1] eq "NA";if($line[-1]<=0.05 && ($line[-4]>=0.15||$line[-4]<=-0.15)){print OUT $_,"\n";}}close IN;\'';
	# $cmd.= "\n";
	# print OUT $cmd,"\n";

	# $cmd = 'cat *_diff_filter_chip | awk \'$NF<=0.05\' | cut -f 1 | sort | uniq > all_diff_clus';
	# $cmd.= "\n";
	# $cmd.= "perl $Bin/exp_add_anno.pl -exp all_diff_clus -anno all_sj_ratio.txt -o all_diff_clus_sjratio -t ";
	# $cmd.= "\n";
	# $cmd.= "perl $Bin/pipeline_stat_sj.pl -i ./ -p _combine_sj_diff -m _combine_sj_diff -c 18 -o all_sj1_diff.txt && perl $Bin/exp_add_anno.pl -exp all_diff_clus -anno all_sj1_diff.txt -o all_diff_clus_sjdiff -t ";
	# $cmd.= "\n";
	# print OUT $cmd,"\n";

	# $cmd = "Rscript $Bin/deg_heatmap_sj.r -f all_diff_clus_sjratio -n all_diff_clus_sjratio -a /users/chengc/dev2016/CS2017/basic/astarget/ALL/basic/result/exp_analysis/anno.txt &";
	# $cmd.= "\n";
	# $cmd.= "Rscript $Bin/deg_heatmap_sj.r -f all_diff_clus_sjdiff -n all_diff_clus_sjdiff -a /users/chengc/dev2016/CS2017/basic/astarget/ALL/basic/result/exp_analysis/anno_sample.txt &";
	# $cmd.= "\nwait\n";
	# print OUT $cmd,"\n";

	# $cmd = "python /users/chengc/work/ablifedev/Basic/MapRegion/mapping_distribution_analyse.py -d /data0/Genome/human/gencode_v23/v26_test/genome/gencode.v26.basic.annotation.db -b sj_clu_all_anno.bed -m sj_clu_all_anno_geneinfo.txt ";
	# print OUT "\n\n",$cmd,"\n";
}




sub Make_all_cmd {
	my ($shfile) = @_;
	my $cmd = "";
	open OUT,">$shfile";
	my $juncs = "";

	print OUT "cd $outdir/outdir/ALL\n";

	# diff
	$cmd= 'cat '.$configfile.' | perl -ne \'BEGIN{open OUT,">cat.sh";}chomp;next if(/^#/);s/\r//g;@tt=split;$f=$_;print OUT "cat '.$outdir.'/outdir/NIR/3_differential_splicing/$tt[0]\_combine_sj_diff '.$outdir.'/outdir/IR/samples_br/$tt[0]\_combine_sj_diff > $tt[0]\_combine_sj_diff \n";\'';
	$cmd.= "\n";
	$cmd.= "sh cat.sh";
	$cmd.= "\n";
	print OUT $cmd,"\n";

	$cmd= 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;$f=$_;open OUT,">$tt[0]\_combine_sj_diff_filter_chip";open IN,"$tt[0]\_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-1] eq "NA";if($line[-1]<=0.05 && ($line[-4]>=0.15||$line[-4]<=-0.15)){print OUT $_,"\n";}}close IN;\'';
	$cmd.= "\n";
	print OUT $cmd,"\n";

	$cmd= 'cat *_diff_filter_chip | awk \'$NF<=0.05\' | cut -f 1 | sort | uniq > all_diff_clus';
	$cmd.= "\n";

	$cmd.= 'cat '.$outdir.'/outdir/NIR/1_combineSJ_and_cluster/sj_clu_all '.$outdir.'/outdir/IR/samples_br/IR_total_cluster > all_raw_cluster';
	$cmd.= "\n";
	$cmd.= 'cat all_raw_cluster | perl -ne \'chomp;@line=split(/\t/);$s=$line[2];$s = $line[6] if $line[6]<$line[2];$e=$line[3];$e = $line[7] if $line[7]>$line[3];print "$line[1]\t$s\t$e\t$line[0]\t$line[9]\t$line[4]\n";\' > all_raw_cluster.bed ';
	$cmd.= "\n";

	if($gff!~/NONE/i){
		$cmd.= 'python3 '.$Bin."/bed_to_gene_NIR.py -g $gff -b all_raw_cluster.bed -o all_raw_cluster.gene ";
		$cmd.= "\n";
	}
	if($gff!~/NONE/i && $geneanno !~/NONE/i){
		$cmd.= "cat all_raw_cluster.gene  | cut -f 4,7 > cluster.anno && exp_add_anno -exp cluster.anno -anno $geneanno -column 1 -o cluster.anno.final";
		$cmd.= "\n";
	}elsif($gff!~/NONE/i){
		$cmd.= 'cat all_raw_cluster.gene  | cut -f 4,7 > cluster.anno.final';
		$cmd.= "\n";
	}

	if($gff!~/NONE/i){
		$cmd.= 'ls *_combine_sj_diff_filter_chip | cat | perl -ne \'BEGIN{open OUT,">anno.sh";$exps="";}chomp;$nn=$_;$exps.="$nn".",";END{print OUT "perl '.$Bin.'/exp_add_anno_batch.pl -exp $exps -anno cluster.anno.final \n";}\'' if $gff!~/NONE/i;
		$cmd.= "\n";
	}

	$cmd.= "sh anno.sh";
	print OUT "\n\n",$cmd,"\n";


	# $cmd.= "perl $Bin/exp_add_anno.pl -exp all_diff_clus -anno all_sj_ratio.txt -o all_diff_clus_sjratio -t ";
	# $cmd.= "\n";
	# $cmd.= "perl $Bin/pipeline_stat_sj.pl -i ./ -p _combine_sj_diff -m _combine_sj_diff -c 18 -o all_sj1_diff.txt && perl $Bin/exp_add_anno.pl -exp all_diff_clus -anno all_sj1_diff.txt -o all_diff_clus_sjdiff -t ";
	# $cmd.= "\n";
	# print OUT $cmd,"\n";

	# $cmd = "Rscript $Bin/deg_heatmap_sj.r -f all_diff_clus_sjratio -n all_diff_clus_sjratio -a /users/chengc/dev2016/CS2017/basic/astarget/ALL/basic/result/exp_analysis/anno.txt &";
	# $cmd.= "\n";
	# $cmd.= "Rscript $Bin/deg_heatmap_sj.r -f all_diff_clus_sjdiff -n all_diff_clus_sjdiff -a /users/chengc/dev2016/CS2017/basic/astarget/ALL/basic/result/exp_analysis/anno_sample.txt &";
	# $cmd.= "\nwait\n";
	# print OUT $cmd,"\n";

	# $cmd = "python /users/chengc/work/ablifedev/Basic/MapRegion/mapping_distribution_analyse.py -d /data0/Genome/human/gencode_v23/v26_test/genome/gencode.v26.basic.annotation.db -b sj_clu_all_anno.bed -m sj_clu_all_anno_geneinfo.txt ";
	# print OUT "\n\n",$cmd,"\n";
}





sub select_index{
	my ($array_ref,$sample) = @_;
	my $index = 0;
	#print scalar(@{$array_ref}),"\t",$sample,"\n";
	for (my $i=0; $i<scalar(@{$array_ref}); $i++) {
		#print $array_ref->[$i],"\t",$sample,"\n";
		if ($array_ref->[$i] eq $sample ) {		#$array_ref->[$i] =~ /$sample/i
			$index = $i;
		}
	}
	#print $index,"\n";
	return $index;
}

sub AbsolutePath{		# Get the absolute path of the target directory or file
	my ($type,$input) = @_;
	my $return;
	if ($type eq "dir"){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	} elsif($type eq 'file'){
		my $pwd = `pwd`;
		chomp $pwd;
		my $dir=dirname($input);
		my $file=basename($input);
		chdir($dir);
		$return = `pwd`;
		chomp $return;
		$return .="\/".$file;
		chdir($pwd);
	}
	return $return;
}

###############Sub_format_datetime
sub sub_format_datetime {    #Time calculation subroutine
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = @_;
	$wday = $yday = $isdst = 0;
	sprintf( "%4d-%02d-%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $day, $hour, $min, $sec );
}
