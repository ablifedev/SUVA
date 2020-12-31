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

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));

open LOG,">log.SUVA.$Time_Start\.txt";

print LOG "\nStart Time :[$Time_Start]\n\n";
print LOG "\nperl $0 ", join(" ",@ARGV),"\n";
###################################


my %opts;
GetOptions( \%opts, "c=s","g=s","a=s","b=s","j=s","t=i","p=f","d=f","o=s","m=s","h" );
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
		-p               pvalue filter                             default is 0.05
		-d               diff ratio filter                         default is 0.15
		-m               mapping software, star or tophat2         default is star
		-o               outdir                                    default is current dir

	Usage End.
	exit;
}



my $current_dir = `pwd`;
chomp($current_dir);

## read sample config
my $configfile = $opts{c};
$configfile = "$current_dir/$configfile" if ( $configfile !~ /^\/|\~/ );
my %config = ();  # ConfigOption
my $rep = 0;
&Get_config_option($configfile,\%config);


## read arguments
my $outdir = $opts{o} || $current_dir;
`mkdir -p $outdir` if ( !-d $outdir );
$outdir = "$current_dir/$outdir" if ( $outdir !~ /^\/|\~/ );
chdir($outdir);
`mkdir -p $outdir/outdir`;
# `mkdir -p $outdir/report`;
`mkdir -p $outdir/outdir/pre`;
`mkdir -p $outdir/outdir/Events`;
`mkdir -p $outdir/outdir/NIR`;
`mkdir -p $outdir/outdir/IR`;
`mkdir -p $outdir/outdir/ALL`;

my $bamname = $opts{b} || "accepted_hits.uniq.bam";
my $sjname = $opts{j} || "suva_sj.txt";

my $thread = $opts{t} || 8;
my $pv = $opts{p} || 0.01;
my $diffratio = $opts{d} || 0.15;

my $gff = $opts{g} || "NONE";
my $geneanno = $opts{a} || "NONE";

my $mm = $opts{m} || "star";

if ($gff ne "NONE"){
	$gff = &AbsolutePath("file",$gff);
}
if ($geneanno ne "NONE"){
	$geneanno = &AbsolutePath("file",$geneanno);
}

### pre analysis
## make cmds
chdir("$outdir/outdir/pre");
my $sh = "$outdir/outdir/pre/run.sh";
&Make_pre_cmd($sh);

### NIR analysis
## make cmds
chdir("$outdir/outdir/NIR");
my $sh1 = "$outdir/outdir/NIR/run.sh";
&Make_nir_cmd($sh1);

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

### AS events analysis
## make cmds
chdir("$outdir/outdir/Events");
my $sh4 = "$outdir/outdir/Events/run.sh";
&Make_events_cmd($sh4);


# run cmds
my $runsh = "$outdir/outdir/runsuva.sh";
open OUT,">$runsh";
print OUT "cd $outdir/outdir/pre/ && sh run.sh \n\n";
print OUT "cd $outdir/outdir/Events/ && sh run.sh &\n";
print OUT "cd $outdir/outdir/NIR/ && sh run.sh &\n";
print OUT "cd $outdir/outdir/IR/ && sh run.sh &\n";
print OUT "wait\n";
print OUT "cd $outdir/outdir/ALL/ && sh run.sh\n";

close OUT;
# `sh $sh`;

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print LOG "\nEnd Time :[$Time_End]\n\n";
close LOG;

##################################

sub Get_config_option {
	my ($file,$hash_ref) = @_;
	#print $file,"\n";
	#sample	control	test	control	test	control_star_outdir	test_star_outdir
	# CSCC_F0209	CSCC_F0209_N_E	CSCC_F0209_T_E	/data12/dev/chengc/CS2017/basic/astarget/ALL/basic/result/mapping/CSCC_F0209_N_E_mapping	/data12/dev/chengc/CS2017/basic/astarget/ALL/basic/result/mapping/CSCC_F0209_T_E_L_mapping
	open (CONFIG,"<$file") or die "Can't open $file!\n";
	my $flag = 0;
	while (<CONFIG>) {
		chomp;
		next if(/^\#/);
		next if(/^\s*$/);
		s/\r//g;		#chomp the \r character


		if(/^\[sample\]/){
			$flag = 1;
			next;
		}elsif(/^\[group\]/){
			$flag = 2;
			next;
		}elsif(/^\[compare\]/){
			$flag = 3;
			next;
		}elsif(/^\[global\]/){
			$flag = 4;
			next;
		}elsif(/^\[/){
			$flag = 0;
			next;
		}

		if($flag == 0){
			next;
		}
		if($flag == 1){
			my @line = split(/\s*=\s*/);
			my @item = split(/:/,$line[1]);
			$hash_ref->{"sample"}->{$line[0]} = $item[0];
			$hash_ref->{"strand"}->{$line[0]} = "strand";
			if(defined($item[1])){
				$hash_ref->{"strand"}->{$line[0]} = $item[1];
			}
		}
		if($flag == 2){
			my @line = split(/\s*=\s*/);
			$hash_ref->{"group"}->{$line[0]} = $line[1];
			$rep = 1;
		}
		if($flag == 3){
			my @line = split(/\s*=\s*/);
			my @item = split(/:/,$line[1]);
			$hash_ref->{"compare"}->{$line[0]}->{"test"} = $item[0];
			$hash_ref->{"compare"}->{$line[0]}->{"ctrl"} = $item[1];
		}

		# my @line = split(/\t/);
		# $hash_ref->{$line[0]}->{"tn"} = $line[1];  # test name
		# $hash_ref->{$line[0]}->{"cn"} = $line[2];  # control name
		# $hash_ref->{$line[0]}->{"td"} = $line[3];  # test STAR outdir
		# $hash_ref->{$line[0]}->{"cd"} = $line[4];  # control STAR outdir


	}
	close(CONFIG);
}


sub Make_pre_cmd {
	my ($shfile) = @_;
	my $cmd = "";
	open OUT,">$shfile";
	my $juncs = "";

	print OUT "mkdir -p $outdir/outdir/pre && cd $outdir/outdir/pre\n";

	foreach my $sample (keys %{$config{"sample"}}){
		# $cmd = "ln -s -f ".$config{$sample}{"cd"}."/$sjname ".$config{$sample}{"cn"}."_sj.txt";
		# print OUT $cmd,"\n";
		# $cmd = "ln -s -f ".$config{$sample}{"td"}."/$sjname ".$config{$sample}{"tn"}."_sj.txt";
		# print OUT $cmd,"\n";
		$cmd = "";
		if($mm eq "star"){
			$cmd = "sh $Bin/starsj2suvasj.sh ".$config{"sample"}{$sample}."/SJ.out.tab ".$sample."_sj.txt &\n";
		}elsif($mm eq "tophat2"){
			$cmd = "sh $Bin/tophatsj2suvasj.sh ".$config{"sample"}{$sample}."/junctions.bed ".$sample."_sj.txt &\n";
		}
		print OUT $cmd;
	}
	print OUT $cmd,"\nwait\n";
	close OUT;

	open OUT,">>$shfile";

	$cmd = "cat *_sj.txt | awk \'\$5>2\' | cut -f 1-3,6 | sort | uniq > sj_uniq";
	$cmd.= "\n";
	print OUT "\n\n",$cmd,"\n";

	open OUT3,">$outdir/outdir/pre/run_cal_boundaryreads.sh";
	
	my $tn = 0;
	foreach my $sample (keys %{$config{"sample"}}){
		my $dir = $config{"sample"}{$sample};
		# print OUT3 "cd $dir && python $Bin/count_boundary_reads_nopa_sample_star.py -b SJ.out.tab.junc -l $bamname -o SJ.out.tab.junc.br -j 0 -t $outdir/outdir/IR/totalsj &&\n";
		if($config{"strand"}{$sample} eq "strand"){
			print OUT3 "cd $outdir/outdir/pre && python3 $Bin/count_boundary_reads_star.py -b ".$sample."_sj.txt -l $dir/$bamname -o ".$sample."_sj.txt.br -j 0 -t $outdir/outdir/pre/sj_uniq &&";
		}else{
			print OUT3 "cd $outdir/outdir/pre && python3 $Bin/count_boundary_reads_star.py -b ".$sample."_sj.txt -l $dir/$bamname -o ".$sample."_sj.txt.br -j 0 -t $outdir/outdir/pre/sj_uniq -u &&";
		}
		
		# print OUT3 "cd $outdir/outdir/pre && ln -s -f $dir/suva_sj.txt.br ".$sample."_sj.txt.br &&";
		my $c = "python3 $Bin/splice_site_sj_stat.py -t $outdir/outdir/pre/sj_uniq -s ".$sample."_sj.txt -o ".$sample."_ss_sjreads.txt && cat ".$sample."_ss_sjreads.txt".' | perl -ne \'BEGIN{%h=();open IN,"'.$sample.'_sj.txt.br";while(<IN>){chomp;@line=split(/\t/);$k="$line[0]\t$line[1]\t$line[3]";$h{$k}=$line[5];$k="$line[0]\t$line[2]\t$line[3]";$h{$k}=$line[6];}close IN;}chomp;@line=split(/\t/);$k="$line[0]\t$line[1]\t$line[2]";$line[3]+=$h{$k};print join("\t",@line),"\n";\' > '.$sample.'_ss_sjreads_final.txt && cat '.$sample.'_sj.txt | perl -ne \'chomp;@line=split(/\t/);$n="$line[0]:$line[1]:$line[2]:$line[5]";print $n,"\t",$line[4],"\n";\' > '.$sample.'_sj_stat.txt';
		print OUT3 $c,"&&\n";
		$tn++;
		# if ($tn % $thread == 0) {
		# 	print OUT3 "wait\n";
		# }
	}

	# print OUT3 "wait\n";
	close OUT3;

	
	# print OUT "\nsh run_cal_boundaryreads.sh\n\n\n";
	## support SGE
	print OUT "\nperl $Bin/qsub-sge.pl ./run_cal_boundaryreads.sh\n\n\n";

}

sub Make_events_cmd {
	my ($shfile) = @_;
	my $cmd = "";
	open OUT,">$shfile";
	my $juncs = "";

	print OUT "mkdir -p $outdir/outdir/Events && cd $outdir/outdir/Events\n";

	$cmd = "sh $Bin/detectASevents.sh $gff $geneanno";
	$cmd.= "\n";
	print OUT "\n\n",$cmd,"\n";
	close OUT;
}


sub Make_nir_cmd {
	my ($shfile) = @_;
	my $cmd = "";
	open OUT,">$shfile";

	# print OUT "cd $outdir/outdir/NIR\n";

	print OUT "## 1_combineSJ_and_cluster\n";
	print OUT "mkdir -p $outdir/outdir/NIR/1_combineSJ_and_cluster && cd $outdir/outdir/NIR/1_combineSJ_and_cluster\n";

	# foreach my $sample (keys %{$config{"sample"}}){
	# 	# $cmd = "ln -s -f ".$config{$sample}{"cd"}."/$sjname ".$config{$sample}{"cn"}."_sj.txt";
	# 	# print OUT $cmd,"\n";
	# 	# $cmd = "ln -s -f ".$config{$sample}{"td"}."/$sjname ".$config{$sample}{"tn"}."_sj.txt";
	# 	# print OUT $cmd,"\n";
	# 	$cmd = "ln -s -f ".$config{"sample"}{$sample}."/$sjname ".$sample."_sj.txt";
	# 	print OUT $cmd,"\n";
	# 	$cmd = "ln -s -f ".$config{"sample"}{$sample}."/suva_ss_sjreads_final.txt ".$sample."_ss_sjreads.txt";
	# 	print OUT $cmd,"\n";
	# }
	# close OUT;

	# `cat $Bin/make_cluster.sh >> $shfile`;


	# open OUT,">>$shfile";

	# $cmd = "cat *_sj.txt | awk \'\$5>2\' | cut -f 1-3,6 | sort | uniq > sj_uniq";
	# $cmd.= "\n";
	$cmd.= "\npython3.6 $Bin/cluster_overlap_sj.py -t $outdir/outdir/pre/sj_uniq -o sj_clu_all\n";
	print OUT "\n\n",$cmd,"\n";

	# $cmd = 'less -S sj_clu_all_anno | perl -ne \'chomp;@line=split(/\t/);$str="+";if($line[4] eq "2"){$str="-";}print "$line[1]\t$line[2]\t$line[3]\t$line[0]\t0\t$str\n";\' > sj_clu_all_anno.bed';
	# print OUT $cmd,"\n";
	# $cmd = "python3 /users/ablife/RNA-Seq/Pipeline/AS_Pipeline/suva/bed_to_gene_NIR.py -g /data0/Genome/human/gencode_v23/v26_test/genome/gencode.v26.basic.annotation.gff3 -b sj_clu_all_anno.bed -o sj_clu_all_anno_geneinfo.txt ";
	# print OUT $cmd,"\n";




	print OUT "\n\n\n## 2_quantification\n";
	print OUT "mkdir -p $outdir/outdir/NIR/2_quantification && cd $outdir/outdir/NIR/2_quantification\n";

	$cmd = 'ls ../../pre/*_sj.txt | cat | perl -ne \'BEGIN{open OUT,">cal_sj_ratio.sh";}chomp;$nn=$_;$ns=$nn;$ns=~s/_sj.txt/_ss_sjreads_final.txt/;print OUT "perl '.$Bin.'/cal_sj_ratio.pl $nn $ns &\n";$i++;if($i%10==0){print OUT "wait\n";}END{print OUT "wait\n";}\'';
	$cmd.= "\n";
	$cmd.= "\nsh cal_sj_ratio.sh\n";
	print OUT "\n\n",$cmd,"\n";

	# $cmd = "perl $Bin/pipeline_stat_sj.pl -i ./ -p _sjclu -m _sjclu -c 12 -o all_sj1_ratio.txt";
	# print OUT "\n\n",$cmd,"\n";

	`mkdir -p $outdir/outdir/NIR/2_quantification/`;
	open CMD,">$outdir/outdir/NIR/2_quantification/combine_sj.sh";
	my $i=0;
	if ($rep==1){
		foreach my $compare (keys %{$config{"compare"}}){
			my $g = $config{"compare"}{$compare}{"test"};
			my $t = $config{"group"}{$g};
			$g = $config{"compare"}{$compare}{"ctrl"};
			my $c = $config{"group"}{$g};
			print CMD "perl ".$Bin."/combine_sj_rep.pl ".$compare." ".$t." ".$c." \&\n";
			$i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}else{
		foreach my $compare (keys %{$config{"compare"}}){
			print CMD "perl ".$Bin."/combine_sj.pl ".$compare." ".$config{"compare"}{$compare}{"test"}." ".$config{"compare"}{$compare}{"ctrl"}." \&\n";
			$i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}

	print CMD "wait\n";
	close CMD;
	# $cmd = 'cat '.$configfile.' | perl -ne \'BEGIN{open OUT,">combine_sj.sh";$flag}chomp;next if(/^#/);s/\r//g;@tt=split;print OUT "perl '.$Bin.'/combine_sj.pl $tt[0] $tt[1] $tt[2] &\n";$i++;if($i%10==0){print OUT "wait\n";}END{print OUT "wait\n";}\'';
	# $cmd.= "\n";
	$cmd= "\nsh combine_sj.sh\n";
	print OUT "\n\n",$cmd,"\n";

	$cmd = "\nperl $Bin/pipeline_stat_exp.pl -i ./ -p _sjclu -m _sjclu -c 12 -o all_NIR_ratio.txt\n";
	print OUT "\n\n",$cmd,"\n";

	print OUT "\n\n\n## 3_differential_splicing\n";
	print OUT "mkdir -p $outdir/outdir/NIR/3_differential_splicing && cd $outdir/outdir/NIR/3_differential_splicing\n";

	`mkdir -p $outdir/outdir/NIR/3_differential_splicing/`;
	open CMD,">$outdir/outdir/NIR/3_differential_splicing/diff.sh";
	my $ii=0;
	if ($rep==1){
		foreach my $compare (keys %{$config{"compare"}}){
			my $g = $config{"compare"}{$compare}{"test"};
			my $t = $config{"group"}{$g};
			my @l = split(",",$t);
			my $tn = 13;
			my $tn_str = "13";
			for(my $i=1;$i<@l;$i++){$tn+=5;$tn_str.=",$tn";}
			$g = $config{"compare"}{$compare}{"ctrl"};
			my $c = $config{"group"}{$g};
			@l = split(",",$c);
			my $cn = $tn+5;
			my $cn_str = "$cn";
			for(my $i=1;$i<@l;$i++){$cn+=5;$cn_str.=",$cn";}
			print CMD "perl $Bin/run_cmd_parallel.pl -i ../2_quantification/".$compare."_combine_sj_ratio -o ".$compare."_combine_sj_diff -cmd".' \'R --slave < '.$Bin.'/sj_diff_ttest.r --args <IN> '.$tn_str.' '.$cn_str.' <OUT> ./\' -p '.$thread.' '," && rm -rf result \n";
			# $i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}else{
		foreach my $compare (keys %{$config{"compare"}}){
			print CMD "R --slave < ".$Bin."/sj_diff_fisher.r --args ../2_quantification/".$compare."\_combine_sj_ratio ".$compare."\_combine_sj_diff ./ \&\n";
			# print CMD "perl $Bin/run_cmd_parallel.pl -i ../2_quantification/".$compare."_combine_sj_ratio -o ".$compare."_combine_sj_diff -cmd".' \'R --slave < '.$Bin.'/sj_diff_fisher.r --args <IN> <OUT> ./\' -p '.$thread.' '," && rm -rf result \n";
			$ii++;if($ii%$thread==0){print CMD "wait\n";}
		}
	}
	print CMD "wait\n";
	close CMD;
	# $cmd = 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;print "R --slave < '.$Bin.'/sj_diff_fisher.r --args ../2_quantification/$tt[0]\_combine_sj_ratio $tt[0]\_combine_sj_diff ./ \&\n";$i++;if($i%10==0){print OUT "wait\n";}END{print OUT "wait\n";}\' >> diff.sh';
	# $cmd.= "\n";
	# $cmd.= 'echo "wait" >> diff.sh';
	# $cmd.= "\n";
	$cmd= "\nsh diff.sh\n";
	print OUT "\n\n",$cmd,"\n";

	# $cmd = 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;$f=$_;open OUT,">$tt[0]\_combine_sj_diff_filter";open IN,"$tt[0]\_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-1] eq "NA";if($line[-1]<=0.05 && ($line[-4]>=0.15||$line[-4]<=-0.15)){print OUT $_,"\n";}}close IN;\'';
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


	# open OUT3,">$outdir/outdir/IR/run_cal_boundaryreads.sh";
	
	# my $tn = 0;
	# foreach my $sample (keys %{$config{"sample"}}){
	# 	my $dir = $config{"sample"}{$sample};
	# 	# print OUT3 "cd $dir && python $Bin/count_boundary_reads_nopa_sample_star.py -b SJ.out.tab.junc -l $bamname -o SJ.out.tab.junc.br -j 0 -t $outdir/outdir/IR/totalsj &&\n";
	# 	print OUT3 "cd $dir && python $Bin/count_boundary_reads_star.py -b $sjname -l $bamname -o $sjname\.br -j 0 -t $outdir/outdir/IR/sj_uniq &\n";
	# 	$tn++;
	# 	if ($tn % $thread == 0) {
	# 		print OUT3 "wait\n";
	# 	}
	# }


	# print OUT3 "wait\n";
	# close OUT3;

	## TODO: for sge
	# print OUT "perl /public/bin/qsub-sge.pl --queue new.q ./run_cal_boundaryreads.sh\n\n\n";
	
	# print OUT "\nsh run_cal_boundaryreads.sh\n\n\n";

	# $cmd = "mkdir -p samples_br";
	# print OUT $cmd,"\n";
	# foreach my $sample (keys %{$config{"sample"}}){
	# 	$cmd = "ln -s -f ".$config{"sample"}{$sample}."/$sjname\.br samples_br/".$sample."_br.txt";
	# 	print OUT $cmd,"\n";
	# 	$cmd = "ln -s -f ".$config{"sample"}{$sample}."/suva_ss_sjreads_final.txt ".$sample."_ss_sjreads.txt";
	# 	print OUT $cmd,"\n";
	# }


	# print OUT "cat $outdir/outdir/pre/*_sj.txt.br | cut -f 1-4 | sort | uniq > br_uniq\n";

	# print OUT "cd samples_br\n";

	print OUT "cat $outdir/outdir/pre/*_sj.txt.br".' | perl -ne \'chomp;@line=split(/\t/);next if $line[5]<5;next if $line[6]<5;print join("\t",@line),"\n";\' | cut -f 1-4 | sort | uniq | perl -ne \'chomp;$i++;print "cluir$i\t",$_,"\t",$_,"\tir\n";\' > IR_total_cluster'."\n";

	$cmd = "ls $outdir/outdir/pre/*_sj.txt.br".' | cat | perl -ne \'BEGIN{open OUT,">cal_ir_ratio.sh";}chomp;$nn=$_;$ns=$nn;$ns=~s/_sj.txt.br/_ss_sjreads_final.txt/;print OUT "perl '.$Bin.'/cal_ir_ratio.pl $nn $ns &\n";END{print OUT "wait\n";}\'';
	$cmd.= "\n";
	$cmd.= "\nsh cal_ir_ratio.sh\n";
	print OUT "\n\n",$cmd,"\n";

	print OUT "perl $Bin/pipeline_stat_sj.pl -i ./ -p _sjclu -m _sjclu -c 12 -o all_ir_ratio.txt\n";

	# `mkdir -p $outdir/outdir/IR/samples_br/`;
	open CMD,">$outdir/outdir/IR/combine_ir.sh";
	my $i=0;
	if ($rep==1){
		foreach my $compare (keys %{$config{"compare"}}){
			my $g = $config{"compare"}{$compare}{"test"};
			my $t = $config{"group"}{$g};
			$g = $config{"compare"}{$compare}{"ctrl"};
			my $c = $config{"group"}{$g};
			print CMD "perl ".$Bin."/combine_ir_rep.pl ".$compare." ".$t." ".$c." \&\n";
			$i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}else{
		foreach my $compare (keys %{$config{"compare"}}){
			print CMD "perl ".$Bin."/combine_ir.pl ".$compare." ".$config{"compare"}{$compare}{"test"}." ".$config{"compare"}{$compare}{"ctrl"}." \&\n";
			$i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}
	# foreach my $compare (keys %{$config{"compare"}}){
	# 	print CMD "perl ".$Bin."/combine_ir.pl ".$compare." ".$config{"compare"}{$compare}{"test"}." ".$config{"compare"}{$compare}{"ctrl"}." \&\n";
	# 	$i++;if($i%$thread==0){print CMD "wait\n";}
	# }
	print CMD "wait\n";
	close CMD;
	# $cmd = 'cat '.$configfile.' | perl -ne \'BEGIN{open OUT,">combine_ir.sh";}chomp;next if(/^#/);s/\r//g;@tt=split;print OUT "perl '.$Bin.'/combine_ir.pl $tt[0] $tt[1] $tt[2] &\n";$i++;if($i%5==0){print OUT "wait\n";}END{print OUT "wait\n";}\'';
	# $cmd.= "\n";
	$cmd = "\nsh combine_ir.sh\n";
	print OUT "\n\n",$cmd,"\n";

	#filters
	# print OUT "sh $Bin/filter_ir.sh\n";
	open CMD,">$outdir/outdir/IR/diff.sh";
	my $ii=0;
	if ($rep==1){
		foreach my $compare (keys %{$config{"compare"}}){
			my $g = $config{"compare"}{$compare}{"test"};
			my $t = $config{"group"}{$g};
			my @l = split(",",$t);
			my $tn = 13;
			my $tn_str = "13";
			for(my $i=1;$i<@l;$i++){$tn+=5;$tn_str.=",$tn";}
			$g = $config{"compare"}{$compare}{"ctrl"};
			my $c = $config{"group"}{$g};
			@l = split(",",$c);
			my $cn = $tn+5;
			my $cn_str = "$cn";
			for(my $i=1;$i<@l;$i++){$cn+=5;$cn_str.=",$cn";}
			# print CMD "R --slave < ".$Bin."/sj_diff_ttest.r --args ".$compare."_combine_sj_ratio $tn_str $cn_str ".$compare."\_combine_sj_diff ./ \&\n";

			print CMD "perl $Bin/run_cmd_parallel.pl -i ".$compare."_combine_sj_ratio -o ".$compare."_combine_sj_diff -cmd".' \'R --slave < '.$Bin.'/sj_diff_ttest.r --args <IN> '.$tn_str.' '.$cn_str.' <OUT> ./\' -p '.$thread.' '," && rm -rf result \n";
			# $i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}else{
		foreach my $compare (keys %{$config{"compare"}}){
			print CMD "R --slave < ".$Bin."/sj_diff_fisher.r --args ".$compare."\_combine_sj_ratio ".$compare."\_combine_sj_diff ./ \&\n";
			# print CMD "perl $Bin/run_cmd_parallel.pl -i ".$compare."_combine_sj_ratio -o ".$compare."_combine_sj_diff -cmd".' \'R --slave < '.$Bin.'/sj_diff_fisher.r --args <IN> <OUT> ./\' -p '.$thread.' '," && rm -rf result \n";
			# print CMD "perl $Bin/run_cmd_parallel.pl -i ../2_quantification/".$compare."_combine_sj_ratio -o ".$compare."_combine_sj_diff -cmd".' \'R --slave < '.$Bin.'/sj_diff_fisher.r --args <IN> <OUT> ./\' -p '.$thread.' '," && rm -rf result \n";
			$ii++;if($ii%$thread==0){print CMD "wait\n";}
		}
	}
	# foreach my $compare (keys %{$config{"compare"}}){
	# 	print CMD "touch $compare\_combine_sj_diff && R --slave < ".$Bin."/sj_diff_fisher.r --args $compare\_combine_sj_ratio $compare\_combine_sj_diff ./ \&\n";
	# 	$i++;if($i%$thread==0){print CMD "wait\n";}
	# }
	print CMD "wait\n";
	close CMD;
	# $cmd = "rm -rf diff.sh \n";
	# $cmd .= 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;print "touch $tt[0]\_combine_sj_diff && R --slave < '.$Bin.'/sj_diff_fisher.r --args $tt[0]\_combine_sj_ratio $tt[0]\_combine_sj_diff ./ \&\n";$i++;if($i%5==0){print "wait\n";}END{print "wait\n";}\' >> diff.sh';
	# $cmd.= "\n";

	$cmd= "\nsh diff.sh\n";
	print OUT "\n\n",$cmd,"\n";

	# diff
	# $cmd = 'cat '.$configfile.' | perl -ne \'chomp;s/\r//g;@tt=split;$f=$_;open OUT,">$tt[0]\_combine_sj_diff_filter";open IN,"$tt[0]\_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-1] eq "NA";if($line[-1]<=0.05 && ($line[-4]>=0.15||$line[-4]<=-0.15)){print OUT $_,"\n";}}close IN;\'';
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
	`mkdir -p $outdir/outdir/ALL/`;
	open CMD,">$outdir/outdir/ALL/cat.sh";
	my $i=0;
	foreach my $compare (keys %{$config{"compare"}}){
		print CMD 'cat '.$outdir."/outdir/NIR/3_differential_splicing/".$compare."_combine_sj_diff ".$outdir."/outdir/IR/".$compare."_combine_sj_diff > ".$compare."_combine_sj_diff \&\n";
		$i++;if($i%$thread==0){print CMD "wait\n";}
	}
	print CMD "wait\n";
	close CMD;
	# $cmd= 'cat '.$configfile.' | perl -ne \'BEGIN{open OUT,">cat.sh";}chomp;next if(/^#/);s/\r//g;@tt=split;$f=$_;print OUT "cat '.$outdir.'/outdir/NIR/3_differential_splicing/$tt[0]\_combine_sj_diff '.$outdir.'/outdir/IR/$tt[0]\_combine_sj_diff > $tt[0]\_combine_sj_diff \n";\'';
	# $cmd.= "\n";
	$cmd.= "\nsh cat.sh\n";
	$cmd.= "\n";
	print OUT $cmd,"\n";

	open CMD,">$outdir/outdir/ALL/filter.sh";
	$i=0;
	if ($rep==1){
		foreach my $compare (keys %{$config{"compare"}}){
			print CMD 'perl -e \'open OUT,">'.$compare.'_combine_sj_diff_filter";open IN,"'.$compare.'_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-2] eq "NA";if($line[-3]<='.$pv.' && ($line[-4]>='.$diffratio.'||$line[-4]<=-'.$diffratio.')){print OUT $_,"\n";}}close IN;\''." \&\n";
			$i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}else{
		foreach my $compare (keys %{$config{"compare"}}){
			print CMD 'perl -e \'open OUT,">'.$compare.'_combine_sj_diff_filter";open IN,"'.$compare.'_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-1] eq "NA";if($line[-1]<='.$pv.' && ($line[-4]>='.$diffratio.'||$line[-4]<=-'.$diffratio.')){print OUT $_,"\n";}}close IN;\''." \&\n";
			$i++;if($i%$thread==0){print CMD "wait\n";}
		}
	}
	print CMD "wait\n";
	close CMD;
	# $cmd= 'cat '.$configfile.' | perl -ne \'chomp;next if(/^#/);s/\r//g;@tt=split;$f=$_;open OUT,">$tt[0]\_combine_sj_diff_filter";open IN,"$tt[0]\_combine_sj_diff";while(<IN>){chomp;@line=split(/\t/);next if $line[-1] eq "NA";if($line[-1]<=0.05 && ($line[-4]>=0.15||$line[-4]<=-0.15)){print OUT $_,"\n";}}close IN;\'';
	# $cmd.= "\n";
	$cmd.= "\nsh filter.sh\n";
	$cmd.= "\n";
	print OUT $cmd,"\n";

	if ($rep==1){
		$cmd= 'cat *_diff_filter | cut -f 1 | sort | uniq > all_diff_clus';
		$cmd.= "\n";
	}else{
		$cmd= 'cat *_diff_filter | cut -f 1 | sort | uniq > all_diff_clus';
		$cmd.= "\n";
	}

	$cmd.= 'cat '.$outdir.'/outdir/NIR/1_combineSJ_and_cluster/sj_clu_all '.$outdir.'/outdir/IR/IR_total_cluster > all_raw_cluster';
	$cmd.= "\n";
	$cmd.= 'cat all_raw_cluster | perl -ne \'chomp;@line=split(/\t/);$s=$line[2];$s = $line[6] if $line[6]<$line[2];$e=$line[3];$e = $line[7] if $line[7]>$line[3];print "$line[1]\t$s\t$e\t$line[0]\t$line[9]\t$line[4]\n";\' > all_raw_cluster.bed ';
	$cmd.= "\n";

	if($gff!~/NONE/i){
		$cmd.= 'python3 '.$Bin."/bed_to_gene.py -g $gff -b all_raw_cluster.bed -o all_raw_cluster.gene ";
		$cmd.= "\n";
	}else{
		$cmd.= 'touch cluster.anno.final.tmp';
		$cmd.= "\n";
	}
	if($gff!~/NONE/i && $geneanno !~/NONE/i){
		$cmd.= "cat all_raw_cluster.gene  | cut -f 4,7 > cluster.anno && exp_add_anno -exp cluster.anno -anno $geneanno -column 1 -o cluster.anno.final.tmp";
		$cmd.= "\n";
	}elsif($gff!~/NONE/i){
		$cmd.= 'cat all_raw_cluster.gene  | cut -f 4,7 > cluster.anno.final.tmp';
		$cmd.= "\n";
	}

	$cmd .= "\n\nsh $Bin/annoas.sh \n\n";

	$cmd.= 'ls *_combine_sj_diff | cat | perl -ne \'BEGIN{open OUT,">anno.sh";$exps="";}chomp;$nn=$_;$exps.="$nn".",";END{print OUT "perl '.$Bin.'/exp_add_anno_batch.pl -exp $exps -anno cluster.anno.final \n";}\'';
	$cmd.= "\n";

	$cmd.= 'ls *_combine_sj_diff_filter | cat | perl -ne \'BEGIN{open OUT,">>anno.sh";$exps="";}chomp;$nn=$_;$exps.="$nn".",";END{print OUT "perl '.$Bin.'/exp_add_anno_batch.pl -exp $exps -anno cluster.anno.final \n";}\'';
	$cmd.= "\n";



	$cmd.= 'ls *_combine_sj_diff_filter | cat | perl -ne \'BEGIN{open OUT,">>anno.sh";$exps="";}chomp;$nn=$_;print OUT "sh '.$Bin.'/getregion.sh $nn\n";\'';
	$cmd.= "\n";

	$cmd.= "\nsh anno.sh\n";
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
	my ($file,$sample) = @_;
	my $index = 0;
	open IN,$file;
	my $title = <IN>;
	close IN;
	chomp($title);
	my @line=split(/\t/,$title);
	for(my $i=0;$i<scalar(@line);$i++){
		if($line[$i] eq $sample){
			$index = $i;
			return $index;
		}
	}
	return -1;
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
	sprintf( "%4d-%02d-%02d-%02d:%02d:%02d", $year + 1900, $mon + 1, $day, $hour, $min, $sec );
}
