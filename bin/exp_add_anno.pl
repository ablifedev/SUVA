#!/usr/bin/perl -w
#
# Copyright (c)   AB_Life 2011
# Writer:         xuxiong <xuxiong19880610@163.com>
# Program Date:   2011.
# Modifier:       xuxiong <xuxiong19880610@163.com>
# Last Modified:  2011.
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions(\%opts,"exp=s","anno=s","column=s","force=s","annonum=s","rc=s","r","c","u","t","o=s","h" );

if(!defined($opts{exp}) || !defined($opts{anno}))
{
	print <<"	Usage End.";

	Description:This programme is used for ~

		Version: $ver

	Usage:perl $0

		-exp               exp file         must be given

		-anno              anno file                 must be given

		-o                 outfile                   option,default is [exp].anno

		-column            column num for annotation                  option,default is 0;

		-force             ignore line if no anno                     option,default is 0;

		-annonum           "auto",when no anno,use "--" replace                  option,default is 0;

		-rc                for no anno,use this value to replace (when -force is 0)                  option,default is "";

		-r                 replace the exp file                       option,default is no replace;

		-c                 check if exp file had be anno,check "Symbol" in title                       option,default is no check;

		-u                 check if exp have multi anno,only print first one                       option,default is no check;

		-t                 check if exp have title                       option,default is auto check(#,Gene...);





	Usage End.

	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################

my $current_dir = `pwd`;chomp($current_dir);

my $exp_file=$opts{exp};
my $out_file=$exp_file.".anno";
$out_file=$opts{o} if defined($opts{o});
my $anno_file=$opts{anno};
my $column=0;
$column=$opts{column} if defined($opts{column});
my $force=0;
$force=$opts{force} if defined($opts{force});
my $annonum=0;
$annonum=$opts{annonum} if defined($opts{annonum});
my $rc="";
$rc=$opts{rc} if defined($opts{rc});
print $column,"\n";


open(EXP,"$exp_file")||die "Can't open $exp_file\n";
while(<EXP>){
	chomp;
	s/\r//g;
	my @line=split(/\t/);
	# Gene 	 Chromosome 	logConc	logFC	p.value	Tair_Gene
	if(/^Gene\s*/ or /^\#/ or /^X\./ or /^\+/ or /^Chr\s*/){
		if($opts{c}){
			if(/\tSymbol/i){
				die "your exp file title include Symbol,when -c is used,anno is stopped.\n";
			}
		}
		next;
	}
}
close EXP;


open(OUT,">".$out_file)||die $!;

my %anno = ();
open(ANNO,"$anno_file")||die "Can't open $anno_file\n";
while(<ANNO>){
	chomp;
	my @line=split(/\t/);
	my $key = shift @line;
	# $anno{$key}="\t".join("\t",@line);
	push @{$anno{$key}},"\t".join("\t",@line);

}
close ANNO;
open(ANNO,"$anno_file")||die "Can't open $anno_file\n";
my $title = <ANNO>;
chomp($title);
$title=~s/^\S+\s+//;
my @annolist = split(/\t/,$title);
$annonum = scalar(@annolist) if (defined($opts{annonum}) and $opts{annonum} eq "auto");
close ANNO;
# print $title,"\n";
#读取deg
# print $sample1,"\n";
my $titleflag = 1;
open(EXP,"$exp_file")||die "Can't open $exp_file\n";
while(<EXP>){
	chomp;
	next if(/^\s+$/);
	s/\r//g;
	my @line=split(/\t/);
	# Gene 	 Chromosome 	logConc	logFC	p.value	Tair_Gene
	if($titleflag==1){
		$titleflag=0;
		if(/^Gene\s+/i or /^X.ID/i or /^miRNA\s+/i or /^\#/ or /^X\./ or /^\+/ or /^CircRNA\s+/ or /^consensus.cluster/ or /^Chr\s+/i or $opts{t}){
			print OUT $_,"\t",$title,"\n";
			next;
		}
	}

	my $tmp=$line[$column];
	if (not defined($anno{$tmp})){
		if($force==0){
			my $str="";
			my $i = 0;
			while($i<$annonum){
				$str.="\t$rc";
				$i++;
			}
			push @{$anno{$tmp}},$str;
		}else{
			next;
		}
	}
	foreach my $a (@{$anno{$tmp}}){
		print OUT $_;
		print OUT $a,"\n";
		last if($opts{u});
	}
}
close EXP;


close OUT;

if($opts{r}){
	`mv $exp_file $exp_file\.beforeAnno`;
	`mv $out_file $exp_file`;
}

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}



