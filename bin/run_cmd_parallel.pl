#!/usr/bin/perl -w
# 

my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions(\%opts,"i=s","o=s","cmd=s","p=s","l=s","h");

if(!defined($opts{i}) || !defined($opts{o}) || !defined($opts{cmd}))
{
	print <<"	Usage End.";

	Description:This programme is used for ~
		
		Version: $ver

	Usage:perl $0 -i /data2/huojiao/Assemble/Corn_earworm/Assemble20120304/Assemble_31/transcripts_cluster/transcripts_unigene.fa -o transcripts_format -d /data2/huojiao/Assemble/Corn_earworm/Assemble20120304/Assemble_31/transcripts_cluster/blast_result -l 5000 -e 1e-5 -a 5

		-i      infile                          must be given
		-o      outfile                         must be given
		-cmd    运行的命令，使用单引号引用，输入文件使用<IN>代替，输出文件用<OUT>代替                    must be given 
		-p      processors to use               default is 10
		-l      输入文件的最小单位行数            default is 2

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################
# >Locus_1_Transcript_1/1_Confidence_1.000_Length_559#
my $in=$opts{i};
my $out=$opts{o};
my $cmd=$opts{cmd};
my $p=$opts{p} || 10;
my $minlines=$opts{l} || 2;

print $cmd ,"\n";


$out=&AbsolutePath("file",$out);
$in=&AbsolutePath("file",$in);
my $outdir = &AbsolutePath("dir","./");
print $outdir,"\n";

`rm -rf result`;

my $outcome =`wc -l $in `;
print "$outcome\n";
my $linecount = (split/\s+/,$outcome)[0];
$linecount=int($linecount/$minlines);
my $a = 3;
my $line= (int($linecount/$p)+1)*$minlines;
my $in_name = basename($in);
`split -l $line -d -a $a $in $outdir/$in_name\_subset `;			#split the fa file into subset fa file
my $raw_postfix = basename($in)."_subset";


chdir($outdir);
my @raw;
&load_raw_data($outdir,\@raw,$raw_postfix);
print join("\n",@raw),"\n";
die("not found the raw data\n") if (!@raw);

my $SH="$outdir/parallel_cmd.sh";
open(SH ,">$SH") || die $!;

mkdir("$outdir/result")unless(-d "$outdir/result");
for (my $i=0;$i<@raw ;$i++) {
	mkdir("$outdir/result/$raw[$i]") unless(-d "$outdir/result/$raw[$i]");
	my $tmp_cmd = $cmd ;
	$tmp_cmd =~s/<IN>/$outdir\/$raw[$i]/;
	$tmp_cmd =~s/<OUT>/$outdir\/result\/$raw[$i]\/out_sub/;

	print SH "$tmp_cmd & \n";
	

}
close SH;

`sh $SH`;


my $command="cat $outdir/result/*/out_sub > $out && rm -rf $outdir/*subset*";
# my $command="cat $outdir/result/*/out_sub > $out && rm -rf $outdir/result/ && rm -rf $outdir/*subset*";
# for (my $i=0;$i<@raw ;$i++) {
# 	$command .= " $outdir/result/$raw[$i]/out_sub ";
# }
# mkdir("$outdir/result/$out_name\_entire/") unless(-d "$outdir/result/$out_name\_entire/");
# $command .= ">$out";
print $command,"\n";
system($command);

#############Time_end#############
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";

my $time_used = time() - $start_time;
my $h = $time_used/3600;
my $m = $time_used%3600/60;
my $s = $time_used%3600%60;
printf("\nAll Time used : %d hours\, %d minutes\, %d seconds\n\n",$h,$m,$s);


#######Sub_format_datetime#######
sub sub_format_datetime {    #Time calculation subroutine
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = @_;
	$wday = $yday = $isdst = 0;
	sprintf(
		"%4d-%02d-%02d %02d:%02d:%02d",
		$year + 1900,
		$mon + 1, $day, $hour, $min, $sec
	);
}


sub load_raw_data{
	my ($INDIR,$filenames_ref,$RAW_POSTFIX)=@_;
	#print "($INDIR,$filenames_ref,$RAW_POSTFIX)\n";die;
	opendir(DIR,$INDIR) or die "Can't open $INDIR: $!";
	my $tmp;
	while ($tmp=readdir(DIR)) {
		chomp $tmp;
		next if ($tmp!~/$RAW_POSTFIX/) ;
		push @{$filenames_ref},$tmp;
	}
	@{$filenames_ref} = sort @{$filenames_ref};
	close(DIR);
}
sub AbsolutePath{		#获取指定目录或文件的决定路径
	my ($type,$input) = @_;
	my $return;
	if ($type eq 'dir'){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	}
	elsif($type eq 'file'){
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
