# %n=();%t=();open OUT,">$ARGV[0]\_combine_sj";print OUT "#chr\tstart\tend\tstrand\tT_sj\tT_left_br\tT_right_br\tN_sj\tN_left_br\tN_right_br\n";open IN,"$ARGV[1]\_br.txt";while(<IN>){chomp;@line=split(/\t/);$value="$line[4]\t$line[5]\t$line[6]";$key="$line[0]\t$line[1]\t$line[2]\t$line[3]";$t{$key}=$value;}close IN;open IN,"$ARGV[2]\_br.txt";while(<IN>){chomp;@line=split(/\t/);$value="$line[4]\t$line[5]\t$line[6]";$key="$line[0]\t$line[1]\t$line[2]\t$line[3]";$n{$key}=$value;}close IN;foreach $k (keys(%t)){next if not defined($n{$k});print OUT $k,"\t",$t{$k},"\t",$n{$k},"\n";}


%n = ();
%t = ();
open OUT, ">$ARGV[0]\_combine_sj";
open IN,  "$ARGV[1]\_sjclu";
while (<IN>) {
    chomp;
    @line    = split(/\t/);
    $value   = "$line[10]\t$line[11]\t$line[12]\t$line[13]";
    $key     = join( "\t", @line[ 0 .. 9 ] );
    $t{$key} = $value;
}
close IN;
open IN, "$ARGV[2]\_sjclu";
while (<IN>) {
    chomp;
    @line    = split(/\t/);
    $value   = "$line[10]\t$line[11]\t$line[12]\t$line[13]";
    $key     = join( "\t", @line[ 0 .. 9 ] );
    $n{$key} = $value;
}
close IN;
foreach $k ( keys(%t) ) {
    next if not defined( $n{$k} );
    print OUT $k, "\t", $t{$k}, "\t", $n{$k}, "\n";
}
close OUT;

open IN,  "$ARGV[0]\_combine_sj";
open OUT, ">$ARGV[0]\_combine_sj_ratio_all";

while (<IN>) {
    chomp;
    my @line = split(/\t/);
    my $r    = $line[12] - $line[16];
    print OUT $_, "\t", $r, "\n";
}
close IN;
close OUT;


use FindBin qw($Bin $Script);

`cat ../../NIR/1_combineSJ_and_cluster/$ARGV[1]\_sj.txt ../../NIR/1_combineSJ_and_cluster/$ARGV[2]\_sj.txt | awk '\$7>1' | cut -f 1-4 > $ARGV[0]\_sj_uniq`;
`python2.7 $Bin/cluster_overlap_sj.py -t $ARGV[0]\_sj_uniq -o $ARGV[0]\_sj_uniq_cluster`;
`less -S $ARGV[0]\_sj_uniq_cluster | awk '\$10!~/contain/' | awk '{print \$2"\\t"\$3"\\t"\$4"\\t"\$5;print \$6"\\t"\$7"\\t"\$8"\\t"\$9;}' | sort | uniq > $ARGV[0]\_sj_for_filter_ir`;




%hash=();
open IN,"$ARGV[0]\_sj_for_filter_ir";
while(<IN>){chomp;@line=split(/\t/);$hash{$_}=1;}
close IN;



open IN,  "$ARGV[0]\_combine_sj";
open OUT, ">$ARGV[0]\_combine_sj_ratio";

while (<IN>) {
    chomp;
    my @line = split(/\t/);
    my $sj=join("\t",@line[1..4]);
    next if $hash{$sj}==1;
    my $r    = $line[12] - $line[16];
    if ( ( $line[12] >= 0.15 || $line[16] >= 0.15 ) && ( $line[13] >= 0.15 || $line[17] >= 0.15 ) && ($line[10]+$line[14]>=2) && ($line[11]+$line[15]>=2) ) {
        print OUT $_, "\t", $r, "\n";
    }
}
close IN;
close OUT;
