# %n=();%t=();open OUT,">$ARGV[0]\_combine_sj";print OUT "#chr\tstart\tend\tstrand\tT_sj\tT_left_br\tT_right_br\tN_sj\tN_left_br\tN_right_br\n";open IN,"$ARGV[1]\_br.txt";while(<IN>){chomp;@line=split(/\t/);$value="$line[4]\t$line[5]\t$line[6]";$key="$line[0]\t$line[1]\t$line[2]\t$line[3]";$t{$key}=$value;}close IN;open IN,"$ARGV[2]\_br.txt";while(<IN>){chomp;@line=split(/\t/);$value="$line[4]\t$line[5]\t$line[6]";$key="$line[0]\t$line[1]\t$line[2]\t$line[3]";$n{$key}=$value;}close IN;foreach $k (keys(%t)){next if not defined($n{$k});print OUT $k,"\t",$t{$k},"\t",$n{$k},"\n";}


%n = ();
%t = ();
%minbr_n = ();
%minbr_t = ();
open OUT, ">$ARGV[0]\_combine_sj";
open IN,  "$ARGV[1]\_sjclu";
while (<IN>) {
    chomp;
    @line    = split(/\t/);
    $value   = "$line[11]\t$line[12]\t$line[13]\t$line[14]\t$line[15]";
    $key     = join( "\t", @line[ 0 .. 9 ] );
    $t{$key} = $value;
    $minbr_t{$key} = $line[10];
}
close IN;
open IN, "$ARGV[2]\_sjclu";
while (<IN>) {
    chomp;
    @line    = split(/\t/);
    $value   = "$line[11]\t$line[12]\t$line[13]\t$line[14]\t$line[15]";
    $key     = join( "\t", @line[ 0 .. 9 ] );
    $n{$key} = $value;
    $minbr_n{$key} = $line[10];
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
    my $r    = $line[12] - $line[17];
    print OUT $_, "\t", $r, "\n";
}
close IN;
close OUT;


use FindBin qw($Bin $Script);

`cat ../pre/$ARGV[1]\_sj.txt ../pre/$ARGV[2]\_sj.txt | awk '\$5>=2' | cut -f 1-3,6 > $ARGV[0]\_sj_uniq`;
`python3.6 $Bin/cluster_overlap_sj.py -t $ARGV[0]\_sj_uniq -o $ARGV[0]\_sj_uniq_cluster`;
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
    my $k=join("\t",@line[0..9]);
    next if $hash{$sj}==1;
    my $r    = $line[12] - $line[17];
    if ( ( $line[12] >= 0.15 || $line[17] >= 0.15 ) && ( $line[13] >= 0.15 || $line[18] >= 0.15 ) && ($line[10]+$line[15]>=2) && ($line[11]+$line[16]>=2) && ($line[10]+$line[11]>=2) && ($line[15]+$line[16]>=2) && ( $minbr_t{$k} >= 2 || $minbr_n{$k} >= 2 )  ) {
        print OUT $_, "\t", $r, "\n";
    }
}
close IN;
close OUT;
