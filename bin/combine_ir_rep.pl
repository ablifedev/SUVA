
use FindBin qw($Bin $Script);

my $sjfiles = "";
foreach my $s (split(/,/,$ARGV[1])){
    $sjfiles.=" ../pre/$s\_sj.txt";
}
foreach my $s (split(/,/,$ARGV[2])){
    $sjfiles.=" ../pre/$s\_sj.txt";
}

`cat $sjfiles | awk '\$5>1' | cut -f 1-3,6 | sort | uniq > $ARGV[0]\_sj_uniq`;
`python3.6 $Bin/cluster_overlap_sj.py -t $ARGV[0]\_sj_uniq -o $ARGV[0]\_sj_uniq_cluster`;
`less -S $ARGV[0]\_sj_uniq_cluster | awk '\$10!~/contain/' | awk '{print \$2"\\t"\$3"\\t"\$4"\\t"\$5;print \$6"\\t"\$7"\\t"\$8"\\t"\$9;}' | sort | uniq > $ARGV[0]\_sj_for_filter_ir`;

%sjhash=();
open IN,"$ARGV[0]\_sj_for_filter_ir";
while(<IN>){chomp;@line=split(/\t/);$sjhash{$_}=1;}
close IN;


%n = ();
%t = ();
%hash=();
open OUT, ">$ARGV[0]\_combine_sj_ratio";
open OUT2, ">$ARGV[0]\_combine_sj_ratio_all";

foreach my $s (split(/,/,$ARGV[1])){
    # print $s,"\n";
    open IN,  "$s\_sjclu";
    while (<IN>) {
        chomp;
        @line    = split(/\t/);
        $value   = "$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[14]";
        $key     = join( "\t", @line[ 0 .. 9 ] );
        $hash{$key}{$s} = $value;
    }
    close IN;
}
foreach my $s (split(/,/,$ARGV[2])){
    # print $s,"\n";
    open IN,  "$s\_sjclu";
    while (<IN>) {
        chomp;
        @line    = split(/\t/);
        $value   = "$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[14]";
        $key     = join( "\t", @line[ 0 .. 9 ] );
        $hash{$key}{$s} = $value;
    }
    close IN;
}


foreach my $k ( keys(%hash) ) {
    my $thisline = $k;
    my $ratio_t=0;
    my $ratio_t_n=0;
    my $ratio_n=0;
    my $ratio_n_n=0;
    my $sj1_reads=0;
    my $sj2_reads=0;
    my $ss1_per=0;
    my $ss2_per=0;

    foreach my $s (split(/,/,$ARGV[1])){
        # print $s,"\n";
        if (not defined ($hash{$k}{$s})){
            $hash{$k}{$s} = "NA\tNA\tNA\tNA\tNA";
        }else{
            my @tt = split(/\t/,$hash{$k}{$s});
            $ratio_t += $tt[2];
            $sj1_reads += $tt[0];
            $sj2_reads += $tt[1];
            $ss1_per += $tt[4];
            $ratio_t_n += 1;
        }
        $thisline .= "\t".$hash{$k}{$s};
    }
    foreach my $s (split(/,/,$ARGV[2])){
        # print $s,"\n";
        if (not defined ($hash{$k}{$s})){
            $hash{$k}{$s} = "NA\tNA\tNA\tNA\tNA";
        }else{
            my @tt = split(/\t/,$hash{$k}{$s});
            $ratio_n += $tt[2];
            $sj1_reads += $tt[0];
            $sj2_reads += $tt[1];
            $ss2_per += $tt[4];
            $ratio_n_n += 1;
        }
        $thisline .= "\t".$hash{$k}{$s};
    }
    next if $ratio_t_n==0 || $ratio_n_n==0;
    my $r = $ratio_t/$ratio_t_n - $ratio_n/$ratio_n_n;
    my $p1 = sprintf("%.2f",$ss1_per/$ratio_t_n);
    my $p2 = sprintf("%.2f",$ss2_per/$ratio_n_n);
    $thisline .= "\t$p1\t$p2\t".$r."\n";
    print OUT2 $thisline;
    my @line = split(/\t/,$k);
    my $sj=join("\t",@line[1..4]);
    next if $sjhash{$sj}==1;
    if ( ( $ratio_t/$ratio_t_n < 0.15 && $ratio_n/$ratio_n_n < 0.15 ) || ( $ratio_t/$ratio_t_n > 0.85 && $ratio_n/$ratio_n_n > 0.85 ) || ($sj1_reads<2) || ($sj2_reads<2) ) {
        next;
    }else{
        print OUT $thisline;
    }
}
close OUT;
close OUT2;
