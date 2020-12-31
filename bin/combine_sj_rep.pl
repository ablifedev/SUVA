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

    

    if ( ( $ratio_t/$ratio_t_n < 0.15 && $ratio_n/$ratio_n_n < 0.15 ) || ( $ratio_t/$ratio_t_n > 0.85 && $ratio_n/$ratio_n_n > 0.85 ) || ($sj1_reads<2) || ($sj2_reads<2) ) {
        # print $ratio_t/$ratio_t_n,"\n";
        next;
    }else{
        print OUT $thisline;
    }
}
close OUT;
close OUT2;
