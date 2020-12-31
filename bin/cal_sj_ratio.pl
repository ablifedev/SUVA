# GL000008.2      154760  154839  .       2       +

chomp;
%sj = ();
%ss = ();
$f  = $ARGV[0];
$f2  = $ARGV[1];
open IN, $f;
while (<IN>) {
    chomp;
    @line     = split(/\t/);
    $key      = "$line[0]\t$line[1]\t$line[2]\t$line[5]";
    $sj{$key} = $line[4];
}
close IN;

open IN, $f2;
while (<IN>) {
    # chr1    16607   -       110
    chomp;
    @line     = split(/\t/);
    $key      = "$line[0]\t$line[1]\t$line[2]";
    $ss{$key} = $line[3];
}
close IN;

$f =~ /\S+\/(\S+)_sj.txt/;
# $f =~ /(\S+)_sj.txt/;
open OUT, ">$1\_sjclu";
open IN2, "../1_combineSJ_and_cluster/sj_clu_all";
# open IN2, "cluchange";
while (<IN2>) {
    chomp;
    @line = split(/\t/);
    $k1   = "$line[1]\t$line[2]\t$line[3]\t$line[4]";
    $k2   = "$line[5]\t$line[6]\t$line[7]\t$line[8]";
    $sj{$k1} = 0 if not defined( $sj{$k1} );
    $sj{$k2} = 0 if not defined( $sj{$k2} );
    $sum     = $sj{$k1} + $sj{$k2};
    next if $sum <= 5;
    $ratio1 = sprintf( "%.2f", $sj{$k1} / $sum );
    $ratio2 = sprintf( "%.2f", $sj{$k2} / $sum );

    my $t = $line[9];
    my $s = $line[4];
    my $r = 0;
    my $k = "";
    if(($t eq "alt5p" && $s eq "-") or ($t eq "alt3p" && $s eq "+") or ($t eq "olp" && $line[2]>$line[6]) or ($t eq "contain" && $line[2]>$line[6])){
        $k = "$line[1]\t$line[2]\t$line[4]";
        if (not defined $ss{$k} or $ss{$k} == 0){
            $r = 0;
        }
        else{
            $r = sprintf("%.2f",($sj{$k1}+$sj{$k2})/$ss{$k}*100);
        }

    }else{
        $k = "$line[1]\t$line[3]\t$line[4]";
        if (not defined $ss{$k} or $ss{$k} == 0){
            $r = 0;
        }
        else{
            $r = sprintf("%.2f",($sj{$k1}+$sj{$k2})/$ss{$k}*100);
        }
    }

    print OUT $_, "\t", $sj{$k1}, "\t", $sj{$k2}, "\t", $ratio1, "\t",
        $ratio2, "\t$r\n";
}
