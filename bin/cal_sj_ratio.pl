chomp;
%sj = ();
$f  = $ARGV[0];
open IN, $f;
while (<IN>) {
    chomp;
    @line     = split(/\t/);
    $s = "+";
    $s = "-" if $line[3] eq "2";
    $key      = "$line[0]\t$line[1]\t$line[2]\t$s";
    $sj{$key} = $line[6];
}
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
    print OUT $_, "\t", $sj{$k1}, "\t", $sj{$k2}, "\t", $ratio1, "\t",
        $ratio2, "\n";
}
