# chomp;%sj=();$f=$ARGV[0];$f=~/(\w+)_br.txt/;open OUT,">$1\_sjclu";open IN,$f;while(<IN>){chomp;@line=split;next if $line[4]+$line[5]+$line[6]==0;$k="$line[0]:$line[1]:$line[2]:$line[3]";$bravg=($line[5]+$line[6])/2;$r=sprintf("%.2f",$bravg/($bravg+$line[4]));print OUT $k,"\t",$r,"\n";}close IN;close OUT;




chomp;
%sj = ();
%ss = ();
%minbr = ();
$f  = $ARGV[0];
$f2  = $ARGV[1];
open IN, $f;
while (<IN>) {
    chomp;
    @line     = split(/\t/);
    $key      = "br\t$line[0]\t$line[1]\t$line[2]\t$line[3]";
    $key2      = "sj\t$line[0]\t$line[1]\t$line[2]\t$line[3]";
    $bravg=($line[5]+$line[6])/2;
    $minbr{$key} = $line[5];
    $minbr{$key} = $line[6] if $line[6]<$line[5];
    $sj{$key} = $bravg;
    $sj{$key2} = $line[4];
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

$f=~s/(\S+)\///;
$f =~ /(\S+)_sj.txt.br/;
# $f =~ /(\w+)_br.txt/;
open OUT, ">$1\_sjclu";
# chr10	100069869	100075910	-	5	1	0
# open IN2, "../1_combineSJ_and_cluster/sj_clu_all_anno";
open IN2, "IR_total_cluster";
while (<IN2>) {
    chomp;
    @line = split(/\t/);
    # $line[2]-=1;
    # $line[6]-=1;
    $k1   = "br\t$line[1]\t$line[2]\t$line[3]\t$line[4]";
    $k2   = "sj\t$line[1]\t$line[2]\t$line[3]\t$line[4]";
    # $k2   = "sj\t$line[5]\t$line[6]\t$line[7]\t$line[8]";

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
    my $ka = "";
    my $kb = "";
    if(($t eq "alt5p" && $s eq "-") or ($t eq "alt3p" && $s eq "+") or ($t eq "olp" && $line[2]>$line[6]) or ($t eq "contain" && $line[2]>$line[6])){
        $k = "$line[1]\t$line[2]\t$line[4]";
        if (not defined $ss{$k} or $ss{$k} == 0){
            $r = 0;
        }
        else{
            # $r = ($sj{$k1}+$sj{$k2})/$ss{$k}*100;
            $r = sprintf("%.2f",($sj{$k1}+$sj{$k2})/$ss{$k}*100);
        }

    }elsif($t eq "ir"){
        $ka = "$line[1]\t$line[2]\t$line[4]";
        $kb = "$line[1]\t$line[3]\t$line[4]";
        if (not defined $ss{$ka}){$ss{$ka}=0;}
        if (not defined $ss{$kb}){$ss{$kb}=0;}
        if ($ss{$ka} + $ss{$kb} == 0){
            $r = 0;
        }
        else{
            $r = sprintf("%.2f",($sj{$k1}+$sj{$k2})*2/($ss{$ka}+$ss{$kb})*100);
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

    print OUT $_, "\t", $minbr{$k1}, "\t", $sj{$k1}, "\t", $sj{$k2}, "\t", $ratio1, "\t",
        $ratio2, "\t$r\n";
}
