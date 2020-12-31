%n = ();
%t = ();
open OUT, ">$ARGV[0]\_combine_sj";
open IN,  "$ARGV[1]\_sjclu";
while (<IN>) {
    chomp;
    @line    = split(/\t/);
    $value   = "$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[14]";
    $key     = join( "\t", @line[ 0 .. 9 ] );
    $t{$key} = $value;
}
close IN;
open IN, "$ARGV[2]\_sjclu";
while (<IN>) {
    chomp;
    @line    = split(/\t/);
    $value   = "$line[10]\t$line[11]\t$line[12]\t$line[13]\t$line[14]";
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
    my $r    = $line[12] - $line[17];
    print OUT $_, "\t", $r, "\n";
}
close IN;
close OUT;

open IN,  "$ARGV[0]\_combine_sj";
open OUT, ">$ARGV[0]\_combine_sj_ratio";

while (<IN>) {
    chomp;
    my @line = split(/\t/);
    my $r    = $line[12] - $line[17];
    if ( ( $line[12] >= 0.15 || $line[17] >= 0.15 ) && ( $line[13] >= 0.15 || $line[18] >= 0.15 ) && ($line[10]+$line[15]>=2) && ($line[11]+$line[16]>=2) && ($line[10]+$line[11]>=2) && ($line[15]+$line[16]>=2)) {
        print OUT $_, "\t", $r, "\n";
    }
}
close IN;
close OUT;
