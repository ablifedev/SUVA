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

open IN,  "$ARGV[0]\_combine_sj";
open OUT, ">$ARGV[0]\_combine_sj_ratio";

while (<IN>) {
    chomp;
    my @line = split(/\t/);
    my $r    = $line[12] - $line[16];
    if ( ( $line[12] >= 0.15 || $line[16] >= 0.15 ) && ( $line[13] >= 0.15 || $line[17] >= 0.15 ) && ($line[10]+$line[14]>=2) && ($line[11]+$line[15]>=2) ) {
        print OUT $_, "\t", $r, "\n";
    }
}
close IN;
close OUT;
