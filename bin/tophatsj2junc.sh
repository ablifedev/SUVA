sjfile=$1
juncfile=$2
# cat $sjfile | perl -ne 'chomp;@line=split(/\t/);$str="+";if($line[3] eq "2"){$str="-";}$line[1]-=1;print "$line[0]\t$line[1]\t$line[2]\t.\t$line[6]\t$str\n";' > $juncfile

cat $sjfile | perl -ne 'chomp;next if ( $_ =~ /^track/i );my ( $chr, $left, $right, $feature, $readsNum, $strand, $thickStart, $thickEnd,$itemRgb, $blockCount, $blockSizes, $blockStarts) = split /\s+/ ;my ( $leftsize, $rightsize ) = ( split /\,/, $blockSizes )[ 0 .. 1 ];$start = $left + $leftsize - 1;$end   = $right - $rightsize + 1;print "$chr\t$start\t$end\t.\t$readsNum\t$strand\t$start\t$end\t0,0,0\t2\t1,1\t0,0\n";' > $juncfile