sjfile=$1
juncfile=$2
cat $sjfile | perl -ne 'chomp;@line=split(/\t/);$str="+";if($line[3] eq "2"){$str="-";}$line[1]-=1;print "$line[0]\t$line[1]\t$line[2]\t.\t$line[6]\t$str\n";' > $juncfile
