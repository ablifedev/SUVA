file=$1

cat $file | perl -ne 'chomp;@line=split(/\t/);$s=$line[2];$s=$line[6] if $line[6]<$s;$e=$line[3];$e=$line[7] if $line[7]>$e;print $line[1],"\t",$s,"\t",$e,"\t",$line[0],"_$line[10]_$line[11]\n";' > $file.region