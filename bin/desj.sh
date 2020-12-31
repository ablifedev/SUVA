perl /users/ablife/ablife-perl/public/pipeline_stat_exp.pl -i ../pre/ -p sj_stat.txt -m sj_stat.txt -c 1 -o sj_exp_reads.txt
cat sj_exp_reads.txt | perl -ne 'chomp;@line=split;shift @line;$i=0;foreach $k(@line){$i=$k if $k>$i;}if($i>2){print $_,"\n";}' > sj_exp_reads_filter.txt

