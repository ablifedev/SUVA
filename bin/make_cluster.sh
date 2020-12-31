# make cluster
cat *_sj.txt | awk '$7>2' | awk '$4!=0' | cut -f 1-4 | sort | uniq > sj_uniq
# less -S sj_uniq | sort -k 1,1 -k 2,2n > sj_uniq_sortbystart
# less -S sj_uniq | sort -k 1,1 -k 3,3n > sj_uniq_sortbyend
# less -S sj_uniq_sortbystart | perl -ne 'BEGIN{@sj=();$i=0;$fs="";}chomp;@line=split(/\t/);if($line[1] ne $fs){@sj=();push @sj,$_;}else{foreach $key (@sj){$i++;print "cls",$i,"\t",$key,"\t",$_,"\n";}push @sj,$_;}$fs=$line[1];' > sj_clu_bystart
# less -S sj_uniq_sortbyend | perl -ne 'BEGIN{@sj=();$i=0;$fs="";}chomp;@line=split(/\t/);if($line[2] ne $fs){@sj=();push @sj,$_;}else{foreach $key (@sj){$i++;print "cle",$i,"\t",$key,"\t",$_,"\n";}push @sj,$_;}$fs=$line[2];' > sj_clu_byend
# cat sj_clu_bystart sj_clu_byend > sj_clu_all_tmp
# less -S sj_clu_all_tmp |  perl -ne 'chomp;@line=split;$a=$line[3]-$line[2];$b=$line[7]-$line[6];if($a>$b){$t=$line[2];$line[2]=$line[6];$line[6]=$t;$t=$line[3];$line[3]=$line[7];$line[7]=$t;}print join("\t",@line),"\n";' > sj_clu_all
# less -S sj_clu_all | perl -ne 'chomp;@line=split(/\t/);$alt="";if(/^cls/ && $line[4] eq "1"){$alt="alt3p";}elsif(/^cls/ && $line[4] eq "2"){$alt="alt5p";}elsif(/^cle/ && $line[4] eq "1"){$alt="alt5p";}elsif(/^cle/ && $line[4] eq "2"){$alt="alt3p";}print $_,"\t$alt\n";' > sj_clu_all_anno
python3.6 /users/ablife/RNA-Seq/Pipeline/AS_Pipeline/suva_v4/bin/cluster_overlap_sj.py -t sj_uniq -o sj_clu_all