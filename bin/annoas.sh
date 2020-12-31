cat all_raw_cluster| perl -ne 'BEGIN{%h=();open IN,"../Events/AS/all_AS_event_sj_anno";while(<IN>){chomp;@line=split;$h{$line[0]}.=",".$line[2];}close IN;}chomp;@line=split;$a="$line[1]:$line[4]:$line[2]:$line[3]";$b="$line[5]:$line[8]:$line[6]:$line[7]";print $_,"\t",$h{$a},"\t",$h{$b},"\n";' > all_raw_cluster_AStype && cat all_raw_cluster_AStype | cut -f 1,11,12 > all_raw_cluster_AStype.anno && cat all_raw_cluster_AStype.anno | perl -ne 'chomp;@line=split(/\t/);if($line[1] eq ""){$line[1]="--"}if($line[2] eq ""){$line[2]="--"}$line[1]=~s/^,//;$line[2]=~s/^,//;print "$line[0]\t$line[1]\t$line[2]\n";' > all_raw_cluster_AStype.anno.final


cat all_raw_cluster| perl -ne 'BEGIN{%h=();open IN,"../Events/Anno/sj.isoform.final";while(<IN>){chomp;@line=split;$h{$line[0]}=$line[1];}close IN;}chomp;@line=split;$a="$line[1]:$line[2]:$line[3]:$line[4]";$b="$line[5]:$line[6]:$line[7]:$line[8]";print $_,"\t",$h{$a},"\t",$h{$b},"\n";' > all_raw_cluster_SJisoform && cat all_raw_cluster_SJisoform | cut -f 1,11,12 > all_raw_cluster_SJisoform.anno && cat all_raw_cluster_SJisoform.anno | perl -ne 'chomp;@line=split(/\t/);$type="known";if($line[1] eq ""){$line[1]="--";$type="novel";}if($line[2] eq ""){$line[2]="--";$type="novel";}$line[1]=~s/^,//;$line[2]=~s/^,//;print "$line[0]\t$line[1];$line[2]\t$type\n";' > all_raw_cluster_SJisoform.anno.final

cat all_raw_cluster| perl -ne 'BEGIN{%hash=();open IN,"../Events/Anno/exon.isoform.final";while(<IN>){chomp;@line=split;$hash{$line[0]}=$line[1];}close IN;}chomp;@line=split;%h=();$h{$line[2]}=1;$h{$line[3]}=1;$h{$line[6]}=1;$h{$line[7]}=1;$s="";foreach $p (sort {$a<=>$b} keys %h) {$k="$line[1]:$p:$line[4]";$hash{$k}="--" if not defined $hash{$k};$s.=";".$hash{$k};}$s=~s/^;//;print $line[0],"\t",$s,"\n";' > all_raw_cluster_exon.anno.final

anno -exp all_raw_cluster_AStype.anno.final -anno all_raw_cluster_SJisoform.anno.final
anno -exp all_raw_cluster_AStype.anno.final.anno -anno all_raw_cluster_exon.anno.final
anno -exp all_raw_cluster_AStype.anno.final.anno.anno -anno cluster.anno.final.tmp
 
mv -f all_raw_cluster_AStype.anno.final.anno.anno.anno cluster.anno.final 
