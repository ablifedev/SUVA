gff=$1
anno=$2
bin=$3


mkdir -p Anno
mkdir -p SJ
mkdir -p AS
mkdir -p KnownAS

cat ../pre/*_sj.txt | awk '$5>2' |perl -ne 'chomp;@line=split(/\t/);$line[4]=100;print join("\t",@line),"\n";'| sort | uniq > suva_sj.bed

sh /users/newablife/ablifebio/SUVA/bin/suvasj2junc.sh suva_sj.bed junctions.bed


cd KnownAS && cat $gff |awk '$3~/transcript|mRNA/'| perl -ne 'chomp;@line=split(/\t/);$line[-1]=~/ID=([\w\.\-]+);\S*Parent=([\w\.\-]+)/;$mrna=$1;$gene=$2;print $mrna,"\t",$gene,"\n";'| sort|uniq > isoform_gene.txt && cat isoform_gene.txt > genemode.txt && sed '1 i#IsoformGene' -i genemode.txt && perl /users/chengc/work/ablifedev/src/AS/bin/Anno/gff2bed.pl -i $gff -o gff2bed.txt && perl /users/chengc/work/ablifedev/src/AS/bin/Anno/genemodel2sj.pl -gff $gff -mod genemode.txt -o model2bed.txt && perl /users/chengc/work/ablifedev/src/AS/bin/AS/nas/detectAS.pl -gff $gff -mod genemode.txt -sj gff2bed.txt -o _known_AS_NIR_all.txt && less -S _known_AS_NIR_all.txt | perl -ne 'chomp;@line=split(/\t/);$tmp="$line[1];$line[2];$line[4];$line[5]";@pos=split(/;/,$tmp);@pos2=sort {$a<=>$b} @pos;print $_,"\t",join(";",@pos2),"\n";' > _known_AS_NIR_all_pos.txt && less -S _known_AS_NIR_all_pos.txt | perl -ne 'chomp;@line=split(/\t/);next if($nir{$line[-1]}==1);$nir{$line[-1]}=1;pop @line;print join("\t",@line),"\n";' > _known_AS_NIR.txt && uniq-c 7 _known_AS_NIR.txt > known_AS_NIR.txt.result &&perl /users/chengc/work/ablifedev/src/AS/bin/AS/nas/detectIR_gff.pl -gff $gff -sj model2bed.txt -allsj model2bed.txt -g 1 -o _known_AS_IR_all.txt && less -S _known_AS_IR_all.txt | perl -ne 'chomp;@line=split(/\t/);$tmp="$line[0]:$line[1]:$line[2]:$line[3]";next if($ir{$tmp}==1);$ir{$tmp}=1;print $_,"\n";' > _known_AS_IR.txt &&  uniq-c 5 _known_AS_IR.txt > known_AS_IR.txt.result 
perl /users/ablife/ablife-perl/public/exp_add_anno.pl -exp genemode.txt -anno /data12/Genome/new/Animals/human/annotation/gencode_v33/human_anno.txt -column 1 -annonum auto -o trans_anno.txt
perl /users/ablife/ablife-perl/public/exp_add_anno.pl -exp _known_AS_NIR.txt -anno trans_anno.txt -column -1 -o known_AS_NIR.txt && perl /users/ablife/ablife-perl/public/exp_add_anno.pl -exp _known_AS_IR.txt -anno trans_anno.txt -column -1 -o known_AS_IR.txt

cd ..

###### add intron to gff #####
perl /users/chengc/work/ablifedev/src/AS/bin//Anno/addintro.pl -g $gff -o ./Anno/gff.addintron
#
###### make Annofile #####
perl /users/chengc/work/ablifedev/src/AS/bin//Anno/make_gff_at.pl -i ./Anno/gff.addintron -o ./Anno/annofile
#
###### genemodel to bed #####
perl /users/chengc/work/ablifedev/src/AS/bin//Anno/genemodel2sj.pl -gff ./Anno/gff.addintron -mod ./KnownAS/genemode.txt -o ./Anno/annoSj.bed
#
###### gff to bed #####
perl /users/chengc/work/ablifedev/src/AS/bin//Anno/gff2bed.pl -i ./Anno/gff.addintron -o ./Anno/gffSj.bed

###### anno #####
cat ./Anno/gff.addintron | perl -ne 'chomp;@line=split(/\t/);next if $line[2] ne "intron";$line[3]-=1;$line[4]+=1;$sj="$line[0]:$line[3]:$line[4]:$line[6]";$line[-1]=~/gene_id=([\w|\.]+);transcript_id=([\w|\.]+);/;print "$sj\t$2\t$1\n";' | sort | uniq > ./Anno/sj.isoform && cat ./Anno/sj.isoform | perl -ne 'chomp;@line=split(/\t/);push @{$h{$line[0]}},$line[1];END{foreach $sj(keys %h){print $sj,"\t",join(",",@{$h{$sj}}),"\n";}}' > ./Anno/sj.isoform.final &

cat ./Anno/gff.addintron | perl -ne 'BEGIN{$en=0;}chomp;@line=split(/\t/);if($line[2] eq "transcript"){$en=0;}next if $line[2] ne "exon";$en++;$a="$line[0]:$line[3]:$line[6]";$b="$line[0]:$line[4]:$line[6]";$exon="$line[0]:$line[3]:$line[4]:$line[6]";$line[-1]=~/gene_id=([\w|\.]+);transcript_id=([\w|\.]+);/;print "$a\t$exon\t$en\t$2\t$1\n";print "$b\t$exon\t$en\t$2\t$1\n";' > ./Anno/exon.isoform  && cat ./Anno/exon.isoform | cut -f 1-2 | sort | uniq | perl -ne 'chomp;@line=split(/\t/);push @{$h{$line[0]}},$line[1];END{foreach $sj(keys %h){print $sj,"\t",join(",",@{$h{$sj}}),"\n";}}' > ./Anno/exon.isoform.final &
wait

#
###### get new sj #####
perl /users/chengc/work/ablifedev/src/AS/bin//SJ/get_known_novel_sj.pl -asj ./Anno/gffSj.bed -allsj ./junctions.bed -o ./SJ/new_sj -o2 ./SJ/known_sj ; cat ./junctions.bed | sed '1d' | wc -l  > ./Anno/allsjNum ; cat ./SJ/new_sj | wc -l  > ./Anno/novelsjNum ; cat ./SJ/known_sj | wc -l  > ./Anno/knownsjNum
#
###### known_as_quantify #####
perl /users/chengc/work/ablifedev/src/AS/bin//AS/nas/as_quantify.pl -as ./KnownAS/_known_AS_NIR.txt -sj ./junctions.bed -o ./AS/known_AS_NIR.txt.quantify ;perl /users/chengc/work/ablifedev/src/AS/bin//AS/nas/as_filter.pl -as ./AS/known_AS_NIR.txt.quantify -num_alt_min 2 -num_model_min 2 -num_sum_alt_model 10 -o ./AS/known_AS_NIR.txt.quantify.filter
#
###### get_novel_AS_events #####
perl /users/chengc/work/ablifedev/src/AS/bin//AS/nas/detectAS.pl -gff ./Anno/gff.addintron -mod ./KnownAS/genemode.txt -sj ./SJ/new_sj -o ./AS/novel_AS_NIR_all.txt ;less -S ./AS/novel_AS_NIR_all.txt | perl -ne 'chomp;@line=split(/\t/);$tmp="$line[1];$line[2];$line[4];$line[5]";@pos=split(/;/,$tmp);@pos2=sort {$a<=>$b} @pos;print $_,"\t",join(";",@pos2),"\n";' > ./ASnovel_AS_NIR_all_pos.txt ;less -S ./ASnovel_AS_NIR_all_pos.txt | perl -ne 'chomp;@line=split(/\t/);next if($nir{$line[-1]}==1);$nir{$line[-1]}=1;pop @line;print join("\t",@line),"\n";' > ./AS/novel_AS_NIR.txt;perl /users/chengc/work/ablifedev/src/AS/bin//AS/nas/as_quantify.pl -as ./AS/novel_AS_NIR.txt -sj ./junctions.bed -o ./AS/novel_AS_NIR.txt.quantify ;perl /users/chengc/work/ablifedev/src/AS/bin//AS/nas/as_filter.pl -as ./AS/novel_AS_NIR.txt.quantify -num_alt_min 2 -num_alt_ratio 0.15 -o ./AS/novel_AS_NIR.txt.quantify.filter
#
###### get_available_AS #####
cat ./AS/known_AS_NIR.txt.quantify.filter ./AS/novel_AS_NIR.txt.quantify.filter > ./AS/available_AS_NIR_events ;cat ./AS/available_AS_NIR_events |awk '{print $7}' |sort|uniq -c > ./AS/available_AS_NIR_events_result ;cat ./AS/known_AS_NIR.txt.quantify.filter |awk '{print $7}' |sort|uniq -c > ./AS/available_known_AS_NIR_events_result ;cat ./AS/novel_AS_NIR.txt.quantify.filter |awk '{print $7}' |sort|uniq -c > ./AS/available_novel_AS_NIR_events_result 
#
###### get_NIR_sj #####
perl /users/chengc/work/ablifedev/src/AS/bin//AS/nas/as_get_sj.pl -as ./AS/available_AS_NIR_events -sj ./junctions.bed -o ./AS/available_AS_NIR_events_sj
#
# ###### do EXP and detectIntroR #####
# perl5.16 /users/chengc/work/ablifedev/src/AS/bin//EXP/detect_exon_intron_exp_base.pl -bam /data11/RBP_data/RBP_Analysis_0406/SRSF6/result/mapping/SRSF6_1st_mapping/accepted_hits.uniq.bam -fa /data0/Genome/human/gencode/human.fa -chrlen /data0/Genome/human/gencode/GRCh38.p2_chrLen.txt -anno ./Anno/annofile -od ./Expression & 
# perl5.16 /users/chengc/work/ablifedev/src/AS/bin//AS/nas/detectIR_strand_animals.pl -gff ./Anno/gff.addintron -mod ./KnownAS/genemode.txt -fa /data0/Genome/human/gencode/human.fa -as ./AS/available_AS_NIR_events_sj -bam /data11/RBP_data/RBP_Analysis_0406/SRSF6/result/mapping/SRSF6_1st_mapping/accepted_hits.uniq.bam -chrlen /data0/Genome/human/gencode/GRCh38.p2_chrLen.txt -o ./AS/available_AS_IR_events_all && less -S ./AS/available_AS_IR_events_all | perl -ne 'chomp;@line=split(/\t/);$tmp="$line[0]:$line[1]:$line[2]:$line[3]";next if($ir{$tmp}==1);$ir{$tmp}=1;print $_,"\n";' > ./AS/available_AS_IR_events & 
#  wait 

# #
# ###### get_known_novel_IntroR #####
# perl /users/chengc/work/ablifedev/src/AS/bin//AS/analyze/get_known_novel_as.pl -i ./AS/available_AS_IR_events -gas ./KnownAS/_known_AS_IR.txt -ok ./AS/available_known_AS_IR_events -on ./AS/available_novel_AS_IR_events
# #
# ###### stat_IntroR #####
# cat ./AS/available_AS_IR_events |awk '{print $5}' |sort|uniq -c > ./AS/available_AS_IR_events_result ;cat ./AS/available_known_AS_IR_events |awk '{print $5}' |sort|uniq -c > ./AS/available_known_AS_IR_events_result ;cat ./AS/available_novel_AS_IR_events |awk '{print $5}' |sort|uniq -c > ./AS/available_novel_AS_IR_events_result 
# #
# ###### get_IR_sj #####
cp ./KnownAS/_known_AS_IR.txt ./AS/known_AS_IR.txt &&  perl /users/chengc/work/ablifedev/src/AS/bin//AS/nas/as_get_sj_IR.pl -as ./AS/known_AS_IR.txt -sj ./junctions.bed -o ./AS/available_AS_IR_events_sj
#
###### get_AS_sj #####
cd ./AS && cat available_AS_NIR_events_sj available_AS_IR_events_sj > AS_events_sj

cat AS_events_sj  | grep -v ">" | perl -ne '@line=split;$n=join(":",@line[0..3]);print $n,"\t$line[5]\t$line[6]\n";'  > all_AS_event_sj_anno
#
###### stat all #####
# cd ./AS&&cat available_AS_NIR_events_result available_AS_IR_events_result > AS_events_result&&cat available_known_AS_NIR_events_result available_known_AS_IR_events_result > AS_events_known_result&&cat available_novel_AS_NIR_events_result available_novel_AS_IR_events_result > AS_events_novel_result
#

