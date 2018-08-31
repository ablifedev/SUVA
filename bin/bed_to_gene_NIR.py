import sys,os,re
import HTSeq

for paramIndex in range(0, len(sys.argv)):
    if sys.argv[paramIndex]=="-g":
        gff_file=sys.argv[paramIndex+1]
    elif sys.argv[paramIndex]=="-b":
        bed_file=sys.argv[paramIndex+1]
    elif sys.argv[paramIndex]=="-o":
        out_file=sys.argv[paramIndex+1]

genes = ('gene','lincRNA_gene','miRNA_gene','mt_gene','processed_pseudogene','pseudogene','rRNA_gene','snoRNA_gene','snRNA_gene')
chr_list=[]
hash1={}
chr_lines=os.popen("less -S "+bed_file+" |cut -f 1|sort -u").readlines()
for chr in chr_lines:
    chr=chr.strip()
    chr_list.append(chr)
    hash1[chr]=HTSeq.GenomicArrayOfSets([chr], stranded=True)

#使用HTSeq载入gff文件
f1=open(gff_file,"r")
while 1:
    line=f1.readline()
    if line=="":
        break
    else:
        if line.startswith("#"):
            continue
        else:
            items=line.split("\t")
            if items[2] in genes and items[0] in chr_list:
                start=items[3]
                end=items[4]
                chr=items[0]
                strand=items[6]
                # gene_name=items[8].split(";")[0].split("=")[1]
                gene_name=""
                match = re.search(r'ID=([\w\.\-]+)', items[8])
                if not match:
                    continue
                gene_name = str(match.group(1))
                gene_iv = HTSeq.GenomicInterval(chr, int(start) - 1, int(end), strand)
                #print(chr + " " + gene_name + " " + start + ":" + end)
                for iv,fs in hash1[chr][gene_iv].steps():
                    if fs==set():
                        hash1[chr][iv] = gene_name
                    else:
                        hash1[chr][iv] = fs + "," + gene_name
print("finish read gff file.")


#将bed文件逐行比对到HTSeq
f3=open(out_file,"w")
f2=open(bed_file,"r")
while 1:
    line=f2.readline()
    if line=="":
        break
    else:
        items = line.split("\t")
        start = items[1]
        end = items[2]
        chr = items[0]
        bed_len=abs(int(start) - int(end))
        strand = items[5].strip()
        gene_iv = HTSeq.GenomicInterval(chr, int(start) - 1, int(end), strand)
        tmp=[]
        hash2 = {}
        for  iv,fs in hash1[chr][gene_iv].steps():
            if fs==set():
                gene="intergenic"
                map_len = abs(iv.start - iv.end)
                percent = "%.2f" % (map_len / (bed_len + 1) * 100)
                tmp.append(gene + "(" + percent + "%)")

            else:
                gene=fs
                map_len=abs(iv.start-iv.end)
                percent="%.2f"%(map_len/(bed_len+1)*100)
                tmp.append(gene + "(" + percent + "%)")
            for ge in gene.split(","):
                if ge =="":
                    pass
                else:
                    if ge not in hash2:
                        hash2[ge] = map_len / (bed_len + 1) * 100
                    else:
                        hash2[ge] += map_len / (bed_len + 1) * 100
        flag=0
        for ge in hash2.keys():
            if hash2[ge]==max(hash2.values()):
                tmp_1=ge
                flag=1
            else:
                pass
        if flag==0:
            tmp_1="error"
            print("error")
        f3.write(line.strip() + "\t" + tmp_1 + "\t" + "---".join(tmp) + "\n")
f1.close()
f2.close()
f3.close()


