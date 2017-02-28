# test001

```sql
SELECT
arrays.array_source_accession, arrays.array_family, arrays.array_score, arrays.array_cas_genes, spacers.array_source_accession,  spacers.array_family, spacers.spacer_seq, spacers.spacer_ID, spacers.spacer_length, spacers.array_index,  spacers.array_ID,  spacers.array_repeat,  spacers.array_source_definition, spacers.array_source_species, spacers.protospacer_IDs, arrays.array_index
FROM
CRISPR_spacers AS spacers
INNER JOIN CRISPR_arrays AS arrays ON spacers.array_ID = arrays.array_ID
WHERE arrays.array_cas_genes LIKE '%Csy1%'
AND arrays.array_score >= 1.5
ORDER BY
spacers.array_source_accession ASC

```

```sql
SELECT
arrays.array_source_accession, arrays.array_family, arrays.array_score, spacers.spacer_seq, spacers.array_source_accession, spacers.spacer_ID, spacers.spacer_length, spacers.array_index, spacers.array_ID, spacers.array_repeat, spacers.array_family, spacers.array_source_definition, spacers.array_source_species, spacers.protospacer_IDs,
arrays.array_cas_genes
FROM
CRISPR_spacers AS spacers
INNER JOIN CRISPR_arrays AS arrays ON spacers.array_ID = arrays.array_ID
WHERE arrays.array_cas_genes LIKE '%Csy1%'
AND arrays.array_cas_genes NOT LIKE '%Cas8a1%'
AND arrays.array_cas_genes NOT LIKE '%Cas8a2%'
AND arrays.array_cas_genes NOT LIKE '%Cas8b%'
AND arrays.array_cas_genes NOT LIKE '%Cas8c%'
AND arrays.array_cas_genes NOT LIKE '%Cas10d%'
AND arrays.array_cas_genes NOT LIKE '%Cse1%'
AND arrays.array_cas_genes NOT LIKE '%Cas9%'
AND arrays.array_cas_genes NOT LIKE '%Cas10%'
AND arrays.array_cas_genes NOT LIKE '%Csm2%'
AND arrays.array_cas_genes NOT LIKE '%Cmr5%'
AND arrays.array_score >= 1.5
ORDER BY
spacers.array_source_accession ASC

```

```bash
Cut –f1 I-Fonly_mySQL > I-Fonly_NCs

sed -e 's/”//g' I-Fonly_NCs > I-FonlyNCs1

Uniq I-FonlyNCs1 > I-Fonly_NCs

```

```bash
Seqret I-FonlyBatchEntrez.gbk > I-FonlyBatchEntrez.fa

```

```bash
Perl CRISPDetect_V.1.0.pl –f I-FonlyBatchEntrez.fa –o I-Fonly_arrays 

```

```bash
python GBK_to_faa.py

```

```bash
makeprofileDB –in Cas_siggenes.pn

```

```bash
rpsblast -db Cas_siggenes.pn \
-query /Volumes/BiochemXsan/scratch/brownlab/bradley/I-Fonly_BatchEntrez.fasta \
-outfmt 6 -out /Volumes/BiochemXsan/scratch/brownlab/bradley/I-Fonly_BatchEntrez.results \
-evalue 0.000001 -max_hsps 1

```

```bash
reformat.sh in=PHAST.fa out=PHAST_Subset1.fa samplerate=0.1

```

```bash
Fasta-dinucleotide-shuffle –f PHAST_Subset1.fa –c 10

```

```R
data<-read.delim(file.choose(),skip=3,as.is=TRUE)
data<-tbl_df(data)  
table<-data%>%select(X.Spacer_ID,Spacer_index,Protospacer_start,Score,Protospacer_Strand,Spacer_description,Protospacer_description)%>%filter(Score>=25)
table<-table%>%filter(Score>=25)%>%mutate(Name=paste(Spacer_description,Protospacer_description))
hits.sig<-table%>%group_by(Name)%>%summarise(Counts=n())%>%arrange(desc(Counts))%>%filter(Counts>=3)
clust.table<-table%>%filter(Name==as.character(hits.sig[1,1]))%>%arrange(Protospacer_start)%>%mutate(Cluster=ifelse((lead(Protospacer_start)-lag(Protospacer_start)<=5000),1,0))
for (i in 2:nrow(hits.sig)){
  dist.table<-table%>%filter(Name==as.character(hits.sig[i,1]))%>%arrange(Protospacer_start)%>%mutate(Cluster=ifelse((lead(Protospacer_start)-lag(Protospacer_start)<=5000),1,0))
  clust.table<-bind_rows(clust.table,dist.table)
}
write.csv(clust.table,"Data.csv")

```

```R
vector<-clust.table%>%select(Name,Protospacer_start)

mean(Vector=n)

Table_redrm<-read.delim(file.choose(),skip=3,as.is=TRUE)
Vector_redrm <- Table_redrm[,2]
names(Vector_redrm) <- Table_redrm[,1]
species <- sapply(strsplit(names(Vector_redrm), " "), function(x) x[[1]])
phage <- sort(unique(Table_redrm[,1]))

Vector_redrm[phage==phage[1]]

null <- runif(length(Vector_redrm))
ks.result <- ks.test(Vector_redrm, null)

# test against null dist(uniform)
ks.result.stat <- rep(NA, length(phage))
ks.result.pval <- rep(NA, length(phage))
names(ks.result) <- phage
null_dist <- runif(length(Vector_redrm))
for(ii in 1:length(phage)){
  ks.result.stat[ii] <- ks.test(Vector_redrm[Table_redrm[,1]==phage[ii]], null_dist)$statistic
  ks.result.pval[ii] <- ks.test(Vector_redrm[Table_redrm[,1]==phage[ii]], null_dist)$p.value
  print(ii)
}
ks.result.pval
ks.result.stat
table(ks.result.pval<0.05)
objct<-table(p.adjust(ks.result.pval, method="fdr")<0.05)

```

```R
p<- ggplot(All, aes(x=Position, y=Spacer_Number)) 
+ geom_point(data = subset(All, Name %in% c("Actinobacillus-pleuropneumoniae L20 serotype 5b complete genome PHAST_3892_2_2172701_2264986_ADOF00000000_PHAGE_Haemop_SuMu_")), aes(colour=Protospacer_Strand, size =3, shape=PAM_status)) 
+ coord_cartesian(xlim = c(-0, 92285)) + ggtitle("A.pleuropneumoniae L20 PHAST_3892_2") + geom_text(data = subset(All, Name %in% c("Actinobacillus-pleuropneumoniae L20 serotype 5b complete genome PHAST_3892_2_2172701_2264986_ADOF00000000_PHAGE_Haemop_SuMu_")), aes(x=Position, y=Spacer_Number, label = Score, hjust =1.5)) 
+ geom_vline(xintercept = 11498.5) + scale_y_discrete("Spacer Number")
print(p)
g1 =  p + facet_grid(Array ~ .,  scales = "free", space="free") 
print(g1)

```

```bash
./mg-download.py –metagenome 4522044.3

```

```bash
velveth /Volumes/BiochemXsan/scratch/brownlab/bradley/BergMiller_filteredfasta/ 25  -short *.fna

velvetg /Volumes/BiochemXsan/scratch/brownlab/bradley/BergMiller_filteredfasta/

```

```bash
velveth /Volumes/BiochemXsan/scratch/brownlab/bradley/Ross_filteredFastQ/ 31 -shortPaired -fastq Bovine*R1.fastq Bovine*R2.fastq

velvetg /Volumes/BiochemXsan/scratch/brownlab/bradley/BergMiller_filteredfasta/

```

```bash
Phython quast.py –o /Volumes/BiochemXsan/scratch/brownlab/Bradley/QUAST/BRV1_assembly/qualityassessment –m 1 BRV1_assembly.fa

```

```bash
Pyfasta split –n 5 BOR.fa

```

```bash
PriceTI -spf Sequences.fa 25 200 -picf 3 Seed.fna 1 1 1 -nc 3 -o MethanoBM_output.fa

```

```bash
blastn -evalue 10 -word_size 7 -db /DB/Test/test.fna -query AssembledContig_BergMillerMethano.fa > AssembledContig_BergMillerMethano.out

```



