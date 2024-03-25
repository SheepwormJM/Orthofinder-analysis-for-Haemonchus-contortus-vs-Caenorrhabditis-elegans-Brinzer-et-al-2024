This contains the code related to the Orthofinder and KAAS/KEGG analysis as published in Brinzer, R. et al 2024 The parasitic nematode Haemonchus contortus lacks molybdenum cofactor synthesis, leading to sulphite sensitivity and lethality in vitro (submitted). 

In brief, proteomes were downloaded for multiple worm species, with sheep as an outgroup. After re-naming fasta headers for consistency, primary transcripts were obtained using OrthoFinder. Next, OrthoFinder was used to produce a rooted species tree (with sheep as outgroup), and to identify putative orthologues within Hierarchical Orthogroups (HOGs). 

Using a mixture of unix/linux commands, and R, results were obtained for the following categories using the N0.tsv file in the HOG directory: 
1. Proteins shared between _Haemonchus contortus_ and _Caenorrhabditis elegans_
2. Proteins present in N0.tsv, but which were not shared between _H. contortus_ and _C. elegans_
3. Proteins that were absent from the N0.tsv file, but which had been present in the primary transcripts

Amino acid sequences for these proteins were then obtained and presented to KAAS/KEGG (https://www.genome.jp/tools/kaas/) for identifcation of enzymes. For _H. contortus_ a bi-directional blast hit was used, while a single-directional blast hit was deemed sufficient for _C. elegans_, as this species is present in the database of KAAS.

The idea for this analysis came from attending a talk given by Dr Pete Steketee, who had compared _Trypanosoma congolense_ vs _Trypanosoma brucei_. See Steketee et al., 2021 https://doi.org/10.1371/journal.ppat.1009734 .

Methods:

1. Install OrthoFinder (https://github.com/davidemms/OrthoFinder#installing-orthofinder-on-linux).

Downloaded ```OrthoFinder.tar.gz v2.5.4``` and transferred to ```~```

Then ```tar xzf OrthoFinder.tar.gz```

Test you can run OrthoFinder: ```./OrthoFinder/orthofinder -h```. OrthoFinder should print its 'help' text.

```
mkdir Orthofinder
cd Orthofinder
mkdir proteomes
cd proteomes
```

2. Get the data:
```
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/caenorhabditis_briggsae/PRJNA10731/caenorhabditis_briggsae.PRJNA10731.WBPS17.protein.fa.gz 
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS17.protein.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/caenorhabditis_remanei/PRJNA53967/caenorhabditis_remanei.PRJNA53967.WBPS17.protein.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS17.protein.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/haemonchus_placei/PRJEB509/haemonchus_placei.PRJEB509.WBPS17.protein.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/necator_americanus/PRJNA72135/necator_americanus.PRJNA72135.WBPS17.protein.fa.gz
```
_No T. axei genome on WBPS. Decided not to use the T. circumcincta genome there (the WASHU) as I know it to be fairly incomplete and fragmented._

3. Make sure the protein fasta files are not gzipped.
```
for f in *fa.gz ; do gunzip $f ; done
```

4. Download the orthofinder file to get the tutorial/tools folder to get the primary_transcript.py script.
```
cd ../
wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz
tar xzvf OrthoFinder.tar.gz
```
Fine.

5. Make the fasta headers suitable for OrthoFinder:

Improve the headers of the protein.fa files to then obtain only primary transcripts - an issue with ```OrthoFinder/tools/primary_transcript.py``` means that some transcripts are identified as multiple genes, and others can be lost completely if header formats vary between different genes.

Note, from a previous run, it was known that certain species files were in the correct header formation for OrthoFinder. If you run primary_transcript.py and get back 'unidentified transcripts' then you want to check and probably correct your headers.

```
sed 's/locus=.*//g' caenorhabditis_elegans.PRJNA13758.WBPS17.protein.fa | sed 's/status=.*//g' > red_caenorhabditis_elegans.PRJNA13758.WBPS17.protein.fa

grep -E '>' red_caenorhabditis_elegans.PRJNA13758.WBPS17.protein.fa > red_transcripts

head red_transcripts
```
```
>2L52.1a wormpep=CE32090 gene=WBGene00007063
>2L52.1b wormpep=CE50569 gene=WBGene00007063
>2RSSE.1a wormpep=CE32785 gene=WBGene00007064
>2RSSE.1b wormpep=CE48524 gene=WBGene00007064
>2RSSE.1c wormpep=CE52653 gene=WBGene00007064
>3R5.1a wormpep=CE24758 gene=WBGene00007065
>3R5.1b wormpep=CE47782 gene=WBGene00007065
>4R79.1a wormpep=CE35820 gene=WBGene00003525
>4R79.1b wormpep=CE39659 gene=WBGene00003525
>4R79.2a wormpep=CE19650 gene=WBGene00007067
```
```
grep -e '>' caenorhabditis_remanei.PRJNA53967.WBPS17.protein.fa > c.remanei.transcripts
grep -E '>' caenorhabditis_briggsae.PRJNA10731.WBPS17.protein.fa > c.briggsae.transcripts
grep -E '>' haemonchus_placei.PRJEB509.WBPS17.protein.fa > h.placei.transcripts
grep -E '>' necator_americanus.PRJNA72135.WBPS17.protein.fa > n.americanus.transcripts
```
```
head c.remanei.transcripts
```
```
>CRE00001 wormpep=RP28513 gene=WBGene00051212 status=Predicted
>CRE00002 wormpep=RP04294 gene=WBGene00051213 locus=Cre-dim-1 status=Partially_confirmed
>CRE00003 wormpep=RP46599 gene=WBGene00051220 status=Predicted
>CRE00004 wormpep=RP30149 gene=WBGene00072597 locus=Cre-sdpn-1 status=Predicted
>CRE00005 wormpep=RP30948 gene=WBGene00077823 status=Predicted
>CRE00006 wormpep=RP32611 gene=WBGene00077870 status=Predicted
>CRE00007 wormpep=RP39037 gene=WBGene00051224 locus=Cre-his-71 status=Predicted
>CRE00008 wormpep=RP45595 gene=WBGene00077894 status=Predicted
>CRE00009 wormpep=RP22848 gene=WBGene00051225 locus=Cre-sem-5 status=Predicted
>CRE00010 wormpep=RP14601 gene=WBGene00051226 locus=Cre-sfxn-2 status=Predicted
```
```
head c.briggsae.transcripts
```
```
>CBG00001 wormpep=CBP43891 gene=WBGene00023521 status=Confirmed uniprot=A8WM42 insdc=CAP21546.2
>CBG00002 wormpep=CBP47711 gene=WBGene00023522 status=Predicted uniprot=A8WM43 insdc=CAP21547.1
>CBG00003 wormpep=CBP49925 gene=WBGene00023523 status=Partially_confirmed uniprot=A8WM44 insdc=CAP21548.2
>CBG00005a wormpep=CBP39399 gene=WBGene00023525 locus=Cbr-fbn-1 status=Confirmed uniprot=A8WM46
>CBG00005b wormpep=CBP41235 gene=WBGene00023525 locus=Cbr-fbn-1 status=Confirmed
>CBG00005c wormpep=CBP43321 gene=WBGene00023525 locus=Cbr-fbn-1 status=Confirmed
>CBG00005d wormpep=CBP45248 gene=WBGene00023525 locus=Cbr-fbn-1 status=Confirmed
>CBG00005e wormpep=CBP38130 gene=WBGene00023525 locus=Cbr-fbn-1 status=Partially_confirmed
>CBG00005f wormpep=CBP44725 gene=WBGene00023525 locus=Cbr-fbn-1 status=Partially_confirmed
>CBG00005g wormpep=CBP46706 gene=WBGene00023525 locus=Cbr-fbn-1 status=Partially_confirmed
```
```
head h.placei.transcripts
```
```
>HPLM_0000000001-mRNA-1 transcript=HPLM_0000000001-mRNA-1 gene=HPLM_0000000001
>HPLM_0000000101-mRNA-1 transcript=HPLM_0000000101-mRNA-1 gene=HPLM_0000000101
>HPLM_0000000201-mRNA-1 transcript=HPLM_0000000201-mRNA-1 gene=HPLM_0000000201
>HPLM_0000000301-mRNA-1 transcript=HPLM_0000000301-mRNA-1 gene=HPLM_0000000301
>HPLM_0000000401-mRNA-1 transcript=HPLM_0000000401-mRNA-1 gene=HPLM_0000000401
>HPLM_0000000501-mRNA-1 transcript=HPLM_0000000501-mRNA-1 gene=HPLM_0000000501
>HPLM_0000000601-mRNA-1 transcript=HPLM_0000000601-mRNA-1 gene=HPLM_0000000601
>HPLM_0000000701-mRNA-1 transcript=HPLM_0000000701-mRNA-1 gene=HPLM_0000000701
>HPLM_0000000801-mRNA-1 transcript=HPLM_0000000801-mRNA-1 gene=HPLM_0000000801
>HPLM_0000000901-mRNA-1 transcript=HPLM_0000000901-mRNA-1 gene=HPLM_0000000901
```
```
head n.americanus.transcripts
```
```
>NAME_00001 transcript=NAME_00001 gene=NAME_00001
>NAME_00002 transcript=NAME_00002 gene=NAME_00002
>NAME_00003 transcript=NAME_00003 gene=NAME_00003
>NAME_00004 transcript=NAME_00004 gene=NAME_00004
>NAME_00005 transcript=NAME_00005 gene=NAME_00005
>NAME_00006 transcript=NAME_00006 gene=NAME_00006
>NAME_00007 transcript=NAME_00007 gene=NAME_00007
>NAME_00008 transcript=NAME_00008 gene=NAME_00008
>NAME_00009 transcript=NAME_00009 gene=NAME_00009
>NAME_00010 transcript=NAME_00010 gene=NAME_00010
```
_Note that N. americanus looks like it doesn't have more that one transcript for any gene._
```
grep -E 'mRNA-2' h.placei.transcripts > hp2 # H. placei also has only a single transcript per gene.
```
Correct the other caenorrhabditis files:
```
sed 's/locus=.*//g' caenorhabditis_remanei.PRJNA53967.WBPS17.protein.fa | sed 's/status=.*//g' > red_caenorhabditis_remanei.PRJNA53967.WBPS17.protein.fa
sed 's/locus=.*//g' caenorhabditis_briggsae.PRJNA10731.WBPS17.protein.fa | sed 's/status=.*//g' > red_caenorhabditis_briggsae.PRJNA10731.WBPS17.protein.fa
```

6. For the outroots:

Sheep! - from the 2021 genome assembly (Rambouliet ewe):
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/772/045/GCF_016772045.1_ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_protein.faa.gz
```
https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giab096/6521876 for the paper on this genome assembly.

Brugia malayi (genome looks very complete and contiguous on WBPS17): 


![image](https://github.com/SheepwormJM/RNAi-LFIT-L4/assets/55552826/12e9d2a6-66d9-44c1-a33a-e1f4420e12d9)

```
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/brugia_malayi/PRJNA10729/brugia_malayi.PRJNA10729.WBPS17.protein.fa.gz
```
```
zcat brugia_malayi.PRJNA10729.WBPS17.protein.fa.gz | grep -E '>' > b.malayi.transcripts # Note this also leaf the file unzipped.

zcat GCF_016772045.1_ARS-UI_Ramb_v2.0_protein.faa.gz | grep -E '>' > sheep.transcripts

zcat GCF_016772045.1_ARS-UI_Ramb_v2.0_protein.faa.gz | head > head_sheep

sed 's/locus=.*//g' brugia_malayi.PRJNA10729.WBPS17.protein.fa | sed 's/status=.*//g' > red_brugia_malayi.PRJNA10729.WBPS17.protein.fa
```
Note - the sheep.transcripts seem to have most with '.1' but there are some '.2' and '.3' I tried pulling out other gene transcripts for one of these '.2' and found only itself. I have emailed the curator to find out whether this file already contains the longest transcripts for each gene, and whether some are simply denoted '.2' or '.3' rather than '.1'. The curator thinks that they are most likely to be the longest transcripts.
```
>NP_001009189.1 UDP-glucuronosyltransferase 1-9 precursor [Ovis aries]
>NP_001009191.1 C-X-C motif chemokine 10 precursor [Ovis aries]
>NP_001009192.1 glycogen phosphorylase, muscle form [Ovis aries]
```
```
grep -E 'NP_001009315' sheep.transcripts > NP_001009315
```
```
more NP_001009315
>NP_001009315.2 somatotropin precursor [Ovis aries]
```

7. Get a single transcript per gene

(first mv the non-reduced files, and the sheep file to another sub-directory, then renamed the red_ files to their original name):
```
for f in *fa ; do python /nfs/users/nfs_j/jm62/lustre118_link/Orthofinder/OrthoFinder/tools/primary_transcript.py $f ; done
```
Perfect! So much better. Note that there are no 'unidentified transcripts' in any file.
```
Looking for "gene=" of "gene:" to identify isoforms of same gene
Found 15203 accessions, 10878 genes, 0 unidentified transcripts
Wrote 10878 genes
/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/brugia_malayi.PRJNA10729.WBPS17.protein.fa

Looking for "gene=" of "gene:" to identify isoforms of same gene
Found 24008 accessions, 20821 genes, 0 unidentified transcripts
Wrote 20821 genes
/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/caenorhabditis_briggsae.PRJNA10731.WBPS17.protein.fa

Looking for "gene=" of "gene:" to identify isoforms of same gene
Found 28447 accessions, 19985 genes, 0 unidentified transcripts
Wrote 19985 genes
/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/caenorhabditis_elegans.PRJNA13758.WBPS17.protein.fa

Looking for "gene=" of "gene:" to identify isoforms of same gene
Found 31457 accessions, 31438 genes, 0 unidentified transcripts
Wrote 31438 genes
/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/caenorhabditis_remanei.PRJNA53967.WBPS17.protein.fa

Looking for "gene=" of "gene:" to identify isoforms of same gene
Found 21320 accessions, 19621 genes, 0 unidentified transcripts
Wrote 19621 genes
/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/haemonchus_contortus.PRJEB506.WBPS17.protein.fa

Looking for "gene=" of "gene:" to identify isoforms of same gene
Found 21928 accessions, 21928 genes, 0 unidentified transcripts
Wrote 21928 genes
/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/haemonchus_placei.PRJEB509.WBPS17.protein.fa

Looking for "gene=" of "gene:" to identify isoforms of same gene
Found 15728 accessions, 15728 genes, 0 unidentified transcripts
Wrote 15728 genes
/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/necator_americanus.PRJNA72135.WBPS17.protein.fa
```


8.  Update the names to, e.g. H_contortus.fa for the primary transcripts, to keep tidy. Include the sheep sequences (need to gunzip):

9. Run Orthofinder

```
#!/bin/bash
#BSUB -q long
#BSUB -o orthofinder.o
#BSUB -e orthofinder.e
#BSUB -J orthofinder
#BSUB -R "rusage[mem=10000] span[hosts=1] select[mem>10000]"
#BSUB -n 2
#BSUB -M 10000
#BSUB -cwd "/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts"


/nfs/users/nfs_j/jm62/lustre118_link/Orthofinder/OrthoFinder/orthofinder -f /lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/
```
It has worked! Took about 4 hours?

10. Let's sanity check the species tree first of all.
 ```
# Transfer the file to computer to be able to visualise:

pscp -P 2227 jm62@localhost:/nfs/users/nfs_j/jm62/lustre118_link/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09/Species_Tree/SpeciesTree_rooted_node_labels.txt ./SpeciesTree_rooted_node_labels.txt
```

![image](https://github.com/SheepwormJM/RNAi-LFIT-L4/assets/55552826/3d95181a-8fa3-4416-a6d2-32762ca200e1)

Good.

11. Get the C. elegans and H. con gene lists

```
awk -F '\t' 'BEGIN {OFS = FS} {print $1, $2, $3, $6, $9}' N0.tsv > CE_HC.N0.tsv #THIS WORKS! The begin is required.

awk -F '\t' "{print NF}" < CE_HC.N0.tsv | uniq > NF_in_CE_HC.N0.tsv # Gave me a total of 5 field per line, for all lines. PERFECT! 
```

Let's try in Rv4.1.0 on the farm5:

```
df<- read.table("CE_HC.N0.tsv", header=T, sep="\t", na.strings ="NA") # Include the sep and the na.strings to load this table with both tabs and commas and spaces, and with blank fields, correctly.
summary(df) # Note that it has just left blank fields blank (from looking at a test)  
```
```
> summary(df)
     HOG                 OG            Gene.Tree.Parent.Clade
 Length:25698       Length:25698       Length:25698
 Class :character   Class :character   Class :character
 Mode  :character   Mode  :character   Mode  :character
  C_elegans         H_contortus
 Length:25698       Length:25698
 Class :character   Class :character
 Mode  :character   Mode  :character
```
```
new_df_Ce<- subset(df, df$H_contortus == "") # Get all lines without any H.con genes
head(new_df_Ce)
summary(new_df_Ce) 

Ce_only<- subset(new_df_Ce, new_df_Ce$C_elegans != "") 
head(Ce_only) # Worked perfectly
summary(Ce_only) 
```
```
> summary(Ce_only)
     HOG                 OG            Gene.Tree.Parent.Clade
 Length:4667        Length:4667        Length:4667
 Class :character   Class :character   Class :character
 Mode  :character   Mode  :character   Mode  :character
  C_elegans         H_contortus
 Length:4667        Length:4667
 Class :character   Class :character
 Mode  :character   Mode  :character
```
```
Ce_genes<-Ce_only[,4] # Get column four
head(Ce_genes) # Perfect, note that it has removed the header though
write.csv(Ce_genes,"Ce_genes_not_Hc.csv", row.names = FALSE, col.names=FALSE, quote=FALSE) 



# To get genes found in H. contortus but not in C. elegans: 
new_df_Hc<- subset(df, df$C_elegans == "") # Get all lines without any c. elegans genes
head(new_df_Hc)
summary(new_df_Hc) 


# To get just the Hc only genes: 

Hc_only<-subset(new_df_Hc, new_df_Hc$H_contortus != "") # Get those rows which have Hcon genes
head(Hc_only) # worked fine. Two rows without either species have been removed.
summary(Hc_only)

```
```
> summary(Hc_only)
     HOG                 OG            Gene.Tree.Parent.Clade
 Length:3076        Length:3076        Length:3076
 Class :character   Class :character   Class :character
 Mode  :character   Mode  :character   Mode  :character
  C_elegans         H_contortus
 Length:3076        Length:3076
 Class :character   Class :character
 Mode  :character   Mode  :character
```
```
Hc_genes<-Hc_only[,5] # Get column five
head(Hc_genes) # Perfect, note that it has removed the header though
write.csv(Hc_genes,"Hc_genes_not_Ce.csv", row.names = FALSE, col.names=FALSE, quote=FALSE) 
```
```
# To get those genes sharing orthogroups in both C. elegans and H. contortus.
new<-subset(df, df$C_elegans != "") # Get those rows which have C.elegans genes
Hc_and_Ce<-subset(new, new$H_contortus != "") # Get those rows which also have H. contortus genes. 
head(Hc_and_Ce)
summary(Hc_and_Ce)
```
```
> summary(Hc_and_Ce)
     HOG                 OG            Gene.Tree.Parent.Clade
 Length:8440        Length:8440        Length:8440
 Class :character   Class :character   Class :character
 Mode  :character   Mode  :character   Mode  :character
  C_elegans         H_contortus
 Length:8440        Length:8440
 Class :character   Class :character
 Mode  :character   Mode  :character
```
```
Hc_and_Ce_genes<-Hc_and_Ce[,4] # Get the C. elegans gene names for the shared orthologues/homologues

Hc_and_Ce_genes_HCON<-Hc_and_Ce[,5] # Get the H.con gene names for the shared orthologues/homologues
write.csv(Hc_and_Ce_genes,"CE_genes_of_Hc_and_CE_HOGs.csv", row.names = FALSE, col.names=FALSE, quote=FALSE) 
write.csv(Hc_and_Ce_genes_HCON,"HCON_genes_of_Hc_and_CE_HOGs.csv", row.names = FALSE, col.names=FALSE, quote=FALSE)
```
All still have a header - the col.names=FALSE didn't work (think it does if a txt file)
```
tail -n +2 Ce_genes_not_Hc.csv > Ce_genes_not_Hc.txt # Removes the header
wc -l Ce_genes_not_Hc.txt

tail -n +2 Hc_genes_not_Ce.csv > Hc_genes_not_Ce.txt # Removes the header
wc -l Hc_genes_not_Ce.txt

tail -n +2 CE_genes_of_Hc_and_CE_HOGs.csv > CE_genes_of_Hc_and_CE_HOGs.txt # Removes the header
wc -l CE_genes_of_Hc_and_CE_HOGs.txt

tail -n +2 HCON_genes_of_Hc_and_CE_HOGs.csv > HCON_genes_of_Hc_and_CE_HOGs.txt # Removes the header
wc -l HCON_genes_of_Hc_and_CE_HOGs.txt
```
To make it easy - ensure the first 20 rows have 3 or fewer genes in them.
```
less Ce_genes_not_Hc.txt # Many have only one or a few genes on them.
wc -l Ce_genes_not_Hc.txt # 4667
tail -n 20 Ce_genes_not_Hc.txt > bottom_Ce_genes_not_Hc.txt
head -n 4647 Ce_genes_not_Hc.txt > top_Ce_genes_not_Hc.txt
cat bottom_Ce_genes_not_Hc.txt top_Ce_genes_not_Hc.txt > Ce_genes_reorder.csv # Puts the last 20 rows, which conveniently have just one or two genes per row, first. 


less Hc_genes_not_Ce.txt # Many have only one or a few genes on them.
wc -l Hc_genes_not_Ce.txt # 3076
tail -n 20 Hc_genes_not_Ce.txt > bottom_Hc_genes_not_Ce.txt
head -n 3056 Hc_genes_not_Ce.txt > top_Hc_genes_not_Ce.txt
cat bottom_Hc_genes_not_Ce.txt top_Hc_genes_not_Ce.txt > Hc_genes_reorder.csv # Puts the last 20 rows, which conveniently have just one or two genes per row, first. 


less CE_genes_of_Hc_and_CE_HOGs.txt # Many have only one or a few genes on them.
wc -l CE_genes_of_Hc_and_CE_HOGs.txt # 8440
tail -n 20 CE_genes_of_Hc_and_CE_HOGs.txt > bottom_CE_genes_of_Hc_and_CE_HOGs.txt
head -n 8420 CE_genes_of_Hc_and_CE_HOGs.txt > top_CE_genes_of_Hc_and_CE_HOGs.txt
cat bottom_CE_genes_of_Hc_and_CE_HOGs.txt top_CE_genes_of_Hc_and_CE_HOGs.txt > CE_AND_HC_genes_reorder.csv # Puts the last 20 rows, which conveniently have just one or two genes per row, first. 



less HCON_genes_of_Hc_and_CE_HOGs.txt # Many have only one or a few genes on them.
wc -l HCON_genes_of_Hc_and_CE_HOGs.txt # 8440
tail -n 20 HCON_genes_of_Hc_and_CE_HOGs.txt > bottom_HCON_genes_of_Hc_and_CE_HOGs.txt
head -n 8420 HCON_genes_of_Hc_and_CE_HOGs.txt > top_HCON_genes_of_Hc_and_CE_HOGs.txt
cat bottom_HCON_genes_of_Hc_and_CE_HOGs.txt top_HCON_genes_of_Hc_and_CE_HOGs.txt > HCON_of_CE_AND_HC_genes_reorder.csv # Puts the last 20 rows, which conveniently have just one or two genes per row, first. 
```
In R, let's get the csv into a single list for use elsewhere.
```
Ce<-read.csv("Ce_genes_reorder.csv", header=F) 
head(Ce)
summary(Ce)

Ce1<-Ce[,1]
Ce2<-Ce[,2]

new<-c(Ce1, Ce2)
write.table(new,"Ce_gene_list.txt",quote = FALSE,row.names = FALSE, col.names=FALSE) # No header! :D 


Hc<-read.csv("Hc_genes_reorder.csv", header=F) 
head(Hc) # Mwah ha ha! Just two columns :D 
summary(Hc)

Hc1<-Hc[,1]
Hc2<-Hc[,2]

new<-c(Hc1, Hc2)

write.table(new,"Hc_gene_list.txt",quote = FALSE,row.names = FALSE, col.names=FALSE) # No header! :D 


Ce<-read.csv("CE_AND_HC_genes_reorder.csv", header=F) 
head(Ce)
summary(Ce)

Ce1<-Ce[,1] # There is only a single column! But doing it this way anyway... paranoid

write.table(Ce1,"Ce_gene_list_of_genes_homologous_to_Hc_genes.txt",quote = FALSE,row.names = FALSE, col.names=FALSE) # No header! :D 

HCON<-read.csv("HCON_of_CE_AND_HC_genes_reorder.csv", header=F) 
head(HCON)
summary(HCON)

HCON1<-HCON[,1] # There is only a single column! But doing it this way anyway... paranoid

write.table(HCON1,"HCON_gene_list_of_genes_homologous_to_Ce_genes.txt",quote = FALSE,row.names = FALSE, col.names=FALSE) # No header! :D 
```
Check lengths, remove blank spaces, check lengths, remove duplicates (hopefully none), check length.
```
awk NF Ce_gene_list.txt | wc -l # This should remove the blank lines and the fields are empty and so not printed. Let's see the number remaining. 

# 8296

sort Ce_gene_list.txt  | uniq | wc -l # Includes a blank field, hence the extra 1. 
# 8297

awk NF Ce_gene_list.txt | sort | uniq > FINAL_Ce_not_Hc_genes.txt # This contains all the genes, but not the blank spaces. Great!


awk NF Hc_gene_list.txt | wc -l # This should remove the blank lines and the fields are empty and so not printed. Let's see the number remaining. 

# 7836

sort Hc_gene_list.txt  | uniq | wc -l # Includes a blank field, hence the extra 1. 
# 7837

awk NF Hc_gene_list.txt | sort | uniq > FINAL_Hc_not_Ce_genes.txt # This contains all the genes, but not the blank spaces. Great! 


awk NF Ce_gene_list_of_genes_homologous_to_Hc_genes.txt | wc -l # This should remove the blank lines and the fields are empty and so not printed. Let's see the number remaining. 

# 10525

sort Ce_gene_list_of_genes_homologous_to_Hc_genes.txt  | uniq | wc -l # Remember that it was already just a single column, so don't expect any blank spaces.
# 10525

awk NF Ce_gene_list_of_genes_homologous_to_Hc_genes.txt | sort | uniq > FINAL_Ce_genes_homologous_to_Hc_genes.txt # This contains all the genes, but not the blank spaces. Great! 


awk NF HCON_gene_list_of_genes_homologous_to_Ce_genes.txt | wc -l # This should remove the blank lines and the fields are empty and so not printed. Let's see the number remaining. 

# 10118

sort HCON_gene_list_of_genes_homologous_to_Ce_genes.txt  | uniq | wc -l # Remember that it was already just a single column, so don't expect any blank spaces.
# 10118

awk NF HCON_gene_list_of_genes_homologous_to_Ce_genes.txt | sort | uniq > FINAL_HCON_genes_homologous_to_Ce_genes.txt # This contains all the genes, but not the blank spaces. Great! 
```
```
HCON, not in C. elegans:
File = FINAL_Hc_not_Ce_genes.txt 
Number of genes = 7836 
Number of HOGs = 3076

CELE, not in H. contortus:
File = FINAL_Ce_not_Hc_genes.txt
Number of genes =  8296
Number of HOGs = 4667

CELE genes, in same HOGs as Hcon genes:
File =  FINAL_Ce_genes_homologous_to_Hc_genes.txt 
Number of genes = 10525
Number of HOGs = 8440

HCON genes, in same HOGs as Cele genes:
File =  FINAL_HCON_genes_homologous_to_Ce_genes.txt 
Number of genes = 10118
Number of HOGs = 8440
```
12. To get the aa seq for KAAS
    
Get the amino acid sequences for the C. elegans gene lists:
```
module load samtools/1.14--hb421002_0
```
```
#!/bin/bash
#BSUB -q normal
#BSUB -o getgenes_Ce_shared_Hc.o
#BSUB -e getgenes_Ce_shared_Hc.e
#BSUB -J getgenes_Ce_shared_Hc
#BSUB -R "rusage[mem=1000] span[hosts=1] select[mem>1000]"
#BSUB -n 1
#BSUB -M 1000
#BSUB -cwd "/nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces"

samtools -v > version_samtools_extract_aa

while read NAME; do samtools faidx /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/C_elegans.fa ${NAME} >> Ce_shared_Hc_N0nospaces.aa.fa; done < /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/FINAL_Ce_genes_homologous_to_Hc_genes.txt_nospace
```
```
#!/bin/bash
#BSUB -q normal
#BSUB -o getgenes_Ce_not_Hc.o
#BSUB -e getgenes_Ce_not_Hc.e
#BSUB -J getgenes_Ce_not_Hc
#BSUB -R "rusage[mem=1000] span[hosts=1] select[mem>1000]"
#BSUB -n 1
#BSUB -M 1000
#BSUB -cwd "/nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces"

samtools -v > version_samtools_extract_aa

while read NAME; do samtools faidx /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/C_elegans.fa ${NAME} >> Ce_NOT_Hc_N0nospaces.aa.fa; done < /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/FINAL_Ce_not_Hc_genes.txt_nospace
```
Get the list of genes not present in N0 but in the celegans proteome:
```
# Get list of all genes
grep -E '>' C_elegans.fa | sed 's/>//g' > Celegenes.txt 
# Get just those that are not in the N0 files: 
cat FINAL_Ce_not_Hc_genes.txt_nospace FINAL_Ce_genes_homologous_to_Hc_genes.txt_nospace Celegenes.txt | sort | uniq -u | wc -l 
#1164 genes 
cat FINAL_Ce_not_Hc_genes.txt_nospace FINAL_Ce_genes_homologous_to_Hc_genes.txt_nospace Celegenes.txt | sort | uniq -u > Cele_genes_not_in_N0.txt
```
```
#!/bin/bash
#BSUB -q normal
#BSUB -o getgenes_Ce_not_N0.o
#BSUB -e getgenes_Ce_not_N0.e
#BSUB -J getgenes_Ce_not_N0
#BSUB -R "rusage[mem=1000] span[hosts=1] select[mem>1000]"
#BSUB -n 1
#BSUB -M 1000
#BSUB -cwd "/nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces"

samtools -v > version_samtools_extract_aa

while read NAME; do samtools faidx /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/C_elegans.fa ${NAME} >> Ce_NOT_N0_nospaces.aa.fa; done < /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/Cele_genes_not_in_N0.txt
```
Get the amino acid sequences of the H. contortus genes in order to BLAST them:
```
#!/bin/bash
#BSUB -q normal
#BSUB -o getgenes_Hc_not_ce.o
#BSUB -e getgenes_Hc_not_ce.e
#BSUB -J getgenes_Hc_not_ce
#BSUB -R "rusage[mem=1000] span[hosts=1] select[mem>1000]"
#BSUB -n 1
#BSUB -M 1000
#BSUB -cwd "/nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces"

samtools -v > version_samtools_extract_aa

while read NAME; do samtools faidx /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/H_contortus.fa ${NAME} >> Hc_NOT_Ce_nospaces.aa.fa; done < /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/FINAL_Hc_not_Ce_genes.txt_nospace
```
```
#!/bin/bash
#BSUB -q normal
#BSUB -o getgenes_Hcshared.o
#BSUB -e getgenes_Hcshared.e
#BSUB -J getgenes_Hcshared
#BSUB -R "rusage[mem=1000] span[hosts=1] select[mem>1000]"
#BSUB -n 1
#BSUB -M 1000
#BSUB -cwd "/nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces"

samtools -v > version_samtools_extract_aa

while read NAME; do samtools faidx /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/H_contortus.fa ${NAME} >> Hc_SHARED_Ce_nospaces.aa.fa; done < /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/FINAL_HCON_genes_homologous_to_Ce_genes.txt_nospace
```
```
# Get list of all genes
grep -E '>' H_contortus.fa | sed 's/>//g' > Hcongenes.txt 
# Get just those that are not in the N0 files: 
cat FINAL_HCON_genes_homologous_to_Ce_genes.txt_nospace FINAL_Hc_not_Ce_genes.txt_nospace Hcongenes.txt | sort | uniq -u | wc -l 
#1667 genes 
cat FINAL_HCON_genes_homologous_to_Ce_genes.txt_nospace FINAL_Hc_not_Ce_genes.txt_nospace Hcongenes.txt | sort | uniq -u > Hcon_genes_not_in_N0.txt
```
```
#!/bin/bash
#BSUB -q normal
#BSUB -o getgenes_HcnotN0.o
#BSUB -e getgenes_HcnotN0.e
#BSUB -J getgenes_HcnotN0
#BSUB -R "rusage[mem=1000] span[hosts=1] select[mem>1000]"
#BSUB -n 1
#BSUB -M 1000
#BSUB -cwd "/nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces"

samtools -v > version_samtools_extract_aa

while read NAME; do samtools faidx /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/H_contortus.fa ${NAME} >> Hc_NOT_N0_nospaces.aa.fa; done < /nfs/users/nfs_j/jm62/scratch/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/Hcon_genes_not_in_N0.txt
```
```
#!/bin/bash
#BSUB -q normal
#BSUB -o orthofinder.o
#BSUB -e orthofinder.e
#BSUB -J orthofinder
#BSUB -R "rusage[mem=1000] span[hosts=1] select[mem>1000]"
#BSUB -n 1
#BSUB -M 1000
#BSUB -cwd "/lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/Orthofinder/Results_Nov09/Phylogenetic_Hierarchical_Orthogroups"

samtools -v > version_samtools_extract_aa

while read NAME; do samtools faidx /lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/C_elegans.fa ${NAME} >> Ce_NOT_Hc.aa.fa; done < /lustre/scratch118/infgen/team333/jm62/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09/Phylogenetic_Hierarchical_Orthogroups/FINAL_Ce_not_Hc_genes.txt

# The output file is the <file_name>.aa.fa
```
```
pscp -P 2227 jm62@localhost:~/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/Hc_NOT_N0_nospaces.aa.fa ./Hc_NOT_N0_nospaces.aa.fa

pscp -P 2227 jm62@localhost:~/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/Hc_SHARED_Ce_nospaces.aa.fa ./Hc_SHARED_Ce_nospaces.aa.fa

pscp -P 2227 jm62@localhost:~/Orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Nov09_Hconvsothers/Phylogenetic_Hierarchical_Orthogroups/FINAL_genelists_celegans_hcon_nospaces/Hc_NOT_Ce_nospaces.aa.fa ./Hc_NOT_Ce_nospaces.aa.fa
```
13. Run KAAS:

**KAAS/KEGG for HCON, no spaces files (March 2023)**

Note that before I seem to have got from BIOMART, but this time I extracted from the .fa file as detailed above:

BBH into hsa, dme, cel, ath, sce, cho, eco, nme, hpy, rpr, bsu, lla, cac, mge, mtu, ctr, bbu, syn, bth, dra, aae, mja, ape, cbr, bmy, loa, nai, tsp, smm, isc

**KAAS/KEGG for Celegans**

Note that as I was using C. elegans genes and the C. elegans gene data set I just did a single directional BLAST (i.e. gene A into genome B).

```
Query info

ID :	1668422708
query name :	Ce_NOT_Hc.aa.fa_14.11.21
query type :	pep
program :	BLAST
method :	SBH
GENES data set :	cel
```
