# ==============================================================================
# A/03_run_geno_smartPCA.sh
# genotype assays --> geno PCs
# ==============================================================================




qlogin -R rusage[mem=12]
module load plink # PLINK v1.90b5 64-bit (14 Nov 2017)

mkdir $analysisdir/genotype/plink
cd $analysisdir/genotype/plink



# step 1: remove control samples -----------------------------------------------

awk 'index($2, "NA1")' ../raw_data/Borengasser_GS_092319_updated_callrate_passing_QC.fam |
  awk '{print $1, $2}' > step01_control_samples.txt

plink --bfile ../raw_data/Borengasser_GS_092319_updated_callrate_passing_QC \
--remove step01_control_samples.txt --chr 1-22 --geno 0.03 --mind 0.1 --maf 0.05 --out step01 --make-bed





# step 2: filter out regions of extended LD ------------------------------------
# Price, et al.,
# Long-Range LD Can Confound Genome Scans in Admixed Populations.
# Am J Hum Genet 2008;83:132–135.

cp step01.fam step02.fam
cp step01.bed step02.bed
cp step01.bim step02.bim

cat > step02_longrangeld.txt
1	48000000	52000000  
2	86000000	100500000 
2	134500000	138000000 
2	183000000	190000000 
3	47500000	50000000  
3	83500000	87000000  
3	89000000	97500000  
5	44500000	50500000  
5	98000000	100500000 
5	129000000	132000000 
5	135500000	138500000 
6	25000000	35000000  
6	57000000	64000000   
6	140000000	142500000 
7	55000000	66000000  
8	8000000	12000000  
8	43000000	50000000 
8	112000000	115000000
10	37000000	43000000
11	46000000	57000000
11	87500000	90500000
12	33000000	40000000
12	109500000	112000000
20	32000000	34500000


while IFS=$'\t', read chr start end
do
  printf "awk \'{ if ( ( \$1 == %s) && ( \$4 >= %s && \$4 <= %s) )
    { print \$1, \$4, \$4 } } \' step01.bim \n " \
  $chr $start $end
done < step02_longrangeld.txt  &> step02_filtering.sh

sh step02_filtering.sh > step02_filtered.bed

plink --bfile step01 --exclude step02_filtered.bed --out step02 --make-bed



# re-format for plink format compatibility  ------------------------------------
# no cM information but didn't impact resulting PCs

module load R/3.6.1
R
library(tidyverse)


bim_file <-
  read_tsv("step02.bim", col_names = F) %>%
  rowwise() %>%
  mutate(X2 = c(X1, X4, X5, X6) %>% paste(collapse = "_")) %>%
  group_by(X1) %>% # chr 
  group_split() %>%
  map(function(x) {
    x$X3 <- cummax(x$X3)
    x
  }) %>%
  bind_rows()

write_tsv(bim_file, path = "step02_short.bim", col_names = F)






# step 3: run eigenstrat -------------------------------------------------------

export PATH="/beevol/home/borengas/software/EIG-6.1.4/bin:$PATH"


cat  >  step03_convert
genotypename: step02.bed
snpname: step02_short.bim
indivname: step02_short.fam
outputformat: EIGENSTRAT
genooutfilename: step02_eigenst.geno
snpoutfilename: step02_eigenst.snp
indoutfilename: step02_eigenst.indiv
pordercheck: NO

mkdir ../eigenstrat

cat  >  step03_smartpca_eigenfmt
genotypename: step02_eigenst.geno
snpname: step02_eigenst.snp
indivname: step02_eigenst.indiv
evecoutname: ../eigenstrat/20210512_ADA_eigenvects.out
evaloutname: ../eigenstrat/20210512_ADA_eigenvals.out
snpweightoutname: ../eigenstrat/20210512_ADA_loadings.out
numoutlieriter: 0
nsnpldregress: 2


convertf -p step03_convert > step03_convert.log
smartpca -p step03_smartpca_eigenfmt  > step03_smartpca.log


