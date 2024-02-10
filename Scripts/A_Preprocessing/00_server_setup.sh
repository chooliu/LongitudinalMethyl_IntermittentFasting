# ==============================================================================
# A/00_server_setup.sh
# setting up software, environment
# ==============================================================================



cat > .Renviron
R_LIBS_USER=$HOME/Software/R_libs/
  
  

qlogin -R rusage[mem=48]
module load R/4.1.0_beta
module unload perl/5.26 
module load perl/5.16.3
module load gcc/7.4.0
cd Analyses/DRIFT2_Methylation


R
cd $HOME
mkdir Analyses
mkdir Analyses/DRIFT2_Methylation
cd Analyses/DRIFT2_Methylation


packrat::init()
# packrat::snapshot()
# packrat::on()

install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("bacon")#, lib = "$HOME/R_libs", force = T)
install.packages("data.table")
install.packages("readxl")
install.packages("khroma")


install.packages("Cairo")
install.packages("writexl")


BiocManager::install("DelayedMatrixStats")

install.packages("https://cran.r-project.org/src/contrib/Archive/sp/sp_1.2-0.tar.gz")
install.packages("https://cran.r-project.org/src/contrib/Archive/akima/akima_0.6-1.tar.gz")
install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.2.tar.gz")
install.packages("https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-2.tar.gz")

module load perl/5.26

BiocManager::install("minfi")
BiocManager::install("sesame")
BiocManager::install("sesameData")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("ExperimentHub")

sesameDataCacheAll()

wget https://sourceforge.net/projects/libpng/files/libpng16/1.6.37/libpng-1.6.37.tar.gz --no-check-certificate
tar libpng-1.6.37.tar.gz
cd libpng-1.6.37
./configure --prefix=/beevol/home/lincuini/Software/libpng-1.6.37/
make
make install

  
  
export PATH=/beevol/home/lincuini/Software/:$PATH
export LD_LIBRARY_PATH=/beevol/home/lincuini/Software/libpng-1.6.37/lib/:$LD_LIBRARY_PATH



# Packrat Autoloader (version 0.7.0)
source("packrat/init.R")

# libpng 1.5 installed system-wide; need link to 1.6 for many Bioconductor pkgs
Sys.setenv(LD_LIBRARY_PATH =
             paste0("/beevol/home/lincuini/Software/libpng-1.6.37/lib/:",
                    Sys.getenv("LD_LIBRARY_PATH")))
dyn.load("/beevol/home/lincuini/Software/libpng-1.6.37/lib/libpng16.so.16")


module load gcc/7.4.0 
env  CXXFLAGS="-std=c++11 -O3" CC=/cluster/software/modules-sw/gcc/gcc-7.4.0/bin/gcc CXX=/cluster/software/modules-sw/gcc/gcc-7.4.0/bin/g++ ./bootstrap
env  CXXFLAGS=-std=c++11 CXXFLAGS=-O3 CC=/cluster/software/modules-sw/gcc/gcc-7.4.0/bin/gcc CXX=/cluster/software/modules-sw/gcc/gcc-7.4.0/bin/g++ ./configure --prefix=/beevol/home/lincuini/Software/cmake-3.23.1/ 



BiocManager::install("RUVSeq")
BiocManager::install("ExperimentHub")
BiocManager::install("annotatr")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
 





# wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1.tar.gz
# tar xvf nloptv2.7.1.tar.gz
# 
# export PATH=/beevol/home/lincuini/Software/cmake-3.23.1-linux-x86_64/bin:$PATH
# 
# wget https://github.com/stevengj/nlopt/archive/v2.7.1.tar.gz
# mv v2.7.1.tar.gz nloptv2.7.1.tar.gz
# tar xvf nloptv2.7.1.tar.gz
install.packages("lme4")
install.packages("lmerTest")

mkdir CL_Methylation




cd $HOME/Software
git clone https://github.com/brentp/combined-pvalues.git

module load anaconda/4.7.10
conda create -n py38 python=3.8.5 anaconda
# conda init bash
conda activate py38

# conda install -yc bioconda combined-pvalues


cd $HOME/Software/combined-pvalues
python ez_setup.py
python setup.py install # --prefix=$HOME/Python_libs

python setup.py install --prefix $HOME/Python_Libs/combp/
export PYTHONPATH=$HOME/Python_Libs/combp/lib/python3.8/site-packages/:$PYTHONPATH
python -m pip install setuptools
python 
python -m pip install --upgrade pip --user
conda install toolshed
conda upgrade numpy
cd combined-pvalues

python setup.py --install-dir=$HOME/Python_libs
conda config --append channels conda-forge
conda install toolshed --use-local

conda install libgfortran



# Genotyping

screen
qlogin -R rusage[mem=12]

cd $HOME/Software

wget https://storage.googleapis.com/broad-alkesgroup-public/EIGENSOFT/EIG-6.1.4.tar.gz
tar xvzf EIG-6.1.4.tar.gz

datadir=$HOME/Data/DRIFT2/
analysisdir=$HOME/Analyses/DRIFT2_Methylation/


mkdir $analysisdir/Genotype
cd $analysisdir/Genotype

mkdir raw_data
cd raw_data

cp $datadir/mega/Borengasser_MEGA2_09232019_Final_Deliverable/cleaned_files/* $PWD
