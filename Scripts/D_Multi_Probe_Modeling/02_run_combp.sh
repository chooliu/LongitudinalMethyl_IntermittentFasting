# ==============================================================================
# D/02_run_combp.R
# apply combp for .bed files from script D01
# ==============================================================================



# usually run interactively, ~12 hours -----------------------------------------

screen -S combp
qlogin -R rusage[mem=16]

module unload anaconda
conda activate py38
module load bedtools/2.26.0

export PATH=$HOME/Software/combined-pvalues/:$PATH

cd $HOME/Analyses/DRIFT2_Methylation/Output/20220425/combp
mkdir results




# baseline only ----------------------------------------------------------------

cat > run_combp_BLonly.sh
for filename in bedfiles/BLonly*.bed; do
out=$(basename "${filename%.bed}")
out="results/${out}"
comb-p pipeline -c 4 --seed 0.05 --dist 500 -p $out --region-filter-n 3 $filename
done

sh run_combp_BLonly.sh &> combp_BLonly.log



# baseline --> 3mo or 6mo delta outcomes ---------------------------------------

cat > run_combp_shortterm.sh

for filename in bedfiles/Methyl0*_30_*.bed; do
out=$(basename "${filename%.bed}")
out="results/${out}"
comb-p pipeline -c 4 --seed 0.05 --dist 500 -p $out --region-filter-n 3 $filename
done

for filename in bedfiles/Methyl0*_60_*.bed; do
out=$(basename "${filename%.bed}")
out="results/${out}"
comb-p pipeline -c 4 --seed 0.05 --dist 500 -p $out --region-filter-n 3 $filename
done

sh run_combp_shortterm.sh &> combp_shortterm.log




# time or trt --> delta 3mo methylation ----------------------------------------

cat > run_combp_trtandtime.sh
for filename in bedfiles/Trt*.bed; do
out=$(basename "${filename%.bed}")
out="results/${out}"
comb-p pipeline -c 4 --seed 0.05 --dist 500 -p $out --region-filter-n 3 $filename
done

for filename in bedfiles/Time*.bed; do
out=$(basename "${filename%.bed}")
out="results/${out}"
comb-p pipeline -c 4 --seed 0.05 --dist 500 -p $out --region-filter-n 3 $filename
done

sh run_combp_trtandtime.sh &> combp_trtandtime.log



