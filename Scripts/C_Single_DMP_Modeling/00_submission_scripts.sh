# ==============================================================================
# C/00_submission_scripts.sh
# helper scripts to submit regression scripts to the scheduler
# ==============================================================================


# in directory to hold scripts
# "cat > [filename].sh" below creates the submission script

# also performed but not shown here: unfinalized model exploration /
# e.g., model goodness of fit checks, alt modeling frameworks (e.g., linear
# mixed model approaches, beta vals), covariate sensitivity analyses 

mkdir $HOME/Analyses/DRIFT2_Methylation/SubmissionScripts
cd $HOME/Analyses/DRIFT2_Methylation/SubmissionScripts



# test if methylation changes btwn timepoint / btwn intervention ---------------

screen
qlogin -R rusage[mem=12]
cd $HOME/Analyses/DRIFT2_Methylation
module load R/4.1.0_beta
R

# then ran 01_intervention_and_time.R interactively
# (fast enough / small number of modules that don't need batching below)



# BL --> short term cardiometabolic outcomes -----------------------------------
# e.g., baseline methylation to 6-mo change in insulin
# run in batches (e.g., batch 1 to 15)

cat > run_models_BLMe_to_ShortTerm.sh
#!/usr/bin/env bash
#BSUB -J BLtoST_placeholderindex_placeholderdate
#BSUB -q normal
#BSUB -R "select[mem>2] rusage[mem=5] span[hosts=1]"
#BSUB -W 48:00 
#BSUB -o BLtoST_placeholderindex_placeholderdate.out 
#BSUB -e BLtoST_placeholderindex_placeholderdate.errout

cd $HOME/Analyses/DRIFT2_Methylation
module load R/4.1.0_beta
R CMD BATCH --no-save '--args placeholderindex' \
Scripts/C/02_BLMethyl_to_ShortTermOutcomes.R \
SubmissionScripts/BLtoST_placeholderindex_placeholderdate.Rout

# note that the end of seq should match nbatch in the script
sed -i 's/placeholderdate/20220425/g' run_models_BLMe_to_ShortTerm.sh
for i in `seq 1 15`; do
echo $(sed s/placeholderindex/$i/g < run_models_BLMe_to_ShortTerm.sh | bsub);
done




# baseline methylation <--> outcomes --------------------------------------------
# note: ultimately not reported in primary DRIFT2 DNAme manuscript,
# since larger sample size, large profiles of common measures e.g., BMI exist
# however, run for rare outcomes (e.g., food/exercise motivation)

cat > run_BLonly_models.sh
#!/usr/bin/env bash
#BSUB -J BLonly_placeholderindex_placeholderdate
#BSUB -q normal
#BSUB -R "select[mem>12] rusage[mem=12] span[hosts=1]"
#BSUB -W 48:00 
#BSUB -o BLonly_placeholderindex_placeholderdate.out 
#BSUB -e BLonly_placeholderindex_placeholderdate.errout

cd $HOME/Analyses/DRIFT2_Methylation
module load R/4.1.0_beta
R CMD BATCH --no-save '--args placeholderindex' \
Scripts/C/03_baseline_only.R \
SubmissionScripts/BLonly_placeholderindex_placeholderdate.Rout

# need end of sequence to match nbatch in the script
sed -i 's/placeholderdate/20220425/g'  run_BLonly_models.sh
for i in `seq 1 10`; do
echo $(sed s/placeholderindex/$i/g < run_BLonly_models.sh | bsub);
done

