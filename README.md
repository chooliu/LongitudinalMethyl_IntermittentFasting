# LongitudinalMethyl_IntermittentFasting 

## Summary

Testing associations between longitudinal measurements of methylation and participant weight/cardiometabolic outcomes in the Daily Caloric Restriction vs. Intermittent Fasting Trial (DRIFT).

We show that **(i) trial participants' baseline, pre-intervention methylation profiles in PBMC are associated with their 3-month and 6-month weight and cardiometabolic outcomes**, and that **(ii) intermittent fasting is associated with different 3-month changes in DNA methylation levels than daily caloric restriction**, controlling for cell composition, ancestry (measured via genotyping array), baseline BMI, age, sex, and technical confounders.


## Files in Repo

```
File Structure
└─ Scripts
│   └─ A_Preprocessing/             # process methylation & participant metadata
│   └─ B_EDA_Descriptive_Stats/     # exploratory data analysis: demographics, outcomes
│   └─ C_Single_DMP_Modeling/       # univariate testing of CpGs
│   └─ D_Multiprobe_Modeling/       # multivariate testing (DMRs, machine learning models)
└─ README.md                        # this readme
└─ sessionInfo.txt                  # version control
```

**Code Sections:**

* **A_Preprocessing**: processing of methylation & genotype array data (incl. ancestry, cell composition, RUVg inference), tidying clinical metadata 
* **B_EDA_Descriptive_Stats**: exploratory data analysis
* **C_Single_DMP_Modeling**: single CpG at a time testing of baseline methylation <--> Δoutcomes, and intervention <--> Δmethylation regression modeling / visualization (where Δ = 3mo - baseline or 6mo - baseline changes)
* **D_Multiprobe_Modeling**: differentially methylated region calling, multiple CpG at a time prediction models for feature selection and association estimates (e.g., variance in Δinsulin levels explained by baseline methylation)

## Links

**Manuscript Links**
  
* This Github Repo: [chooliu/LongitudinalMethylation_IntermittentFasting](http://www.github.com/chooliu/LongitudinalMethylation_IntermittentFasting)
* Manuscript Link: TBD

**Data Sharing**

Unfortunately, we cannot share individual-level clinical metadata or methylation measurements due to IRB consent, but can share aggregated information. Hoping to share summary statistics in tidy computer-readable format.

## Related Publications

Primary manuscript under preparation, but some preview abstracts below (no substantive code changes since these abstracts)

1. Borengasser S, Liu C, Kechris-Mays K, Siebert J, Stanislawski M, Litowski E, Ostendorf D, Bing K, Wayland L, Scorsone J, Maclean P, Melanson E, Bessesen D, Catenacci V. Baseline DNA Methylation Predicts Changes in Cardiometabolic Health Post Behavioral Weight Loss. Obesity. 2022 Nov 1;30:235-6.

2. Shogan L, Borengasser S, Cooper E, Liu C, Kechris-Mays K, Lafata E, Heinsbroek J, Stanislawski M, Creasy S, Catenacci V, Ostendorf D. DNA Methylation Is Associated With Exercise-Related Factors Among Adults With Overweight/Obesity. Obesity. 2022 Nov 1;30:235-6.


**Other Related Work with DRIFT2 Methylation Data**

1. Siebert JC, Stanislawski MA, Zaman A, Ostendorf DM, Konigsberg IR, Jambal P, Ir D, Bing K, Wayland L, Scorsone JJ, Lozupone CA. Multiomic Predictors of Short‐Term Weight Loss and Clinical Outcomes During a Behavioral‐Based Weight Loss Intervention. Obesity. 2021 May;29(5):859-69. [PMID:33811477](https://pubmed.ncbi.nlm.nih.gov/33811477/), [doi:10.1002/oby.23127](https://doi.org/10.1002/oby.23127)

2. Ostendorf, Danielle M., et al. "Comparison of weight loss induced by daily caloric restriction versus intermittent fasting (DRIFT) in individuals with obesity: study protocol for a 52-week randomized clinical trial." Trials 23.1 (2022): 718. [PMID: 36038881](https://pubmed.ncbi.nlm.nih.gov/36038881/),  [doi:10.1186/s13063-022-06523-2](https://doi.org/10.1186/s13063-022-06523-2)

3. Hill EB, Konigsberg IR, Ir D, Frank DN, Jambal P, Litkowski EM, Lange EM, Lange LA, Ostendorf DM, Scorsone JJ, Wayland L. The Microbiome, Epigenome, and Diet in Adults with Obesity during Behavioral Weight Loss. Nutrients. 2023 Aug 16;15(16):3588. [PMID:37630778](https://pubmed.ncbi.nlm.nih.gov/37630778/), [doi:10.3390/nu15163588](https://doi.org/10.3390/nu15163588)
