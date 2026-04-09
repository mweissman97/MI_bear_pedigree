# Performing pedigree inference on Michigan’s Upper Peninsula Black Bears

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About This Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
        <li><a href="#files">File List</a></li>
      </ul>
    </li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## Abstract <a name="about-the-project"></a>
The goal of this project was to establish a Close-Kin-Mark-Recapture (CKMR) framework for Michigan’s Upper Peninsula (UP) black bear population. Key to this endeavor was the development of a streamlined genotyping method capable of confidently assigning familial relationships of harvested bears, from which DNA samples were collected. This bioinformatic pipeline uses the R Programming Language to filter both nuclear and mitochondrial SNP data, determinesthe panel’s confidence in assigning relatedness amongst these data, and, finally, construct a pedigree to output a list of all Parent-Offspring pairs.

This work was conducted with the [Brzeski Lab](https://www.mtu.edu/forest/about/faculty-staff/faculty/brzeski/) at Michigan Technological University, and in partnership with the Michigan [Department of Natural Resources](https://www.michigan.gov/dnr).

For more, see the [R Markdown Guide](https://github.com/mweissman97/MI_bear_pedigree/blob/c679d81282330d025a48b15a3a953efbd045f683/bb_kinship_inference-compressed.pdf).

### Built With <a name="built-with"></a>

Code is written using the following languages and packages:
* [R Programming Language](https://www.r-project.org/)
  * [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
  * [tidyr](https://tidyr.tidyverse.org/)
  * [ggplot2](https://ggplot2.tidyverse.org/)
  * [viridis](https://github.com/eriqande/CKMRsim)
  * [readxl](https://readxl.tidyverse.org/)
  * [CKMRsim](https://github.com/eriqande/CKMRsim)
  * [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html)
  * [hierfstat](https://cran.r-project.org/web/packages/hierfstat/index.html)
  * [pegas](https://cran.r-project.org/web/packages/pegas/index.html)
  * [smartsnp](https://cran.r-project.org/web/packages/smartsnp/index.html)
  * [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html)
  * [ape](https://cran.r-project.org/web/packages/ape/index.html)
  * [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html)
  * [sequoia](https://cran.r-project.org/web/packages/sequoia/index.html)

### File List <a name="files"></a>

[R code files](https://github.com/mweissman97/MI_bear_pedigree/tree/100185359ba796c4bf1fad4b714d74868d14afa2/r_code_files):

| R File Name                                | Description                                                                                                                                 |
|----------------------------------------|--------------------------------|
| bb_functions_datafiltering.R    | Custom built data filtering functions                                                    |
| bb_nDNA_qcfilter.R              | Filters samples based on nuclear DNA                                                     |
| bb_mtDNA_qcfilter.R             | Filters and analyzes mitochondrial DNA                                                   |
| bb_functions_popgen.R           | Custom built functions for popgen analyses                                               |
| bb_popgen_nuc.R                 | Performs popgen analyses on nuclear DNA                                                  |
| bb_mtDNA_popgen.R               | Performs popgen analyses on mitochondrial DNA                                            |
| bb_pedigree_diagnostics.R       | Uses CKMRsim package to test inference power of nuclear DNA in resolving pedigree        |
| bb_functions_pedigree_pairs.R   | Custom built functions for analyzing sequoia pedigrees                                   |
| bb_sequoia_pedigree_inference.R | Creates pedigree using sequoia package                                                   |

[Input data files](https://github.com/mweissman97/MI_bear_pedigree/tree/100185359ba796c4bf1fad4b714d74868d14afa2/input_data_files):
* **GT_BB2122_nuc.csv** - CSV with nDNA snps.
  * Rownames - long sample IDs (e.g. BB21_1000_MTU20Oct23_R1)
  * Columns NW_025576331_1_3826608:NW_025578505_1_8669226 - SNP data, where 0 = no copies of variant, 1 = heterozygous, 2 = homozygous with variant, -9 = missing read.
* **2122_CKMR_Sample_Data.csv** - table with life history data for samples collected in 2021 and 2022
  * Sequoia_ID - transformed sample ID to use underscores to match GT_BB2122_nuc.csv and be compatible with the sequoia package (e.g. BB21_1000)
  * Season_Year - year sample was harvested
  * Species - common name for species (i.e. "Black Bear")
  * Harvest_Date - date sample was harvested (e.g. 8-Sep-21)
  * How_Taken - how sample was harvested (e.g. "Hunting")
  * Harvest_BMU - Bear Management Unit where sample was harvested
  * Harvest_County - County where sample was harvested
  * Registration_Sex - sex identified during registration, may differ from genetic sex
  * BirthYear - estimated year the bear was born based on sample age
  * Latitude - approximate decimal latitude coordinates where sample was harvested
  * Longitude - approximate decimal longitude coordinates where sample was harvested
* **BB2122_compiled_R1_genotypes.csv** - CSV of raw nuclear DNA. Includes meta data (i.e. IFI) and sequencing of 5 genetic sex markers
  * Sample - long sample IDs (e.g. BB21_1000_MTU20Oct23_R1)
  * Raw Reads - number of raw reads during genotyping
  * On-Target Reads - number of reads on target 
  * %On-Target - On-Target Reads divided by Raw Reads
  * %GT - percent of targeted loci sequenced
  * IFI
  * NW_025576331.1_3826608:NW_025578505.1_8669226 - unphased genotype matrix, where missing data = 0
  * Uam_SEXY1:Uam_sry2 - genotypes for sex markers; XX = female, XY = male, 0 = missing
* **BB_2022_001_calls_annotated_maf0.1_sequoia.vcf** - VCF (variant call format) file with mitochondrial DNA for all individuals

[Output Files](https://github.com/mweissman97/MI_bear_pedigree/tree/100185359ba796c4bf1fad4b714d74868d14afa2/output_files):
* **AllBears_fulldata_final.csv** and **QCBears_fulldata_final.csv** - merged data for samples, SNPs, and life history data. **QCBears_fulldata_final.csv** is only the subset of samples that passed filtering steps (i.e. filter_remove == "pass")
  * Sequoia_ID - shortened sample ID (e.g. BB21_1000)
  * Sample - long sample IDs (e.g. BB21_1000_MTU20Oct23_R1)
  * NW_025576331.1_3826608:NW_025578505.1_8669226 - SNP data, where 0 = no copies of variant, 1 = heterozygous, 2 = homozygous with variant, -9 = missing read.
  * Uam_SEXY1:Uam_sry2 - genotypes for sex markers; XX = female, XY = male, 0 = missing
  * Genetic_Sex - most common sex from sex marker columns (Uam_SEXY1:Uam_sry2)
  * sex_agreement_prop - fraction of sex markers (Uam_SEXY1:Uam_sry2) that agree with Genetic_Sex consensus
  * Season_Year - year sample was harvested
  * Species - common name for species (i.e. "Black Bear")
  * Harvest_Date - date sample was harvested (e.g. 8-Sep-21)
  * How_Taken - how sample was harvested (e.g. "Hunting")
  * Harvest_BMU - Bear Management Unit where sample was harvested
  * Harvest_County - County where sample was harvested
  * Registration_Sex - sex identified during registration, may differ from genetic sex
  * BirthYear - estimated year the bear was born based on sample age
  * Latitude - approximate decimal latitude coordinates where sample was harvested
  * Longitude - approximate decimal longitude coordinates where sample was harvested
  * filter_remove - whether the sample passed quality control filters (i.e. "pass" or "fail")            
* **QCBears_mDNA.csv** - Pruned mitochondrial DNA snps for QC bears
  * NC_003426.1_107:NC_003426.1_16604 - genotype matrix for mitochondrial DNA, where 0 = wild-type, 1 = variant, NA = missing
* **BB_POpairs_final.csv** - csv of all parent offspring pairs, with life-history data. *File used in CKMR model.*
  * bear1 - long sample ID corresponding to the younger bear (aka offspring) in the relationship pair (e.g. BB21_1000_MTU20Oct23_R1)
  * bear2 - long sample ID corresponding to the older bear (aka parent) in the relationship pair (e.g. BB21_1000_MTU20Oct23_R1)
  * LLR - log likelihood ratio, or log10 transformed likelihood the pair have the assigned relationship divided by the likehlihood the pair have the next most likely relationship type; higher values of LLR correspond to higher confidence
  * rel_type - relationship type; PO_M = mother-offspring, PO_P = father-offspring
  * pair_id - concatenates bear1 ID and bear2 ID to create a unique identifier for the pair
  * bear1_birthyear - birth year of bear1
  * bear2_birthyear - birth year of bear2
  * age_diff - age of parent at offspring's birth, bear1_birthyear - bear2_birthyear
  * bear1_sex - sex of bear1
  * bear2_sex - sex of bear2
  * bear1_county - harvest county of bear1
  * bear2_county - harvest county of bear2
  * county_difference - spatial relationship between bear1_county and bear2_county; "same county" = harvested in the same county, "neighboring county" = bear1_county borders bear2_county, "further" = bears were not harvested in the same or bordering counties
  * module - sequoia model used; "full ped" = full pedigree model, rather than parent-offspring only pedigree
  * age_prior - whether sequoia model incorporated age priors; "yes ap" = age priors were used
  * LLR_bin - turns LLR into a discrete bin; "negative" = LLR was less than 0 and thus there is low confidence in relationship assignment, ">0" = LLR was greater than 0 and thus there is high confidence in relationship assignment

<!-- CONTACT -->
## Contact <a name="contact"></a>

[Personal Website](https://sciencemaya.com) - mweissman97@gmail.com
