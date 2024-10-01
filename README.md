# Performing pedigree inference in order to develop replacement population index for Michigan’s Upper Peninsula Black Bears

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

[Input data files](https://github.com/mweissman97/MI_bear_pedigree/tree/100185359ba796c4bf1fad4b714d74868d14afa2/input_data_files):

| 21/22 File Name                                | Description                                                                                                                                 |
|----------------------------------------|--------------------------------|
| GT_BB2122_nuc.csv                              | CSV with nDNA snps. Columns correspond to SNPs. 0 = no copies of variant, 1 = heterozygous, 2 = homozygous with variant, -9 = missing read. |
| 2021_CKMR_Sample_Data_age.xlsx                 | Excel sheet with life history data for samples collected in 2021                                                                            |
| 2022_CKMR_Sample_Data_age.xlsx                 | Excel sheet with life history data for samples collected in 2022                                                                            |
| BB2122_compiled_R1_genotypes.csv               | CSV of raw nDNA. Includes meta data (i.e. IFI) and sequencing of 5 genetic sex markers                                                      |
| BB_2022_001_calls_annotated_maf0.1_sequoia.vcf | VCF (variant call format) file with mitochondrial DNA for all individuals                                                                   |

[R code files](https://github.com/mweissman97/MI_bear_pedigree/tree/100185359ba796c4bf1fad4b714d74868d14afa2/r_code_files):

| R File                          | Description                                                                              |
|------------------------------------|------------------------------------|
| bb_functions_datafiltering.R    | Custom built data filtering functions                                                    |
| bb_nDNA_qcfilter.R              | Filters samples based on nuclear DNA                                                     |
| bb_mtDNA_qcfilter.R             | Filters and analyzes mitochondrial DNA                                                   |
| bb_functions_popgen.R           | Custom built functions for popgen analyses                                               |
| bb_popgen_nuc.R                 | Performs popgen analyses on nuclear DNA                                                  |
| bb_mtDNA_popgen.R               | Performs popgen analyses on mitochondrial DNA                                            |
| bb_mDNA_haplotypes.R.           | Assigns mitochondrial haplotypes and analyzes haplotype network. *not working right now* |
| bb_pedigree_diagnostics.R       | Uses CKMRsim package to test inference power of nuclear DNA in resolving pedigree        |
| bb_functions_pedigree_pairs.R   | Custom built functions for analyzing sequoia pedigrees                                   |
| bb_sequoia_pedigree_inference.R | Creates pedigree using sequoia package                                                   |

[Output Files](https://github.com/mweissman97/MI_bear_pedigree/tree/100185359ba796c4bf1fad4b714d74868d14afa2/output_files):

| Output File                 | Description                                                                                                          |
|------------------------------------|------------------------------------|
| QCBears_fulldata_final.csv  | Pruned QC bears, snps, and life history data                                                                         |
| AllBears_fulldata_final.csv | All bears, with SNPS and life history                                                                                |
| QCBears_mDNA.csv            | Pruned mitochondrial DNA snps for QC bears                                                                           |
| BB_POpairs_final.csv        | CSV of all PO pairs, with LH data. **File used in CKMR model.**                                                      |

<!-- CONTACT -->
## Contact <a name="contact"></a>

[Personal Website](https://sciencemaya.com) - mweissman97@gmail.com
