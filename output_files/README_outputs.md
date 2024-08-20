| Output File                 | Created In                      | Description                                                                                                          |
|-------------------|-----------------|------------------------------------|
| QCBears_fulldata_final.csv  | bb_nDNA_qcfilter.R              | Pruned QC bears, snps, and life history data                                                                         |
| AllBears_fulldata_final.csv | bb_nDNA_qcfilter.R              | All bears, with SNPS and life history                                                                                |
| nuclear_popgen_clusters.csv | bb_popgen_nuc.R                 | QC bear dataframe, plus columns for PC and k-means cluster data based on nuclear DNA                                 |
| nuc_gdist_fst.csv           | bb_popgen_nuc.R                 | CSV of genetic distance (Fst) between Harvest BMUs                                                                   |
| nuc_popgen_ibd.csv          | bb_popgen_nuc.R                 | CSV of pairwise genetic distances and geographic distances for all individuals.                                      |
| QCBears_mDNA.csv            | bb_mtDNA_qcfilter.R             | Pruned mitochondrial DNA snps for QC bears                                                                           |
| mito_popgen_clusters.csv    | bb_popgen_nuc.R                 | QC bear dataframe, plus columns for PC and k-means cluster data based on mitochondrial DNA                           |
| CKMRsim_error_rates.csv     | bb_pedigree_diagnostics.R       | False positive and false negative error rates found using simulated genotypes to test inference power of nuclear DNA |
| seqped_ap_po.rds            | bb_sequoia_pedigree_inference.R | R Data object of full pedigree, with age priors and PO pairs only                                                    |
| seqped_ap_full.rds          | bb_sequoia_pedigree_inference.R | R Data object of full pedigree, with age priors and full pedigree                                                    |
| seqped_noap_po.rds          | bb_sequoia_pedigree_inference.R | R Data object of full pedigree, with no age priors and PO pairs only                                                 |
| seqped_noap_full.rds        | bb_sequoia_pedigree_inference.R | R Data object of full pedigree, with no age priors and full pedigree                                                 |
| BB_POpairs_final.csv        | bb_sequoia_pedigree_inference.R | CSV of all PO pairs, with LH data. **File used in CKMR model.**                                                      |
