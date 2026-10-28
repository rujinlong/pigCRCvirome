# Processed inputs

These files are fixed inputs to the [R analyses](../../../README.md). The qproj directory name `00-raw` denotes the read-only input layer; the data have already been processed. SHA-256 hashes and file sizes are listed in [checksums.tsv](checksums.tsv).

| File | Contents |
|---|---|
| [samples.tsv.gz](samples.tsv.gz) | 73 libraries: analysis sample ID, virome library ID, pig ID, diet, visit, site, sex, and litter |
| [clinical.tsv.gz](clinical.tsv.gz) | 16 endpoint pigs: sample and pig IDs, diet, pre-intervention distal-colon polyp count, and recorded progressing-polyp count |
| [counts.rds](counts.rds) | Numeric matrix of mapped counts: 16,934 vOTUs by 73 libraries; row names are vOTU IDs and column names are analysis sample IDs |
| [endpoint_tpm.rds](endpoint_tpm.rds) | TPM matrix: the same 16,934 vOTUs by 16 endpoint fecal samples |
| [votu_annotation.tsv.gz](votu_annotation.tsv.gz) | Catalog length, CheckV completeness, geNomad and VirSorter2 scores, predicted lifestyle, vConTACT3 group, and assembly library |
| [linda_coefficients.tsv.gz](linda_coefficients.tsv.gz) | Diet-versus-CTRL log2 fold changes and unadjusted P values from the primary and depth-adjusted LinDA fits in feces and gut contents |
| [host_family_membership.tsv.gz](host_family_membership.tsv.gz) | Unique vOTU-family links from genome-level iPHoP predictions with confidence at least 75 |
| [family_permutation_pvalues.tsv.gz](family_permutation_pvalues.tsv.gz) | Original permutation P values and permutation counts for all 279 eligible host families |
| [desulfo_genera.tsv.gz](desulfo_genera.tsv.gz) | Highest-confidence host genus for each of the 69 Desulfovibrionaceae-linked vOTUs |
| [bai_candidates.tsv.gz](bai_candidates.tsv.gz) | 167 candidate phages linked to 16 reference genomes with qualifying baiE neighborhoods; host predictions, sequence quality, and lifestyle calls |
| [bai_hosts.tsv.gz](bai_hosts.tsv.gz) | The 16 qualifying reference genomes, host families, and baiE-neighborhood flag |
| [candidate_detection.tsv.gz](candidate_detection.tsv.gz) | TPM and covered fraction for the two illustrated candidates in the 45 available fecal libraries |
| [candidate_genes.tsv.gz](candidate_genes.tsv.gz) | 150 predicted genes: vOTU and protein IDs, gene ordinal, nucleotide coordinates, strand, functional class, and highlighted label |
| [baseline_covered_fraction.rds](baseline_covered_fraction.rds) | Covered fractions for the 33 pre-intervention-origin Desulfovibrionaceae vOTUs in 16 endpoint fecal samples |
| [endpoint_depth.tsv.gz](endpoint_depth.tsv.gz) | Mapped and total read counts for the 16 endpoint fecal libraries |
| [permutation_states.rds](permutation_states.rds) | Five R random states for the community, paired-change, diet, pairwise abundance, and catalog-origin tests |
| [submitted_source_data.xlsx](submitted_source_data.xlsx) | The submitted workbook, with Contents and D1-D19; used only to check the regenerated results |

## Identifiers and units

`sample_id` joins sample metadata, clinical data, depth, and matrix columns. `id_virome` identifies the corresponding sequencing library; `host_id` identifies the pig. `Month_raw` preserves the input visit labels: Aug, Sep, and Dec correspond to study months 3, 4, and 7. `Location` is fecal, cecum, or colon. Diet codes are CTRL (control feed), RL (red meat and lard), and COL (RL plus colestyramine). Missing litter values remain NA.

Counts are mapped read counts. TPM is length-normalized read abundance scaled per million. Set abundance is summed member TPM divided by total catalog TPM; plots express this fraction as a percentage. Covered fraction ranges from 0 to 1. Candidate detection requires TPM > 0 and covered fraction >= 0.20. Gene coordinates are 1-based nucleotide positions; strand is 1 (forward) or -1 (reverse). CheckV completeness is a percentage.

`assembly_library` records the source of a representative sequence. Detection and abundance come from read mapping. For the catalog-origin comparison, the denominator is the 10,109 vOTUs assembled from month-3 or month-4 feces, including 33 Desulfovibrionaceae-linked members.

## Processed statistical and annotation results

LinDA v0.2.0 used the 13,782 vOTUs with geNomad scores >= 0.8 and VirSorter2 scores >= 0.9, retaining non-zero counts in at least two samples per compartment. Models included diet and sex; gut-content models also included site and a pig random intercept. Depth-adjusted models added log10 mapped counts. Settings were prevalence cutoff 0.05, library cutoff 1, winsorization quantile 0.97, and default adaptive zero handling. The public analysis applies BH correction jointly across both diet coefficients within each fit and defines diet-associated vOTUs at q < 0.10.

Family sets may overlap when a vOTU has several predicted hosts. Of 443 represented host families, 279 have non-zero TPM in at least eight endpoint pigs. Their supplied P values come from 4,999 Freedman-Lane residual permutations within diet, adjusting for diet and ranked baseline polyp count. The analysis recalculates each partial rank correlation, its approximate Fisher-z interval, and BH correction over the 279 P values.

The baiE-neighborhood flag summarizes the host-proteome screen: baiE, distinct baiCD and baiH homologs, and at least two auxiliary enzyme profile types on the same contig within 10 kb. The 16 qualifying genomes are the highest-confidence predicted hosts of the 167 candidate phages. Gene functional classes and highlighted annotations are processed DRAM-v and PHROG assignments. Predicted temperate status is the union of VIBRANT lysogenic, geNomad provirus, CheckV provirus, and BACPHLIP temperate calls.
