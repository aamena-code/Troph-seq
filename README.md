# Mouse trophoblast sequencing
Analysis of sequencing data from mouse trophoblasts

This work is comprised of bulk RNA-seq and ATAC-seq datasets. 
We also compare our RNA-seq data to single-cell and single-nuclei RNA-seq data taken from published work by Jiang _et al.,_ 2023 and Marsh and Blelloch _et al.,_ 2020.

Comparison of ATAC-seq with histone marks of regulatory elements in TSCs and day 4 differentiated trophoblasts was carried out in seqMINER using k-means normalised linear enrichment. We used the application and, therefore, no script was required. This ChIP-seq data was taken from Chuong _et al.,_ 2013.

Following alignment, and generation of counts files, we carried out subsequent RNA-seq analysis in R v4.3.1. Similarly, following initial processing of raw files, ATAC-seq peak calling and annotation was carried out in R v4.3.1.

Experimental work:
We treated mouse TSCs _in vitro_ with PROTAC degrader ACBI1, or bromodomain inhibitor PFI-3 and used _cis_-ACBI1 and DMSO as control treatments. 
Treatments were added to media at day -1 of differentiation and RNA harvested from trophoblasts at day 0 (in TSC stage), day 2 and day 4 of differentiation. ATAC-seq was carried out in parallel in day 0 and day 2 differentiated trophoblasts. 
Activin A was also added to differentiation media of half of the day 2 and 4 differentiated trophoblasts, in combination with ACBI1 or PFI-3 treatments, to determine whether BAF complex dysfunction affected response to Activin A.

For any further questions, please contact me directly at: aamenapatel@live.co.uk
