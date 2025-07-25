Notes on the reference.fasta

This reference file is a modified version of the reference.fasta downloaded from the IGenotyper github page (https://github.com/oscarlr/IGenotyper) November 2023.

Modifications to light chain loci:

chr22 was exchanged for T2T (CHM13v2.0) chr22. The IGL locus is on chromosome 22. ~28.3 kb of sequence (chr22:23,317,584-23,345,947) is from a HPRC haplotype from sample HG00621, which includes a 16093 bp insertion relative to the CHM13v2.0 reference, reflecting 3 additional copies of IGLJ3-C3. To lift-over a chr22 coordinate from CHM13v2.0 to the chr22 in the reference.fasta in this repo, use the formula: "if start_coord > 23324706, then add 16093".

chr2, which includes the IGK locus, was modified to include a common ~24 kb insertion in the distal IGK region that includes the gene "IGKV1-NL1" - described in our publication: https://www.nature.com/articles/s41435-024-00279-2 To lift-over from hg38 to the IGK locus in reference.fasta in this repo, use the formula "if start_coord > 90225183, add 24729"

Modifications to heavy chain loci (igh, ighc): 

No modifications to igh were made in this repo relative to the igh reference at https://github.com/oscarlr/IGenotyper. For a description of this igh reference, see our publications: i) https://www.frontiersin.org/articles/10.3389/fimmu.2020.02136/full ii) https://www.nature.com/articles/s41467-023-40070-x which describe the igh franken reference.

The ighc reference sequence is from CHM13v2.0, without modification.

#to-do: add ancestry-informative markers BED file and "non-IGL_regions_to_mask" BED file
