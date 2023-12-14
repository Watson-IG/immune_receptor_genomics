#!/bin/bash
set -e -x

function run_make_chr2_fasta {
    awk '/^>chr2($|[^0-9])/{p=1;print;next} /^>/{p=0} p' /home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/reference_for_VDJbase/reference.fasta > /home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr2.fasta
    samtools faidx /home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr2.fasta
}

function run_get_IGK_alleles {
    IGK_coords=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/IGK_coords
    mgf=/home/egenge01/VDJbase-Genomics/python/make_gene_file_EEmod.py
    ifa=/home/egenge01/VDJbase-Genomics/python/import_from_assemblies.py
    mkdir -p $PWD/IGK_alleles
    chr2_fasta=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr2.fasta
    cat filtered_2023-10-25_fofn.txt | cut -f1 | grep -v '2209419009' | grep -v '191181105C' | \
	grep -v '220672802C' | grep -v '220683005C' | grep -v '2208409001' | grep -v '220872901C' | grep -v '2202422009' | while read sample
    do
	mkdir -p $PWD/IGK_alleles/${sample}
	# conda activate vdjb-genomics
	bamf=/home/egenge01/projects/CW50/curated_samples/merged_alignments/${sample}/${sample}_merged.sorted.bam
	python ${mgf} ${bamf} IGK chr2 +- ${IGK_coords} ${chr2_fasta} $PWD/IGK_alleles/${sample}/all_genes.csv  
    done
}

function run_get_IGK_alleles_p2 {
    IGK_coords=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/IGK_coords
    ifa=/home/egenge01/VDJbase-Genomics/python/import_from_assemblies.py
    seqs_ref=/home/egenge01/VDJbase-Genomics/ref 
    chr2_fasta=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr2.fasta
    cat filtered_2023-10-25_fofn.txt | cut -f1 | grep -v '2209419009' | grep -v '191181105C' | \
	grep -v '220672802C' | grep -v '220683005C' | grep -v '2208409001' | grep -v '220872901C' | grep -v '2202422009' | while read sample
    do
	python ${ifa} IGK chr2 +- $PWD/IGK_alleles/${sample}/all_genes.csv \
	    ${chr2_fasta} ${IGK_coords} ${seqs_ref} $PWD/IGK_alleles/${sample}/annotations.csv 
    done
}

#run_make_chr2_fasta
#run_get_IGK_alleles
run_get_IGK_alleles_p2
