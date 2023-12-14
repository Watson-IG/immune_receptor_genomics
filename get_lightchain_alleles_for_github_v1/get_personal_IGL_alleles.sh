#!/bin/bash
set -e -x

function run_make_chr22_fasta {
    awk '/^>chr22($|[^0-9])/{p=1;print;next} /^>/{p=0} p' /home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/reference_for_VDJbase/reference.fasta > /home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr22.fasta
    samtools faidx /home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr22.fasta
}

function run_get_IGL_alleles {
    IGL_coords=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/IGL_coords
    mgf=/home/egenge01/VDJbase-Genomics/python/make_gene_file_EEmod.py
    ifa=/home/egenge01/VDJbase-Genomics/python/import_from_assemblies.py
    mkdir -p $PWD/IGL_alleles
    chr22_fasta=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr22.fasta
    cat /home/egenge01/projects/CW50/IGL_curated_contigs/IGL_samples.txt | while read sample 
    do
	mkdir -p $PWD/IGL_alleles/${sample}
	# conda activate vdjb-genomics
	bamf=/home/egenge01/projects/CW50/IGL/curated_samples/merged_alignments/${sample}/${sample}_merged.sorted.bam
	python ${mgf} ${bamf} IGL chr22 + ${IGL_coords} ${chr22_fasta} $PWD/IGL_alleles/${sample}/all_genes.csv  
    done
}

function run_get_IGL_alleles_p2 {
    IGL_coords=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/IGL_coords
    ifa=/home/egenge01/VDJbase-Genomics/python/import_from_assemblies.py
    seqs_ref=/home/egenge01/VDJbase-Genomics/ref
    chr22_fasta=/home/egenge01/projects/IGL_ref_mod/files_for_VDJbase/putting_together/annotations_freeze1/chr22.fasta
    cat /home/egenge01/projects/CW50/IGL_curated_contigs/IGL_samples.txt | while read sample
    do
	python ${ifa} IGL chr22 + $PWD/IGL_alleles/${sample}/all_genes.csv \
	    ${chr22_fasta} ${IGL_coords} ${seqs_ref} $PWD/IGL_alleles/${sample}/annotations.csv 
    done
}

#run_make_chr22_fasta
#run_get_IGL_alleles
run_get_IGL_alleles_p2
