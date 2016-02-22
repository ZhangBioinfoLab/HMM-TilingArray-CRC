#!/bin/sh

#argument1: feature file with name in bed format, 
#argument2: output prefix

python ~/workspace_git/working_repo/cirDNA_project/scripts/analysis/addNameToBed.py $1 $1'temp'

bedtools intersect -a $1'temp' -b /home/zhanglab/junjianglin/project_data/cir_cancer/thesis_project/Bio_analysis/genomic_annotation/hg18_exons_uniq.txt -wa -wb > exon.txt
bedtools intersect -a $1'temp' -b /home/zhanglab/junjianglin/project_data/cir_cancer/thesis_project/Bio_analysis/genomic_annotation/hg18_introns_uniq.txt -wa -wb > introns.txt
bedtools intersect -a $1'temp' -b /home/zhanglab/junjianglin/project_data/cir_cancer/thesis_project/Bio_analysis/genomic_annotation/hg18_promotor_uniq.txt -wa -wb > promotor.txt
bedtools intersect -a $1'temp' -b /home/zhanglab/junjianglin/project_data/cir_cancer/thesis_project/Bio_analysis/genomic_annotation/hg18_all_uniq.txt -v -wa -wb > intergenic.txt

python ~/workspace_git/working_repo/cirDNA_project/scripts/analysis/genomic_annotation.py $1 $2 -e exon.txt -i introns.txt -g intergenic.txt -p promotor.txt

rm exon.txt introns.txt intergenic.txt promotor.txt $1'temp'
