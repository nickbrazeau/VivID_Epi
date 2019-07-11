#!/usr/bin/env bash

vcf=/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/WetLabWork/PvampliconDesign/vcf/snps.VQSR900.human_only.vcf.gz
bed=/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/WetLabWork/PvampliconDesign/primerbed_fastas/pv_p01_prop-primers.bed
out=/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/WetLabWork/PvampliconDesign/vcf/prop-primers_snps.VQSR900.human_only.vcf.gz

bedtools intersect -a $vcf -b $bed -header | bgzip > $out
