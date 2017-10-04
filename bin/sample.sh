#!/usr/bin/bash
#
# sample.sh - silly little script to create a file containing a sampling of MGI.gff3
#
# Basic categories:
# MGI:3528744 Xkr4 protein coding
# MGI:98585 Trbv16 gene segment
# MGI:5530644 Mir6341 miRNA gene
# MGI:3643257 Gm7357 pseudogene
# MGI:2387894 Snord87 snoRNA gene
# MGI:4437728 D030040B21Rik lincRNA_gene
# MGI:5591663 Gm32504 lncRNA_gene
# MGI:4413856 n-TEttc1 tRNA_gene
#
# Gene with two ENSEMBL ids:
# MGI:105105 Rprl1  RNase P RNA gene, has 2 ensembl IDs, BIOTYPE CONFLICT
#
# Two genes with same NCBI id: 108168131
# MGI:4413970 n-TKttt23 tRNA gene
# MGI:5623715 Gm40830 lncRNA gene, BIOTYPE CONFLICT
#
source config.sh

MGIIDS="MGI:3528744\|MGI:98585\|MGI:5530644\|MGI:3643257\|MGI:2387894\|MGI:4437728\|MGI:5591663\|MGI:4413856\|MGI:105105\|MGI:4413970\|MGI:5623715"
echo "##gff-version 3"
${GREP} "${MGIIDS}" 


