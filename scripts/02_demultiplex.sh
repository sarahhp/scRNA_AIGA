# ══════════════════════════════════════════════════════════════════════════════
# 02_demultiplex.sh — cellsnp-lite + vireo (run on Gadi as PBS job)
# ══════════════════════════════════════════════════════════════════════════════

#!/bin/bash
#PBS -P pq84
#PBS -q normal
#PBS -l ncpus=28
#PBS -l mem=128GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/pq84
#PBS -N demultiplex

cd /g/data/pq84/single_cell/10X/scRNA_auto/all_data/GIMR_GGP_CIRCA_20250513_5P

CELLSNP=/g/data/pq84/mano/programs/cellsnp-lite/bin/cellsnp-lite
VCF_ORIGINAL=/path/to/your/donors.vcf.gz

# ── Step 1: Add chr prefix to VCF (hg19 → hg38 chromosome naming) ──────────
cat > chr_rename.txt << 'EOF'
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
Y chrY
EOF

VCF_CHR="${VCF_ORIGINAL%.vcf.gz}_chr.vcf.gz"

if [ ! -f "$VCF_CHR" ]; then
  bcftools annotate --rename-chrs chr_rename.txt \
    $VCF_ORIGINAL -Oz -o $VCF_CHR
  tabix -p vcf $VCF_CHR
fi

# ── Step 2: cellsnp-lite + vireo per pool ────────────────────────────────────
for POOL in A B; do

  BAM="Pool_${POOL}/cellranger_count/Pool_${POOL}/outs/possorted_genome_bam.bam"
  BARCODES="Pool_${POOL}/cellranger_count/Pool_${POOL}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  CELLSNP_OUT="cellSNP_Pool_${POOL}_chr_out"
  VIREO_OUT="vireo_Pool_${POOL}_out"

  # cellsnp-lite
  if [ ! -d "$CELLSNP_OUT" ]; then
    echo "Running cellsnp-lite for Pool_${POOL}..."
    $CELLSNP -s $BAM \
      --barcode $BARCODES \
      --cellTAG CB \
      --UMItag UB \
      --regions $VCF_CHR \
      --gzip \
      --outDir $CELLSNP_OUT \
      -p 28
  fi

  # vireo
  if [ ! -d "$VIREO_OUT" ]; then
    echo "Running vireo for Pool_${POOL}..."
    vireo \
      --cellData $CELLSNP_OUT \
      --donorFile $VCF_CHR \
      --outDir $VIREO_OUT \
      --genoTag GT
  fi

  echo "Pool_${POOL} demultiplexing complete"
  cat ${VIREO_OUT}/summary.tsv

done
