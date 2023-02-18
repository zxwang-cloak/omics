# omics
The pipeline for upstream analysis of WGS/WGBS/RNA-seq

## WGS mutation calling
### Data pre-processing
_Step 1: Alignment_
```
bwa mem -t 20 -M /path/reference/hg38/hg38.fa \ 
    /path/Rawdata/SAMPLE_ ID/SAMPLE_ID_1.fq.gz \ 
    /path/Rawdata/SAMPLE_ ID/SAMPLE_ID_2.fq.gz \
    > /path/out/bwa_mem/SAMPLE_ ID.bwa.mem.sam
```
_Step 2: Sort_
```
samtools sort -@ 20 -m 5G -O bam -T SAMPLE_ID \
    -o /path/out/bwa_mem/SAMPLE_ID_sorted.bam \
    /path/out/bwa_mem/SAMPLE_ID.bwa.mem.sam
```
_Step 3: Mark duplication_
```
gatk --java-options "-Xmx50G" MarkDuplicates \
    --TMP_DIR /path/tmp \
    -I /path/out/bwa_mem/SAMPLE_ID _sorted.bam \
    -M /path/out/markdup/SAMPLE_ID _metrics.txt \
    -O /path/out/markdup/SAMPLE_ID _sorted_markdup.bam
```
_Step 4: Add label name (optional, not used if grouped)_
```
picard AddOrReplaceReadGroups \
    -I /path/out/markdup/SAMPLE_ID_sorted_markdup.bam \
    -O /path/out/add_RG/SAMPLE_ID_sorted_markdup.R.bam \
    --RGID SAMPLE_ID \
    --RGLB lib1 \
    --RGPL illumina \
    --RGPU unit1 \
    --RGSM SAMPLE_ID \
    --VALIDATION_STRINGENCY LENIENT \
    -Xms10g -Xmx50g -XX:ParallelGCThreads=10
```
_Step 5: Builds a model of covariation_
```
gatk --java-options “-Xmx50G” BaseRecalibrator \
    -I /path/out/add_RG/SAMPLE_ID_sorted_markdup.R.bam \
    -R /path/reference/hg38/hg38.fa \
    --known-sites /path/reference/hg38/gatk/1000G_phase1.snps.high_confidence.hg38.vcf \
    --known-sites /path/reference/hg38/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    --known-sites /path/reference/hg38/gatk/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O /path/out/pre-BQSR/SAMPLE_ID.recal_data.table
```
_Step 6: Base quality score recalibration (BQSR)_
```
gatk --java-options "-Xmx50G" ApplyBQSR \
    --bqsr-recal-file /path/out/pre-BQSR/SAMPLE_ID.recal_data.table \
    -I /path/out/add_RG/SAMPLE_ID_sorted_markdup.R.bam \
    -R /path/reference/hg38/hg38.fa \
    -O /path/out/BQSR/SAMPLE_ID_sorted_markdup_R_BQSR.bam
```

### Somatic SNV/indel    
**Generate the PON (panel of normal) file**

_Step 1: Run Mutect2 in tumor-only mode for each normal sample_
```
gatk Mutect2 \
    -R /path/reference/hg38/hg38.fa \
    -I /path/out/BQSR/NORMAL_sorted_markdup_R_BQSR.bam \
    -max-mnp-distance 0 \
    --independent-mates \
    -O /path/out/toPON/NORMAL.vcf
```
_Step 2: Create a GenomicsDB from the normal Mutect2 calls_
```
gatk GenomicsDBImport \
    -R /path/out/BQSR/NORMAL_sorted_markdup_R_BQSR.bam \
    -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y \
    --genomicsdb-workspace-path /path/out/toPON \
    --tmp-dir /path/tmp \
    --sample-name-map /path/out/toPON/sample-name.list
```
_Step 3: Combine the normal calls using CreateSomaticPanelOfNormals_
```
gatk CreateSomaticPanelOfNormals \
    -R /path/reference/hg38/hg38.fa \
    --germline-resource /path/reference/hg38/gatk/af-only-gnomad.hg38.vcf.gz \
    -V gendb://./out/PON_db \
    -O /path/out/pon.vcf.gz
```
**somatic SNV/indel calling & filter**

_Step 1: Run Mutect2 in tumor-normal mode_
```
gatk Mutect2 \
    -R /path/reference/hg38/hg38.fa \
    -I /path/out/BQSR/NORMAL_sorted_markdup_R_BQSR.bam \
    -normal NORMAL \
    -I /path/out/BQSR/TUMOR_sorted_markdup_R_BQSR.bam \
    -tumor TUMOR \
    --germline-resource /path/reference/hg38/gatk/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals /path/out/pon.vcf.gz \
    --f1r2-tar-gz /path/out/MuTect2/TUMOR_f1r2.tar.gz \
    -O /path/out/MuTect2/TUMOR_mutect2.vcf.gz
```
_Step 2: Get pileup summaries_
```
gatk --java-options "-Xmx50G" GetPileupSummaries \
    -R /path/reference/hg38/hg38.fa \
    -L /path/reference/hg38/gatk/wgs_calling_regions.hg38.interval_list \
    -V /path/reference/hg38/gatk/af-only-gnomad.hg38.SNP_biallelic.vcf.gz \
    -I /path/out/BQSR/SAMPLE_ID_sorted_markdup_R_BQSR.bam \
    -O /path/out/Pileup/SAMPLE_ID-GetPileupSummaries.table
```
_Step 3: Calculate contamination_
```
gatk --java-options "-Xmx50G" CalculateContamination \
    -I /path/out/Pileup/TUMOR-GetPileupSummaries.table \
    -matched /path/out/Pileup/NORMAL-GetPileupSummaries.table \
    -O /path/out/Contamination/TUMOR-contamination.table
```
_Step 4: Learn the parameters of a model for orientation bias_
```
gatk --java-options "-Xmx50G" LearnReadOrientationModel \
    -I /path/out/MuTect2/TUMOR_f1r2.tar.gz \
    -O /path/out/read-orientation/TUMOR_read-orientation-model.tar.gz
```
_Step 5: Filter_
```
gatk --java-options "-Xmx50G" FilterMutectCalls \
    -R /path/reference/hg38/hg38.fa \
    -V /path/out/MuTect2/TUMOR_mutect2.vcf.gz \
    --ob-priors /path/out/read-orientation/ TUMOR_read-orientation-model.tar.gz \
    --contamination-table /path/out/Contamination/TUMOR-contamination.table \
    -O /path/out/Filter/TUMOR_filtered.vcf.gz
```
### Somatic SV
**calling using Manta**:
https://github.com/Illumina/manta

_Step 1: Generate config file of Manta_
```
configManta.py \
    --normalBam /path/out/BQSR/NORMAL_sorted_markdup_R_BQSR.bam \
    --tumorBam /path/out/BQSR/TUMOR_sorted_markdup_R_BQSR.bam \
    --referenceFasta /path/reference/hg38/hg38.fa \
    --runDir /path/out/Manta/candidateSmallIndels/TUMOR
```
_Step 2: Run the workflow script_
```
cd /path/out/Manta/candidateSmallIndels/TUMOR && python runWorkflow.py
```

### Somatic CNV




## WGBS CpG calling




## RNA-seq gene/transcript expression calling







