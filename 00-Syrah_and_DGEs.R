######################################################
# SYRAH AND DGEs
#
# GOAL: 
######################################################

# ================= IMPORTS ==========================

# ================= PARAMS ===========================

# ================= METHODS ==========================

# ================= INIT DATA ========================

# ================= ANALYSIS =========================

# ================= SAVE DATA ========================
batches <- c('L46130','L46131','L46144','L46145')

pucks <- c('L46130'='/n/analysis/slideseq/210113_14/puck.txt',
           'L46131'='/n/analysis/slideseq/210113_13/puck.txt',
           'L46144'='/n/analysis/slideseq/Puck_210119_29/puck.txt',
           'L46145'='/n/analysis/slideseq/Puck_210119_27/puck.txt')

newDir <- '/home/cb2350/slide-seq/smed/data/'
tempDir <- '/home/cb2350/slide-seq/TEMP/'

codeDir <- '/home/cb2350/slide-seq/smed/code/'
coreBamDir <- c('/n/core/Bioinformatics/tools/slideseq-tools/runs/HGFG3BGXJ/libraries/2021-5-22_')

for (batch in batches[2:length(batches)]) {
  puck <- pucks[[batch]]
  bamDir <- paste0(coreBamDir,batch,'/')
  cat(batch,puck,bamDir,'\n')
  
  # Create bead barcode mapping
  system(paste0('nohup Rscript ',codeDir,'create_bead_barcode_mapping.R ',puck,' ',batch))
  
  # get BAM names
  bams <- list.files(bamDir)[which(endsWith(list.files(bamDir),'.bam'))]
  r1Bams <- bams[which(endsWith(bams,'unmapped.bam'))]
  r2Bam <- bams[which(!endsWith(bams,'unmapped.bam'))]
  
  # R2 = MAPPED, NEED ALIGNMENT DATA
  system(paste0('samtools view -h -F 4 -O bam -o ',newDir,batch,'_r2_filtered.bam ',bamDir,r2Bam))
  system(paste0('samtools sort -n -T /home/cb2350/slide-seq/TEMP/ -O bam -o ',newDir,batch,'_r2_filtered_sorted.bam ',newDir,batch,'_r2_filtered.bam'))
  system(paste0('samtools view -H ',newDir,batch,'_r2_filtered_sorted.bam | grep -v ^@PG > ',newDir,batch,'_r2_header_no_PG.txt'))
  system(paste0('samtools view ',newDir,batch,'_r2_filtered_sorted.bam > ',newDir,batch,'_r2_reads_only.sam'))
  system(paste0('cat ',newDir,batch,'_r2_reads_only.sam | cut -f1 > ',newDir,batch,'_r2_qnames.txt'))
  
  # merge r1 bams
  n1Reads <- sum(sapply(r1Bams,function(x){ as.integer(system(paste0('samtools view -c ',bamDir,x),intern=TRUE)) }))
  # R1 = UNMAPPED, NEED R1 42 NT SEQUENCE FOR BC/UMI RE-MAPPING
  system(paste0('samtools merge ',newDir,batch,'_r1_merged.bam ',paste(paste0(bamDir,r1Bams),collapse=' ')))
  
  # CHECK TO MAKE SURE WE DIDN'T LOSE READS WHEN MERGING, AS SOMETIMES SEEMS TO HAPPEN (WHY??)
  newN1Reads <- as.integer(system(paste0('samtools view -c ',newDir,batch,'_r1_merged.bam'),intern=TRUE))
  if (newN1Reads!=n1Reads) { paste0('merged ',batch,' read1 BAM should have ',n1Reads,' reads but only has ',newN1reads,'\n')}
  
  system(paste0('samtools view -f 77 -o ',newDir,batch,'_r1_merged_filtered.bam ',newDir,batch,'_r1_merged.bam'))
  system(paste0('samtools sort -n -o ',newDir,batch,'_r1_merged_filtered_sorted.sam ',newDir,batch,'_r1_merged_filtered.bam'))
  system(paste0('samtools view -N ',newDir,batch,'_r2_qnames.txt -o ',newDir,batch,'_r1_merged_filtered_sorted_qname_filtered.sam ',newDir,batch,'_r1_merged_filtered_sorted.sam'))
  system(paste0('cat ',newDir,batch,'_r1_merged_filtered_sorted_qname_filtered.sam | cut -f1,10 > ',newDir,batch,'_r1_qname_seq_only_r2_qname_filtered.txt'))
  
  
  # CORRECT BEAD (XC) AND MOLECULAR (XM) BARCODES
  system(paste0('nohup Rscript ',codeDir,'correct_XC_XM.R ',newDir,batch,' ',newDir,'puck_barcode_mappings/',batch,'_barcode_mapping.fa'))
  # filter by linker alignment quality
  system(paste0('cat ',newDir,batch,"_r2_XC_XM_corrected.txt | grep \'XX:Z:[012345]d[89]s\' > ",newDir,batch,'_r2_XC_XM_corrected_filtered.txt'))
  # re-create SAM --> BAM
  system(paste0('cat ',newDir,batch,'_r2_header_no_PG.txt ',newDir,batch,'_r2_XC_XM_corrected_filtered.txt > ',newDir,batch,'_r2_XC_XM_corrected.sam'))
  system(paste0('nohup samtools view -O bam -o ',newDir,batch,'_corrected.bam ',newDir,batch,'_r2_XC_XM_corrected.sam'))
  # get bead barcodes from puck
  #system(paste0('cat ',puck,' | cut -f1 > ',newDir,batch,'_beads.txt'))
  # create DGE
  system(paste0('/n/core/Bioinformatics/tools/Drop-seq_tools-2.3.0/DigitalExpression',
                ' -m 120g',
                ' I=',newDir,batch,'_corrected.bam',
                ' O=',newDir,batch,'.dge.txt.gz',
                ' CELL_BARCODE_TAG=XC',
                ' MOLECULAR_BARCODE_TAG=XM',
                ' GENE_NAME_TAG=gn',
                ' MAX_RECORDS_IN_RAM=500000',
                ' USE_STRAND_INFO=true',
                ' RARE_UMI_FILTER_THRESHOLD=0.0',
                ' SUMMARY=',newDir,batch,'_expression_summary.txt',
                ' MIN_NUM_TRANSCRIPTS_PER_CELL=10',
                ' TMP_DIR=/home/cb2350/slide-seq/TEMP/',
                ' VALIDATION_STRINGENCY=SILENT',
                ' OUTPUT_READS_INSTEAD=false',
                ' COMPRESSION_LEVEL=5',
                ' CREATE_INDEX=false',
                ' CREATE_MD5_FILE=false',
                ' GA4GH_CLIENT_SECRETS=client_secrets.json',
                ' USE_JDK_DEFLATER=false',
                ' USE_JDK_INFLATER=false'))
}