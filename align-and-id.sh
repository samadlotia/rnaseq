GTF_PATH="/work/Common/Data/Annotation/human/2011_Archive/Homo_sapiens.GRCh37.56.chr.gtf"
ALIGN_BAM_PATH="tophat_out/accepted_hits.bam"
ALIGN_NAME=$(echo ${ALIGN_BAM_PATH%.*})
ALIGN_SAM_PATH="${ALIGN_NAME}.sam"
ALIGN_ANNOTATED_SAM_PATH="${ALIGN_NAME}-annotated.sam"
FASTQC_DIR_NAME="accepted_hits_fastqc"

function call_tophat() {
	nice tophat --no-novel-juncs --min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search --GTF="${GTF_PATH}" --num-threads 24 hg19 ${1}
}

function call_htseq() {
	samtools view "${ALIGN_BAM_PATH}" > "${ALIGN_SAM_PATH}"
	nice htseq-count --stranded=no -o "${ALIGN_ANNOTATED_SAM_PATH}" "${ALIGN_SAM_PATH}" "${GTF_PATH}" &> htseq-count-out
}

function call_fastqc() {
	echo /work/Apps/Bio/FastQC_RNASeq_QC_Java/fastqc_0.10.1/fastqc -o . "${ALIGN_BAM_PATH}"
}

for fastqname in `ls *.fastq`; do
	name=$(echo ${fastqname%.*})
	[ ! -d "$name" ] && mkdir "$name"
	cd "$name"
	echo
	echo "======================================"
	echo "$name"
	echo "======================================"
	[ ! -f "${ALIGN_BAM_PATH}" ] && call_tophat "../$fastqname"
	[ ! -d "${FASTQC_DIR_NAME}" ] && call_fastqc
	[ ! -f "${ALIGN_ANNOTATED_SAM_PATH}" ] && call_htseq
	cd ..
done
