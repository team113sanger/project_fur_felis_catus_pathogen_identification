
#This script takes the cram manifest and generates the SH file with the farm jobs to import and transform to fastq all of the cram files from iRODs

PATH_TO_SEARCH="/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis"
for manifest in $PATH_TO_SEARCH/*/metadata/*INFO_from_iRODS.txt; do
    # Extract the STUDY number (first number in the directory name)
    STUDY=$(basename "$manifest" | cut -d'_' -f1)
    # Define PROJECTDIR as the full path of the current directory
    PROJECTDIR=$(dirname $( dirname "$manifest"))
    echo $PROJECTDIR
    export TMPDIR=${PROJECTDIR}/tmp
    /software/team113/dermatlas/R/R-4.2.2/bin/Rscript \
    scripts/dermatlas-starfusion/scripts/cramtofastq_from_iRODs_based_cram_manifest.R \
    --manifest $manifest \
    --projectdir ${PROJECTDIR} \
    --studyID ${STUDY} --mem 16000
done

module load samtools-1.14/python-3.12.0 

/bin/sh /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6555_2711/scripts/6555_cramtofastq_from_iRODs_jobs.sh
/bin/sh /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6711_2820/scripts/6711_cramtofastq_from_iRODs_jobs.sh
# /bin/sh /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6712_2822/scripts/6712_cramtofastq_from_iRODs_jobs.sh
/bin/sh /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6713_2821/scripts/6713_cramtofastq_from_iRODs_jobs.sh
/bin/sh /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6841_2964/scripts/6841_cramtofastq_from_iRODs_jobs.sh
/bin/sh /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6864_2965/scripts/6864_cramtofastq_from_iRODs_jobs.sh