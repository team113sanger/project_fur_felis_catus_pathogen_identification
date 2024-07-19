screen -S 3033_cp
bsub -q long -n 6 -M16000 -R"select[mem>16000] rusage[mem=16000] span[hosts=1]" -Is /bin/bash
 

STUDY=6982
CANAPPS_ID=
PROJECTDIR="/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/${STUDY}_${CANAPPS_ID}"
#Set Sequencescape project ID
mkdir $PROJECTDIR
mkdir -p $PROJECTDIR/metadata
 
module load IRODS/1.0
iinit
 
#To call the script to build the manifests
unset R_LIBS
export R_LIBS=/software/team113/dermatlas/R/R-4.2.2/lib/R/library
 
/software/team113/dermatlas/R/R-4.2.2/bin/Rscript \
./scripts/Build_manifest_from_irods_cram_information.R \
--seqscape_proj_id ${STUDY} \
--outdir $PROJECTDIR/metadata

 
mkdir ${PROJECTDIR}/tmp
export TMPDIR=${PROJECTDIR}/tmp
#This script takes the cram manifest and generates the SH file with the farm jobs to import and transform to fastq all of the cram files from iRODs
 

 PROJECTDIR=/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/KRAKEN/6982_Feline_lymphoma

/software/team113/dermatlas/R/R-4.2.2/bin/Rscript \
./scripts/cramtofastq_from_iRODs_based_cram_manifest.R \
--manifest ${STUDY}_cram_manifest_INFO_from_iRODS.txt \
--projectdir ${PROJECTDIR} \
--studyID ${STUDY} --mem 16000

module load IRODS
module load samtools-1.14/python-3.12.0

/bin/sh scripts/