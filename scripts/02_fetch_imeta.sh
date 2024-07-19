bsub -q long -n 6 -M16000 -R"select[mem>16000] rusage[mem=16000] span[hosts=1]" -Is /bin/bash
module load IRODS/1.0
iinit
 
#To call the script to build the manifests
unset R_LIBS
export R_LIBS=/software/team113/dermatlas/R/R-4.2.2/lib/R/library


PATH_TO_SEARCH="/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis"
for dir in $PATH_TO_SEARCH/*/; do
    # Extract the STUDY number (first number in the directory name)
    STUDY=$(basename "$dir" | cut -d'_' -f1)
    echo $STUDY
    # Define PROJECTDIR as the full path of the current directory
    PROJECTDIR=$(realpath "$dir")
    echo $PROJECTDIR
    /software/team113/dermatlas/R/R-4.2.2/bin/Rscript \
    scripts/dermatlas-starfusion/scripts/Build_manifest_from_irods_cram_information.R \
    --seqscape_proj_id ${STUDY} \
    --outdir $PROJECTDIR/metadata;
done
