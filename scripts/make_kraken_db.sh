# cd /lustre/scratch124/casm/team113/secure-lustre/resources/FUR/kraken2
# DBNAME=fur_kraken
# kraken2-build --download-taxonomy --db $DBNAME
# kraken2-build --download-library bacteria --db $DBNAME
# kraken2-build --download-library viral --db $DBNAME
# kraken2-build --download-library arachea --db $DBNAME
# kraken2-build --download-library human --db $DBNAME

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/350/175/GCF_018350175.1_F.catus_Fca126_mat1.0/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz

# # gunzip GCF_000002285.3_CanFam3.1_genomic.fna.gz
# # gunzip GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz


# bsub -q normal -M 16000 \
# -R'select[mem>16000] rusage[mem=16000] span[hosts=1]' \
# -n 6 \
# -o '/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/add_dog.o' \
# -e '/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/add_dog.e' \
# 'kraken2-build --add-to-library /lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_plus_catdog_july2024/GCF_000002285.3_CanFam3.1_genomic.fna --db kraken2_plus_catdog_july2024'


# bsub -q normal -M 16000 \
# -R'select[mem>16000] rusage[mem=16000] span[hosts=1]' \
# -n 6 \
# -o '/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/add_cat.o' \
# -e '/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/add_cat.e' \
# 'kraken2-build --add-to-library /lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_plus_catdog_july2024/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna --db kraken2_plus_catdog_july2024'


cat pathogen_identification/scripts/pathogens_plus_db.txt | xargs -I {} sh -c './datasets download genome taxon "{}" --assembly-source RefSeq --assembly-version latest --reference --annotated --filename "{}".zip'