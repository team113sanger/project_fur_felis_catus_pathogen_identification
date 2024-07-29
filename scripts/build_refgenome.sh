# docker pull quay.io/biocontainers/bowtie2:2.5.4--he20e202_2
# docker-compose up 
bowtie2-build /lustre/scratch124/resources/FUR/GCF_000002285.3_CanFam3.1_genomic.fna CanFam3.1 -t 8 
bowtie2-build /lustre/scratch124/resources/FUR/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna FelisCatus -t 8 