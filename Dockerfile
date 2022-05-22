# Pull miniconda3 image
FROM continuumio/miniconda3:4.7.12

# Create /genomika-snnb directory and set it as work directory
WORKDIR /biotest

# Install tabix
RUN conda install -c bioconda tabix

# Install fastqc
RUN conda install conda-forge::openjdk && \
    conda install -c bioconda fastqc && \
    conda install -c bioconda fastp=0.20.1


# Install bwa
RUN conda install -c bioconda bwa

# Install samtools
RUN conda install -c bioconda samtools=1.9 --force-reinstall

# Install Freebayes
RUN conda install -c bioconda freebayes

# Install vcftools
RUN conda install -c bioconda vcftools

# Install vcftools
RUN conda install -c bioconda bcftools

# Install snpEff
RUN conda install -c bioconda snpeff=5.0

#INSTALL snpEff hg19 database
RUN snpEff download hg19

CMD ["/bin/bash"]
