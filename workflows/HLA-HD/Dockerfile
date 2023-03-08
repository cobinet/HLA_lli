FROM biocontainers/bowtie2:v2.4.1_cv1

USER root

WORKDIR /usr/local/share/
COPY hlahd.1.7.0.tar.gz /usr/local/share/
RUN tar -zxvf hlahd.1.7.0.tar.gz
WORKDIR /usr/local/share/hlahd.1.7.0
RUN export PATH=$PATH:target_path/hlahd.1.7.0/bin/
# Select the newest dictionary
RUN mv HLA_gene.split.txt HLA_gene.split.default.txt
RUN mv HLA_gene.split.3.50.0.txt HLA_gene.split.txt
# Update the dictionary
RUN sh update.dictionary.sh

WORKDIR /work

USER biodocker