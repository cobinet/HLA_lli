FROM biocontainers/bowtie2:v2.4.1_cv1

USER root

WORKDIR /usr/local/
COPY hlahd.1.7.0.tar.gz /usr/local/
RUN tar -zxvf hlahd.1.7.0.tar.gz
RUN rm hlahd.1.7.0.tar.gz

WORKDIR /usr/local/hlahd.1.7.0
RUN sh install.sh
ENV PATH=$PATH:/usr/local/hlahd.1.7.0/bin/
# Select the newest dictionary
RUN mv HLA_gene.split.txt HLA_gene.split.default.txt
RUN mv HLA_gene.split.3.50.0.txt HLA_gene.split.txt
# Update the dictionary
RUN sh update.dictionary.sh

WORKDIR /work
