FROM biocontainers/bowtie2:v2.4.1_cv1
COPY hlahd.1.6.1.tar.gz .
RUN tar -zxvf hlahd.1.6.1.tar.gz
RUN cd ~/.local/share/hlahd.1.6.1
RUN export PATH=$PATH:target_path/hlahd.x.y.z/bin
# Select the newest dictionary
RUN mv HLA_gene.split.txt HLA_gene.split.default.txt
RUN mv HLA_gene.split.3.50.0.txt HLA_gene.split.txt
# Update the dictionary
RUN sh update.dictionary.sh
# You should check in dictionary/ that all the genes are there
