FROM alpine/git AS clone
RUN git clone https://github.com/ANHIG/IMGTHLA.git /IMGTHLA


FROM intelliseqngs/hisat-genotype:test

COPY --from=clone /IMGTHLA/fasta/ /hisatgenotype/indices/hisatgenotype_db/fasta/
COPY --from=clone /IMGTHLA/msf/ /hisatgenotype/indices/hisatgenotype_db/msf/
COPY --from=clone /IMGTHLA/hla.dat /hisatgenotype/indices/hisatgenotype_db/
