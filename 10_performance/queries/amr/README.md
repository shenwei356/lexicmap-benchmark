## Data

Download from: https://www.ncbi.nlm.nih.gov/bioproject/313047


    seqkit stats sequence.fasta  -a
    file            format  type  num_seqs    sum_len  min_len  avg_len  max_len   Q1   Q2     Q3  sum_gap    N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)
    sequence.fasta  FASTA   DNA      6,915  7,100,828      374  1,026.9    3,353  825  961  1,152        0  1,052      367       0       0        0  50.22


Count genes

    cat sequence.fasta \
        | seqkit seq -n \
        | csvtk mutate -Ht -p '(\S+) gene' \
        | csvtk freq -Ht -f 2 -nr \
        | tee genes.tsv \
        | head -n 5 \
        | csvtk pretty -Ht
    blaOXA     1197
    blaPDC     573
    blaADC     322
    blaCTX-M   254
    blaKPC     208

Reformat

    cat sequence.fasta \
        | seqkit fx2tab \
        | csvtk mutate -Ht -p '(\S+) gene' \
        | csvtk mutate -Ht -p '^\S+\s+(.+)' \
        | csvtk replace -Ht -p '\s+.+' \
        | awk -F'\t' '{print ">"$4"__"$1" "$5"\n"$2}' \
        > amr.fasta

Deduplication, keep one record for each gene.

    seqkit rmdup --id-regexp '^(.+?)__' amr.fasta -o b.amr.fasta

    seqkit stats -a b.amr.fasta
    file            format  type  num_seqs    sum_len  min_len  avg_len  max_len   Q1     Q2     Q3  sum_gap    N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)
    amr.uniq.fasta  FASTA   DNA      1,033  1,108,417      374    1,073    3,335  809  1,001  1,214        0  1,093      222       0       0        0  48.18

