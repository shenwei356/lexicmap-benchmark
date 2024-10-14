## AllTheBacteria v0.2 High-quality genomes

## Data

    #.genomes: 1,858,610
    #.bases: 7,493,622,021,123
    #.gzip_size: 3.1 TB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -b 25000 -O atb_hq.lmi --force" > atb_hq.lmi.log 2>&1

    elapsed time: 2.0days 0h:07m:03s
    peak rss: 88.56 GB

    atb_hq.lmi: 3.88 TB
       2.11 TB      seeds
       1.77 TB      genomes
      39.22 MB      genomes.map.bin
     312.53 KB      masks.bin
      332.00 B      info.toml


Searching

    db=atb_hq.lmi

    ls b.*.fasta | tac |  while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # --------------------------------------------------------------------------

    # hits
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk uniq -t -f query,sgenome {} | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       21288000
    b.gene_E_coli_16S.fasta.lexicmap.tsv            1857761
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         27963
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv     468821

    # -----------------------------------

    # hits with qcovGnm > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovGnm >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       12274740
    b.gene_E_coli_16S.fasta.lexicmap.tsv            1856212
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         27954
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv      11047


    # -----------------------------------

    # hits with qcovHSP > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovHSP >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       12148642
    b.gene_E_coli_16S.fasta.lexicmap.tsv            1740000
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         27953
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv       3618


    # --------------------------------------------------------------------------

    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.amr.fasta.lexicmap.tsv.log
    elapsed time: 2h:18m:55s
    peak rss: 49.92 GB

    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 9m:36s
    peak rss: 14.87 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 30.638s
    peak rss: 3.41 GB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 15m:55s
    peak rss: 15.71 GB




Search with a serial numbers of AMR genomes.

    mkdir -p amrs; cd amrs

    seqkit range -r     1:1 ../b.amr.fasta -o b.amr_0001.fasta
    seqkit range -r     2:3 ../b.amr.fasta -o b.amr_0002.fasta
    seqkit range -r    6:10 ../b.amr.fasta -o b.amr_0005.fasta
    seqkit range -r   11:20 ../b.amr.fasta -o b.amr_0010.fasta
    seqkit range -r   21:40 ../b.amr.fasta -o b.amr_0020.fasta
    seqkit range -r  51:100 ../b.amr.fasta -o b.amr_0050.fasta
    seqkit range -r 101:200 ../b.amr.fasta -o b.amr_0100.fasta
    seqkit range -r 201:400 ../b.amr.fasta -o b.amr_0200.fasta
    seqkit range -r 401:900 ../b.amr.fasta -o b.amr_0500.fasta


    seqkit stats b.amr_*.fasta

    file              format  type  num_seqs  sum_len  min_len  avg_len  max_len
    b.amr_0001.fasta  FASTA   DNA          1    1,010    1,010    1,010    1,010
    b.amr_0002.fasta  FASTA   DNA          2    2,675    1,289  1,337.5    1,386
    b.amr_0005.fasta  FASTA   DNA          5    4,874      861    974.8    1,227
    b.amr_0010.fasta  FASTA   DNA         10   10,115      720  1,011.5    1,346
    b.amr_0020.fasta  FASTA   DNA         20   21,431      645  1,071.6    1,829
    b.amr_0050.fasta  FASTA   DNA         50   49,732      398    994.6    1,940
    b.amr_0100.fasta  FASTA   DNA        100  108,418      437  1,084.2    2,107
    b.amr_0200.fasta  FASTA   DNA        200  223,403      378    1,117    3,335
    b.amr_0500.fasta  FASTA   DNA        500  541,243      374  1,082.5    3,275


    db=../atb_hq.lmi
    ls b.amr_*.fasta | while read q; do \
        echo $q; memusg -H -t -s "lexicmap search -d $db $q -n 0 -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # extract result
    ls b.amr_*.lexicmap.tsv.log \
        | rush -k 't=$(tail -n 4 {} | grep elapsed | sed -E "s/.+: //"); \
                    m=$(tail -n 4 {} | grep peak | sed -E "s/.+: //"); \
                    echo -e "{}\t$t\t$m";' \
        | csvtk replace -Ht -p '.+amr_0+(\d+)\..+' -r '$1' \
        | csvtk add-header -t -n n,time,mem \
        > amrs.tsv

    Rscript amrs.R

## BLASTN

Indexing

    memusg -t -s "rush -i files.txt 'zcat {}' -j 12 -k --eta | makeblastdb -dbtype nucl -title blastdb -in - -out blastdb/blastdb" > blastdb.log 2>&1

    time: 14h:03m:21s
     mem: 2.88 GB
    size: 1.76 TB

Searching

    blastdb=blastdb/blastdb
    ls b.*fasta | while read q; do \
        echo $q; \
        memusg -t -s "blastn -num_threads 48 -max_target_seqs 10000000 -db $blastdb -query $q \
              -outfmt '7 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
            | gzip -c > $q.blastn.tsv.gz" > $q.blastn.tsv.gz.log 2>&1; \
    done


## Phylign


Seaching

    # --------------------------------------
    # b.amr.fasta

    elapsed time: 2h:36m:08s
    peak rss: 85.92 GB

    queries                       1033
    cumul_length_bps              1108417
    matched_queries               1031
    aligned_queries               666
    aligned_segments              15640595
    distinct_genome_query_pairs   11741842
    target_genomes                1135215
    target_batches                587
    nonalignments                 757770

    # cluster
    real    133m49.206s
    user    37m25.037s
    sys     7m29.366s

    queries                       1033
    cumul_length_bps              1108417
    matched_queries               1031
    aligned_queries               645
    aligned_segments              9862303
    distinct_genome_query_pairs   9278046
    target_genomes                1133995
    target_batches                587
    nonalignments                 3221566


    # --------------------------------------
    # b.gene_E_coli_16S.fasta

    elapsed time: 2h:10m:33s
    peak rss: 77.02 GB


    queries                       1
    cumul_length_bps              1542
    matched_queries               1
    aligned_queries               1
    aligned_segments              1835351
    distinct_genome_query_pairs   1017765
    target_genomes                1017765
    target_batches                418
    nonalignments                 15184

    # cluster
    real    86m41.041s
    user    18m34.817s
    sys     6m36.640s

    # --------------------------------------
    # b.gene_E_faecalis_SecY.fasta

    elapsed time: 30m:48s
    peak rss: 77.62 GB

    queries                       1
    cumul_length_bps              1299
    matched_queries               1
    aligned_queries               1
    aligned_segments              7938
    distinct_genome_query_pairs   7936
    target_genomes                7936
    target_batches                3
    nonalignments                 1957

    # cluster
    real    28m33.340s
    user    6m48.370s
    sys     2m25.836s

    # --------------------------------------
    # b.plasmid_pCUVET18-1784.4.fasta

    elapsed time: 47m:33s
    peak rss: 82.56 GB

    queries                       1
    cumul_length_bps              52830
    matched_queries               1
    aligned_queries               1
    aligned_segments              348710
    distinct_genome_query_pairs   46822
    target_genomes                46822
    target_batches                275
    nonalignments                 0

    # cluster

    real    39m33.540s
    user    7m51.464s
    sys     2m47.193s
