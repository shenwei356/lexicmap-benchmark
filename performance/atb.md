## AllTheBacteria v0.2 High-quality genomes

## Data

    #.genomes: 1,858,610
    #.bases: 7,493,622,021,123
    #.gzip_size: 3.1 TB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -b 50000 -O atb_hq.lmi --force" > atb_hq.lmi.log 2>&1

    elapsed time: 11h:51m:03s
    peak rss: 42.62 GB

    atb_hq.lmi: 2.32 TB
       1.77 TB      genomes
     569.44 GB      seeds
      39.22 MB      genomes.map.bin
     312.53 KB      masks.bin
      270.00 B      info.toml

    # 500bp
    elapsed time: 18h:21m:20s
    peak rss: 42.81 GB

    atb_hq.lmi: 2.37 TB
       1.77 TB      genomes
     614.58 GB      seeds
      39.22 MB      genomes.map.bin
     312.53 KB      masks.bin
      332.00 B      info.toml

Searching

    db=atb_hq.lmi

    ls b.*fasta | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -n 0 -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done


    # hits
    ls b.*.tsv | rush -k 'echo -ne "{}\t" ; csvtk head -n 1 -t {} | csvtk cut -t -f hits -U;'

    b.gene_E_coli_16S.fasta.lexicmap.tsv    1854616
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       11205
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    427544


    # hits with qcovGnm > 50
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; \
        csvtk filter2 -t -f "\$qcovGnm >= 50" {} | csvtk uniq -t -f sgenome | csvtk nrow -t '

    b.gene_E_coli_16S.fasta.lexicmap.tsv    1849742
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       11199
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    11082


    # hits with qcovHSP > 50
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; \
        csvtk filter2 -t -f "\$qcovHSP >= 50" {} | csvtk uniq -t -f sgenome | csvtk nrow -t '

    b.gene_E_coli_16S.fasta.lexicmap.tsv    1735763
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       11199
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    3620


    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 6m:43s
    peak rss: 14.7 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 3.226s
    peak rss: 1.83 GB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 3m:42s
    peak rss: 19.24 GB



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
    # b.gene_E_faecalis_SecY.fasta

    queries                       1
    cumul_length_bps              1299
    matched_queries               1
    aligned_queries               1
    aligned_segments              7940
    distinct_genome_query_pairs   7937
    target_genomes                7937
    target_batches                4
    nonalignments                 1956


    # --------------------------------------
    # b.gene_E_coli_16S.fasta
    queries                       1
    cumul_length_bps              1542
    matched_queries               1
    aligned_queries               1
    aligned_segments              2081870
    distinct_genome_query_pairs   1032948
    target_genomes                1032948
    target_batches                498
    nonalignments                 1


    # --------------------------------------
    # b.plasmid_pCUVET18-1784.4.fasta

    queries                       1
    cumul_length_bps              52830
    matched_queries               1
    aligned_queries               1
    aligned_segments              348710
    distinct_genome_query_pairs   46822
    target_genomes                46822
    target_batches                275
    nonalignments                 0

