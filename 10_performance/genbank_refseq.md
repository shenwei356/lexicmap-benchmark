## Genbank+RefSeq

## Data

    #.genomes: 2,340,672
    #.bases: 9,192,667,695,899
    #.gzip_size: 3.5 TB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -b 25000 -O genbank_refseq.lmi -G files.txt.big.txt" > genbank_refseq.lmi.log 2>&1

    elapsed time: 2.0days 6h:32m:42s
    peak rss: 178.23 GB

    genbank_refseq.lmi: 4.94 TB
       2.77 TB      seeds
       2.17 TB      genomes
      55.81 MB      genomes.map.bin
     312.53 KB      masks.bin
      332.00 B      info.toml

Searching

    db=genbank_refseq.lmi

    ls b.*fasta | tac | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # --------------------------------------------------------------------------

    # hits
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk uniq -t -f query,sgenome {} | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       25702419
    b.gene_E_coli_16S.fasta.lexicmap.tsv            1949496
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         37164
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv     544619

    # -----------------------------------

    # hits with qcovGnm > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovGnm >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       14865098
    b.gene_E_coli_16S.fasta.lexicmap.tsv            1666548
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         37085
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv      18039

    # -----------------------------------

    # hits with qcovHSP > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovHSP >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       14692624
    b.gene_E_coli_16S.fasta.lexicmap.tsv            1381974
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         37082
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv       6563

    # --------------------------------------------------------------------------

    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.amr.fasta.lexicmap.tsv.log
    elapsed time: 3h:07m:40s
    peak rss: 55.43 GB

    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 10m:41s
    peak rss: 14.06 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 36.162s
    peak rss: 4.04 GB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 19m:20s
    peak rss: 19.28 GB


## BLASTN

Indexing

    memusg -t -s "rush -i files.txt 'zcat {}' -j 12 -k --eta | makeblastdb -dbtype nucl -title blastdb -in - -out blastdb/blastdb" > blastdb.log 2>&1

    time: 14 h 4 m
     mem: 4.29 GB
    size: 2.15 TB

Searching

    blastdb=blastdb/genbank_refseq
    ls b.*fasta | while read q; do \
        echo $q; \
        memusg -t -s "blastn -num_threads 48 -max_target_seqs 10000000 -db $blastdb -query $q \
              -outfmt '7 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
            | gzip -c > $q.blastn.tsv.gz" > $q.blastn.tsv.gz.log 2>&1; \
    done
