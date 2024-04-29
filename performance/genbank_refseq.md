## Genbank+RefSeq

## Data

    #.genomes: 2,340,672
    #.bases: 9,192,667,695,899
    #.gzip_size: 3.5 TB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -b 50000 -O genbank_refseq.lmi -G files.txt.big.txt" > genbank_refseq.lmi.log 2>&1

    elapsed time: 16h:39m:51s
    peak rss: 79.23 GB

    genbank_refseq.lmi: 2.91 TB
       2.17 TB      genomes
     762.72 GB      seeds
      55.81 MB      genomes.map.bin
     312.53 KB      masks.bin
      270.00 B      info.toml


Searching

    db=genbank_refseq.lmi

    ls b.*fasta | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -n 0 -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # hits
    ls b.*.tsv | rush -k 'echo -ne "{}\t" ; csvtk head -n 1 -t {} | csvtk cut -t -f hits -U;'
    b.gene_E_coli_16S.fasta.lexicmap.tsv    832161
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       11773
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    17793

    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'
    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 21m:54s
    peak rss: 8.15 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 16.224s
    peak rss: 1.26 GB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 5m:34s
    peak rss: 12.46 GB

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
