## Genbank+RefSeq

## Data

    #.genomes: 2,340,672
    #.bases: 9,192,667,695,899
    #.gzip_size: 3.5 TB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -b 50000 -O genbank_refseq.lmi -G files.txt.big.txt" > genbank_refseq.lmi.log 2>&1

    elapsed time: 21h:28m:46s
    peak rss: 82.14 GB

    genbank_refseq.lmi: 2.91 TB
       2.17 TB      genomes
     764.22 GB      seeds
      55.81 MB      genomes.map.bin
     312.53 KB      masks.bin
      270.00 B      info.toml

Searching

    db=genbank_refseq.lmi

    ls b.*fasta | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -n 0 -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # hits
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; csvtk head -n 1 -t {} | csvtk cut -t -f hits -U;'

    b.gene_E_coli_16S.fasta.lexicmap.tsv    1894943
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       16788
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    495915


    # hits with qcovGnm > 50
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; \
        csvtk filter2 -t -f "\$qcovGnm >= 50" {} | csvtk uniq -t -f sgenome | csvtk nrow -t '


    b.gene_E_coli_16S.fasta.lexicmap.tsv    1632443
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       16759
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    18069


    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; \
        csvtk filter2 -t -f "\$qcovHSP >= 50" {} | csvtk uniq -t -f sgenome | csvtk nrow -t '
    b.gene_E_coli_16S.fasta.lexicmap.tsv    1375561
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       16756
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    6555


    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 5m:32s
    peak rss: 17.25 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 3.119s
    peak rss: 2.12 GB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 3m:50s
    peak rss: 22.46 GB


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
