# GTDB r214 representative genomes

## Data

    #.genomes: 85205
    #.bases: 273,864,531,576
    #.gzip_size: 75 GB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -O gtdb_repr.lmi --force" > gtdb_repr.lmi.log 2>&1

    elapsed time: 1h:52m:36s
    peak rss: 105.29 GB

    gtdb_repr.lmi: 229.00 GB (228,999,914,466)
     157.29 GB      seeds
      71.71 GB      genomes
       2.13 MB      genomes.map.bin
     160.03 kB      masks.bin
         613 B      info.toml
          48 B      genomes.chunks.bin

Searching

    db=gtdb_repr.lmi

    ls b.*fasta | tac | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # --------------------------------------------------------------------------

    # hits
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk uniq -t -f query,sgenome {} | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       59772
    b.gene_E_coli_16S.fasta.lexicmap.tsv           42363
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         71
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv     413

    # -----------------------------------

    # hits with qcovGnm > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovGnm >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       34399
    b.gene_E_coli_16S.fasta.lexicmap.tsv           36331
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         65
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv       0

    # -----------------------------------

    # hits with qcovHSP > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovHSP >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk summary -t -f hits:count -w 0 | sed 1d' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       34152
    b.gene_E_coli_16S.fasta.lexicmap.tsv           35058
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         60
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv

    b.amr.fasta.lexicmap.tsv                       33882
    b.gene_E_coli_16S.fasta.lexicmap.tsv           35113
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         65
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv

    # --------------------------------------------------------------------------

    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.amr.fasta.lexicmap.tsv.log
    elapsed time: 19m:30s
    peak rss: 5.25 GB

    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 31.943s
    peak rss: 1.57 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 5.101s
    peak rss: 968.42 MB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 46.147s
    peak rss: 1.36 GB


## BLASTN

Indexing

    memusg -t -s "rush -i files.txt 'zcat {}' -j 12 -k --eta | makeblastdb -dbtype nucl -title blastdb -in - -out blastdb/blastdb" > blastdb.log 2>&1

    time: 1 h
     mem: 215 MB
    size: 67 GB

Searching (in different nodes, because blastn cache the index data in memory)

    blastdb=blastdb/gtdb

    q=b.amr.fasta
    q=b.gene_E_coli_16S.fasta
    q=b.gene_E_faecalis_SecY.fasta
    q=b.plasmid_pCUVET18-1784.4.fasta

    memusg -t -s "blastn -num_threads 48 -max_target_seqs 10000000 -db $blastdb -query $q \
            -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
        | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
        > $q.blastn.tsv" > $q.blastn.tsv.log 2>&1;

    # resource
    ls b.*.blastn.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.amr.fasta.blastn.tsv.log
    elapsed time: 12m:01s
    peak rss: 70.52 GB

    b.gene_E_coli_16S.fasta.blastn.tsv.log
    elapsed time: 10m:59s
    peak rss: 70.91 GB

    b.gene_E_faecalis_SecY.fasta.blastn.tsv.log
    elapsed time: 7m:04s
    peak rss: 67.09 GB

    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.log
    elapsed time: 2m:54s
    peak rss: 67.23 GB

Compute the number of genome hits

    # prepare sseqid2ass.tsv.gz
    cat files.txt \
        | rush --eta -k 'seqkit seq -ni {} | awk "{print \$1\"\t{%..}\"}"' \
        | gzip -c > sseqid2ass.tsv.gz

    # preprare a subseq of sseqid2ass.tsv.gz
    for f in b.*.blastn.tsv; do
        csvtk grep -t -P <(csvtk cut -tU -f sseqid $f) \
            sseqid2ass.tsv.gz -o $f.sseqid2ass.tsv.gz
    done

    # filter and add species
    for f in b.*.blastn.tsv; do
        cat $f \
            | csvtk mutate -t -n sgenome -f sseqid \
            | csvtk replace -t -f sgenome -p '(.+)' -r '{kv}' -k $f.sseqid2ass.tsv.gz \
                -o $f.with_sgenome.gz
    done


    # --------------------------------------------------------------------------

    # hits
    ls b.*.blastn.tsv.with_sgenome.gz \
        | rush -k 'echo -ne "{}\t"; \
            csvtk uniq -t -f qseqid,sgenome {} | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.blastn.tsv.with_sgenome.gz                       104240
    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz            40957
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz         121
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz      734


    # hits with qcovs >= 50
    ls b.*.blastn.tsv.with_sgenome.gz \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovs >= 50" {} \
                | csvtk uniq -t -f qseqid,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.blastn.tsv.with_sgenome.gz                       31810
    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz           34650
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz         70
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz       0


    # hits with qcovHSP >= 50
    ls b.*.blastn.tsv.with_sgenome.gz \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$length / \$qlen * 100 >= 50" {} \
                | csvtk uniq -t -f qseqid,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.blastn.tsv.with_sgenome.gz                       30744
    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz           34410
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz         70
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz       0
    
### MMseqs2

Indexing

    memusg -t -s "rush -i testall.txt  'zcat {} ' -j 12 -k --eta | mmseqs createdb --dbtype 2 stdin mmseqs2"

    elapsed time: 23m:28s
    peak rss: 2.15 GB

Searching

    ls b.*.fasta | while read q; do \
        echo $q;
        memusg -t -s "mmseqs easy-search --search-type 3 $q mmseqs2 $q.mmseqs2.tsv tmp/"
    done

    elapsed time: 1h:11m:58s
    peak rss: 888.85 GB
