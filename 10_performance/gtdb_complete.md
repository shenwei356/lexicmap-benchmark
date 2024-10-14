# GTDB r214 complete genomes

## Data summary

    #.genomes: 402,538
    #.bases: 1,501,880,958,179
    #.gzip_size: 578 GB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -O gtdb_complete.lmi --force" > gtdb_complete.lmi.log 2>&1

    elapsed time: 10h:35m:16s
    peak rss: 63.24 GB

    gtdb_complete.lmi: 906.04 GB
     543.06 GB      seeds
     362.98 GB      genomes
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      330.00 B      info.toml

Searching

    db=gtdb_complete.lmi

    ls b.*fasta | tac | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done



    # --------------------------------------------------------------------------

    # hits
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk uniq -t -f query,sgenome {} | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       3867003
    b.gene_E_coli_16S.fasta.lexicmap.tsv            303925
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         5170
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv     63108

    # -----------------------------------

    # hits with qcovGnm > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovGnm >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       2261989
    b.gene_E_coli_16S.fasta.lexicmap.tsv            288173
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         5143
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv      3235

    # -----------------------------------

    # hits with qcovHSP > 50
    ls b.*.lexicmap.tsv \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovHSP >= 50" {} \
                | csvtk uniq -t -f query,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.lexicmap.tsv                       2228339
    b.gene_E_coli_16S.fasta.lexicmap.tsv            278141
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv         5143
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv      1190

    # --------------------------------------------------------------------------

    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.amr.fasta.lexicmap.tsv.log
    elapsed time: 1h:12m:39s
    peak rss: 16.33 GB

    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 3m:55s
    peak rss: 4.43 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 16.367s
    peak rss: 1.4 GB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 6m:19s
    peak rss: 4.57 GB


## BLASTN

Indexing

    memusg -t -s "rush -i files.txt 'zcat {}' -j 12 -k --eta | makeblastdb -dbtype nucl -title blastdb -in - -out blastdb/blastdb" > blastdb.log 2>&1

    time: 3 h 11m
     mem: 718 MB
    size: 360.17 GB

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
    elapsed time: 1h:18m:06s
    peak rss: 442.07 GB

    b.gene_E_coli_16S.fasta.blastn.tsv.log
    elapsed time: 45m:60s
    peak rss: 378.39 GB

    b.gene_E_faecalis_SecY.fasta.blastn.tsv.log
    elapsed time: 36m:11s
    peak rss: 351.19 GB

    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.log
    elapsed time: 37m:42s
    peak rss: 364.67 GB


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

    b.amr.fasta.blastn.tsv.with_sgenome.gz                       5357772
    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz            301197
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz         7121
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz     69311


    # hits with qcovs >= 50
    ls b.*.blastn.tsv.with_sgenome.gz \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovs >= 50" {} \
                | csvtk uniq -t -f qseqid,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.blastn.tsv.with_sgenome.gz                       2281968
    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz            278066
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz         6177
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz      2997


    # hits with qcovHSP >= 50
    ls b.*.blastn.tsv.with_sgenome.gz \
        | rush -k 'echo -ne "{}\t"; \
            csvtk filter2 -t -f "\$qcovhsp >= 50" {} \
                | csvtk uniq -t -f qseqid,sgenome | csvtk nrow -t' \
        | csvtk pretty -Ht -r 2

    b.amr.fasta.blastn.tsv.with_sgenome.gz                       2240766
    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz            277042
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz         6177
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz      2308

## Data

    time genome_updater.sh -d "refseq,genbank" -g "archaea,bacteria" \
        -f "genomic.fna.gz" -o "GTDB_complete" -M "gtdb" -t 12 -m -L curl -i

    cd GTDB_complete/2024-01-30_19-34-40/

    # ----------------- just in case, check the file integrity -----------------
    genomes=files
    # corrupted files
    # find $genomes -name "*.gz" \
    fd ".gz$" $genomes \
        | rush --eta 'seqkit seq -w 0 {} > /dev/null; if [ $? -ne 0 ]; then echo {}; fi' \
        > failed.txt
    # empty files
    find $genomes -name "*.gz" -size 0 >> failed.txt
    # delete these files
    cat failed.txt | rush '/bin/rm {}'
    # redownload them:
    # run the genome_updater command again
    # ----------------- just in case, check the file integrity -----------------

    # for convenience, creat a new directory containing renamed symbol links pointing to the orginal files
    mkdir gtdb; cd gtdb;
    find ../files -name "*.fna.gz" \
        | rush --eta -v 'id={%@^(..._.........\.\d+)}' \
                     -v 'dir1={%@^(..._...)}' -v 'dir2={%@^..._...(...)}' \
            'mkdir -p {dir1}/{dir2}; cd {dir1}/{dir2}; ln -s ../../{} {id}.fna.gz'
    cd ..

    find gtdb/ -name "*.fna.gz" > gtdb.files.txt

    # mapping file
    cut -f 1,8 assembly_summary.txt > name.map
    cut -f 1,6 assembly_summary.txt > taxid.ncbi.map
    # stats
    cat taxid.ncbi.map | taxonkit reformat -I 2 -f '{s}' > taxid.ncbi.map.species
    csvtk freq -Ht -f 3 -nr taxid.ncbi.map.species > taxid.ncbi.map.species.freq
