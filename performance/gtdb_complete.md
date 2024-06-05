# GTDB r214 complete genomes

## Data summary

    #.genomes: 402538
    #.bases: 1,501,880,958,179
    #.gzip_size: 578 GB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -O gtdb_complete.lmi --force" > gtdb_complete.lmi.log 2>&1

    elapsed time: 3h:26m:43s
    peak rss: 35.21 GB

    gtdb_complete.lmi: 510.01 GB
     362.98 GB      genomes
     147.03 GB      seeds
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      269.00 B      info.toml


Searching

    db=gtdb_complete.lmi

    ls b.*fasta | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -n 0 -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # hits
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; csvtk head -n 1 -t {} | csvtk cut -t -f hits -U;'
    b.gene_E_coli_16S.fasta.lexicmap.tsv    294285
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       3588
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    58930

    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'
    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 1m:19s
    peak rss: 2.98 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 1.371s
    peak rss: 597.7 MB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 39.758s
    peak rss: 3.11 GB

## BLASTN

Indexing

    memusg -t -s "rush -i files.txt 'zcat {}' -j 12 -k --eta | makeblastdb -dbtype nucl -title blastdb -in - -out blastdb/blastdb" > blastdb.log 2>&1

    time: 3 h 11m
     mem: 718 MB
    size: 360.17 GB

Searching (in different nodes, because blastn cache the index data in memory)

    blastdb=blastdb/gtdb

    q=b.gene_E_coli_16S.fasta
    q=b.gene_E_faecalis_SecY.fasta
    q=b.plasmid_pCUVET18-1784.4.fasta
    memusg -t -s "blastn -num_threads 48 -max_target_seqs 10000000 -db $blastdb -query $q \
            -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
        | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
        > $q.blastn.tsv" > $q.blastn.tsv.log 2>&1;

    # resource
    ls b.*.blastn.tsv.log | rush -k 'echo {} ; tail -n 3 {};'


    # this result is not accurate, the index seems being cached.
    ls b.*.blastn.tsv.log | rush -k 'echo {} ; tail -n 3 {};'
    b.gene_E_coli_16S.fasta.blastn.tsv.log
    elapsed time: 39m:13s
    peak rss: 378.37 GB

    b.gene_E_faecalis_SecY.fasta.blastn.tsv.log
    elapsed time: 36m:11s
    peak rss: 351.19 GB

    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.log
    elapsed time: 37m:42s
    peak rss: 364.67 GB


Compute the number of genome hits


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

    # count
    for f in b.*.blastn.tsv.with_sgenome.gz; do
        echo -ne "$f\t";
        zcat $f | csvtk filter2 -t -f '$qcovhsp >= 50' | csvtk uniq -t -f sgenome | csvtk nrow -t;
    done

    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz     301197
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz        7121
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz     69311


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

    # creat a new directory containing symbol links to the orginal files
    mkdir gtdb; cd gtdb;
    find ../files -name "*.fna.gz" | rush --eta 'ln -s {}'
    brename -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz'
    cd ..

    find gtdb/ -name "*.fna.gz" > gtdb.files.txt

    # mapping file
    cut -f 1,8 assembly_summary.txt > name.map
    cut -f 1,6 assembly_summary.txt > taxid.ncbi.map
    # stats
    cat taxid.ncbi.map | taxonkit reformat -I 2 -f '{s}' > taxid.ncbi.map.species
    csvtk freq -Ht -f 3 -nr taxid.ncbi.map.species > taxid.ncbi.map.species.freq
