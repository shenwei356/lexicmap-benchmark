# GTDB r214 complete genomes

## Data summary

    #.genomes: 402538
    #.bases: 1,501,880,958,179
    #.gzip_size: 578 GB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -O gtdb_complete.lmi --force" > gtdb_complete.lmi.log 2>&1

    elapsed time: 3h:39m:33s
    peak rss: 37.37 GB

    gtdb_complete.lmi: 509.99 GB
     362.98 GB      genomes
     147.00 GB      seeds
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      269.00 B      info.toml

    # 500 bp
    elapsed time: 5h:15m:04s
    peak rss: 32.0 GB

    gtdb_complete.lmi: 522.58 GB
     362.98 GB      genomes
     159.59 GB      seeds
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      331.00 B      info.toml

Searching

    db=gtdb_complete.lmi

    ls b.*fasta | while read q; do \
        echo $q; memusg -t -s "lexicmap search -d $db $q -n 0 -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; \
    done

    # hits
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; csvtk head -n 1 -t {} | csvtk cut -t -f hits -U;'

    b.gene_E_coli_16S.fasta.lexicmap.tsv    295628
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       3649
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    59211


    # hits with qcovGnm > 50
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; \
        csvtk filter2 -t -f "\$qcovGnm >= 50" {} | csvtk uniq -t -f sgenome | csvtk nrow -t '

    b.gene_E_coli_16S.fasta.lexicmap.tsv    284806
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       3646
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    3245


    # hits with qcovHSP > 50
    ls b.*.lexicmap.tsv | rush -k 'echo -ne "{}\t" ; \
        csvtk filter2 -t -f "\$qcovHSP >= 50" {} | csvtk uniq -t -f sgenome | csvtk nrow -t '

    b.gene_E_coli_16S.fasta.lexicmap.tsv    275454
    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv       3645
    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv    1186

    # resource
    ls b.*.lexicmap.tsv.log | rush -k 'echo {} ; tail -n 3 {};'

    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 1m:17s
    peak rss: 3.98 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 1.555s
    peak rss: 905.57 MB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 35.269s
    peak rss: 4.46 GB

    # 500 bp
    b.gene_E_coli_16S.fasta.lexicmap.tsv.log
    elapsed time: 3m:36s
    peak rss: 4.05 GB

    b.gene_E_faecalis_SecY.fasta.lexicmap.tsv.log
    elapsed time: 4.606s
    peak rss: 773.94 MB

    b.plasmid_pCUVET18-1784.4.fasta.lexicmap.tsv.log
    elapsed time: 1m:05s
    peak rss: 4.38 GB


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
        zcat $f | csvtk filter2 -t -f '$qcovhsp > 0' | csvtk uniq -t -f sgenome | csvtk nrow -t;
    done

    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz     301197
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz        7121
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz     63792


    # hits with qcovs >= 50
    for f in b.*.blastn.tsv.with_sgenome.gz; do
        echo -ne "$f\t";
        zcat $f | csvtk filter2 -t -f '$qcovs >= 50' | csvtk uniq -t -f sgenome | csvtk nrow -t;
    done

    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz      278066
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz 6177
    b.plasmid_pCUVET18-1784.4.fasta.blastn.tsv.with_sgenome.gz      2997


    # hits with qcovHSP >= 50
    for f in b.*.blastn.tsv.with_sgenome.gz; do
        echo -ne "$f\t";
        zcat $f | csvtk filter2 -t -f '$qcovhsp >= 50' | csvtk uniq -t -f sgenome | csvtk nrow -t;
    done

    b.gene_E_coli_16S.fasta.blastn.tsv.with_sgenome.gz      277042
    b.gene_E_faecalis_SecY.fasta.blastn.tsv.with_sgenome.gz 6177
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
