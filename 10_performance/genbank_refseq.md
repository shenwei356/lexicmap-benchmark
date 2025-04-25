## Genbank+RefSeq

## Data

    #.genomes: 2,340,672
    #.bases: 9,192,667,695,899
    #.gzip_size: 3.5 TB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -j 48 -S -X files.txt -b 25000 -O genbank_refseq.lmi" > genbank_refseq.lmi.log 2>&1

    elapsed time: 2.0days 10h:52m:11s
    peak rss: 174.35 GB


    genbank_refseq.lmi: 5.45 TB (5,454,659,703,138)
       3.07 TB      seeds
       2.38 TB      genomes
      58.52 MB      genomes.map.bin
     160.03 kB      masks.bin
       3.68 kB      genomes.chunks.bin
         619 B      info.toml

         
Searching

    db=genbank_refseq.lmi

    ls b.*fasta | tac | while read q; do \
        if [[ $q =~ "amr" ]]; then debug=""; else debug="--debug"; fi; \
        echo $q; memusg -t -s "lexicmap search -j 48 -d $db $q -o $q.lexicmap.tsv $debug" > $q.lexicmap.tsv.log 2>&1; \
    done
    
    # --------------------------------------------------------------------------
    # resource
    
    ls b.*.lexicmap.tsv.log \
        | rush -v 'query={@b\.(.+)\.fasta}' -v 'tool={@fasta\.(\w+)}' \
            'time=$(tail -n3 {} | head -n 1 | cut -d " " -f3-); \
              mem=$(tail -n2 {} | head -n 1 | cut -d " " -f3-); \
              echo -e "{query}\t{tool}\t$time\t$mem"; ' \
        | csvtk add-header -t -n query,tool,time,memory \
        | csvtk sort -t -k query -k tool \
        | tee genbank_refseq.search_time.tsv \
        | csvtk pretty -t -r 3-

    # --------------------------------------------------------------------------
    # hits
    
    qcov1gene=90
    qcov2gene=50
    qcov1plasmid=70
    qcov2plasmid=30
    pident1=90
    pident2=80
    
    # check similarity
    ls b.*.lexicmap.tsv \
        | rush -v qcov1gene=$qcov1gene       -v qcov2gene=$qcov2gene \
               -v qcov1plasmid=$qcov1plasmid -v qcov2plasmid=$qcov2plasmid \
               -v pident1=$pident1           -v pident2=$pident2 \
            'if [[ "{%}" =~ "plasmid" ]]; then qcov1={qcov1plasmid}; qcov2={qcov2plasmid}; \
                else qcov1={qcov1gene}; qcov2={qcov2gene}; fi; \
            csvtk uniq -t -f query,sgenome {} \
                | csvtk mutate2 -t -n type -e "\$qcovHSP>=$qcov1 && \$pident>={pident1} ? \
                    \"high\" : (\$qcovHSP>=$qcov2 && \$pident>={pident2} ? \"medium\" : \"low\")" -o {}.type; \
            csvtk freq -t -f type -k {}.type > {}.type.freq'


    ls b.*.type.freq \
        | rush -v 'query={@b\.(.+)\.fasta}' -v 'tool={@fasta\.(\w+)}' \
            'n=$(csvtk summary -t -f frequency:sum -w 0 {} | cut -f 2 | sed 1d); \
            h=$(csvtk grep -t -p high   {} | cut -f 2 | sed 1d); \
            m=$(csvtk grep -t -p medium {} | cut -f 2 | sed 1d); \
            l=$(csvtk grep -t -p low    {} | cut -f 2 | sed 1d); \
            echo -e "{query}\t{tool}\t$n\t$h\t$m\t$l" ' \
        | csvtk add-header -t -n query,tool,sum,high-similarity,medium-similarity,low-similarity \
        | csvtk sort -t -k query -k tool \
        | tee genbank_refseq.search_hits.tsv \
        | csvtk pretty -t -r 3-
        
    csvtk join -t -f 1,2 genbank_refseq.search_hits.tsv genbank_refseq.search_time.tsv \
        | csvtk comma -t -F -f sum,*simi* \
        | csvtk pretty -t -r 3-

    query                     tool              sum   high-similarity   medium-similarity   low-similarity          time     memory
    -----------------------   --------   ----------   ---------------   -----------------   --------------   -----------   --------
    amr                       lexicmap   30,967,882         7,636,386           4,858,063       18,473,433   15h:52m:08s   24.86 GB
    gene_E_coli_16S           lexicmap    1,955,167           245,884             501,691        1,207,592       32m:59s   11.09 GB
    gene_E_faecalis_SecY      lexicmap       41,718            11,746                 115           29,857        3m:06s    3.97 GB
    plasmid_pCUVET18-1784.4   lexicmap      560,330                96              15,370          544,864       52m:22s   14.48 GB


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

## MMseqs2

Indexing

    # The final index would produce a file > 4TB (not supported by our file system)
    
    # split the file list into n parts
    parts=6
    split -n r/$parts -d files.txt files.txt.n$parts-
    
    for f in files.txt.n$parts-*; do
        echo $f;
        i=$(echo $f | cut -d "-" -f 2);
        memusg -t -s "rush -i $f  'zcat {}' -j 12 -k --eta \
            | mmseqs createdb --shuffle 0 --dbtype 2 stdin mmseqs2-$i" > mmseqs2-$i.log 2>&1;
    done
    
    # performance
    tail -n 3 mmseqs2-*.log
    
    ==> mmseqs2-00.log <==
    elapsed time: 1h:07m:18s
    peak rss: 8.02 GB

    ==> mmseqs2-01.log <==
    elapsed time: 1h:06m:57s
    peak rss: 8.02 GB

    ==> mmseqs2-02.log <==
    elapsed time: 1h:06m:35s
    peak rss: 8.01 GB

    ==> mmseqs2-03.log <==
    elapsed time: 1h:06m:34s
    peak rss: 8.04 GB

    ==> mmseqs2-04.log <==
    elapsed time: 1h:06m:48s
    peak rss: 8.02 GB

    ==> mmseqs2-05.log <==
    elapsed time: 1h:06m:50s
    peak rss: 8.04 GB

    
    # size
    ls -l  mmseqs2* | csvtk space2tab | csvtk summary -Ht -f 5:sum -w 0 | csvtk comma -Ht
    9,259,045,546,795

## Minimap2

Indexing

    parts=6
    for f in files.txt.n$parts-*; do
        echo $f;
        i=$(echo $f | cut -d "-" -f 2);
        memusg -t -s "rush -i $f 'zcat {}' -j 12 -k --eta \
            | minimap2 -t 48 -x map-ont -d minimap2-$i.mmi -" > minimap2-$i.mmi.log 2>&1;
    done

    # performance
    ==> minimap2-00.mmi.log <==
    elapsed time: 7h:03m:31s
    peak rss: 71.66 GB


    ==> minimap2-01.mmi.log <==
    elapsed time: 7h:14m:24s
    peak rss: 76.21 GB


    ==> minimap2-02.mmi.log <==
    elapsed time: 7h:03m:18s
    peak rss: 69.76 GB


    ==> minimap2-03.mmi.log <==
    elapsed time: 7h:04m:12s
    peak rss: 71.68 GB


    ==> minimap2-04.mmi.log <==
    elapsed time: 7h:02m:35s
    peak rss: 70.93 GB


    ==> minimap2-05.mmi.log <==
    elapsed time: 7h:01m:57s
    peak rss: 74.0 GB
    
    # size
    ls -l  minimap2-*.mmi | csvtk space2tab | csvtk summary -Ht -f 5:sum -w 0 | csvtk comma -Ht
    20,282,900,812,679
    
