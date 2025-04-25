## AllTheBacteria v0.2 High-quality genomes

## Data

    #.genomes: 1,858,610
    #.bases: 7,493,622,021,123
    #.gzip_size: 3.1 TB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -j 48 -S -X files.txt -b 25000 -O atb_hq.lmi --force" > atb_hq.lmi.log 2>&1

    elapsed time: 1.0days 16h:21m:51s
    peak rss: 97.68 GB

    atb_hq.lmi: 4.30 TB (4,304,515,140,156)
       2.36 TB      seeds
       1.94 TB      genomes
      41.12 MB      genomes.map.bin
     160.03 kB      masks.bin
         619 B      info.toml
          24 B      genomes.chunks.bin

Searching

    db=atb_hq.lmi

    ls b.*fasta | tac | while read q; do \
        if [[ $q =~ "amr" ]]; then debug=""; else debug="--debug"; fi; \
        echo $q; memusg -t -s "lexicmap search -j 48 -d $db $q -o $q.lexicmap.tsv $debug" > $q.lexicmap.tsv.log 2>&1; \
    done


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
            csvtk freq -t -f type -k {}.type > {}.type.freq' --eta

    
### Search with a serial numbers of AMR genomes.

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

Counting

    # preprare a subseq of sseqid2ass.tsv.gz
    ls *.sam_summary.gz \
        | rush --eta 'csvtk grep -t -P <(zcat {} | grep -v ^@ {} | grep -v ^= | cut -f 3) sseqid2ass.tsv.gz -o {}.sseqid2ass.tsv.gz'

    # filter and add species
    ls *.sam_summary.gz \
        | rush --eta 'zcat {} \
                | grep -v ^= \
                | sam2tsv.py \
                | csvtk mutate -t -n sgenome -f target \
                | csvtk replace -t -f sgenome -p "(.+)" -r "{kv}" -k {}.sseqid2ass.tsv.gz -o {}.with_sgenome.gz'

    # --------------------------------------------------------------------------
    # hits
    
    qcov1gene=90
    qcov2gene=50
    qcov1plasmid=70
    qcov2plasmid=30
    pident1=90
    pident2=80
    
    # check similarity
    ls *.sam_summary.gz.with_sgenome.gz \
        | rush -v qcov1gene=$qcov1gene       -v qcov2gene=$qcov2gene \
               -v qcov1plasmid=$qcov1plasmid -v qcov2plasmid=$qcov2plasmid \
               -v pident1=$pident1           -v pident2=$pident2 \
            'if [[ "{%}" =~ "plasmid" ]]; then qcov1={qcov1plasmid}; qcov2={qcov2plasmid}; \
                else qcov1={qcov1gene}; qcov2={qcov2gene}; fi; \
            csvtk uniq -t -f query,sgenome {} \
                | csvtk mutate2 -t -n type -e "\$alignment_length/\$query_length*100>=$qcov1 && \$identity*100>={pident1} ? \
                    \"high\" : (\$alignment_length/\$query_length*100>=$qcov2 && \$identity*100>={pident2} ? \"medium\" : \"low\")" -o {}.type; \
            csvtk freq -t -f type -k {}.type > {}.type.freq' --eta

## Summarise

    # time and memory
    ls b.*.{lexicmap,phylign}.tsv.log \
        | rush -v 'query={@b\.(.+)\.fasta}' -v 'tool={@fasta\.(\w+)}' \
            'time=$(tail -n3 {} | head -n 1 | cut -d " " -f3-); \
              mem=$(tail -n2 {} | head -n 1 | cut -d " " -f3-); \
              echo -e "{query}\t{tool}\t$time\t$mem"; ' \
        | csvtk add-header -t -n query,tool,time,memory \
        | csvtk sort -t -k query -k tool \
        | tee atb.search_time.tsv \
        | csvtk pretty -t -r 3-

    # hits
    ls b.*.type.freq \
        | rush -v 'query={@b\.(.+)\.fasta}' -v 'tool={@fasta\.(\w+)}' \
            'n=$(csvtk summary -t -f frequency:sum -w 0 {} | cut -f 2 | sed 1d); \
            h=$(csvtk grep -t -p high   {} | cut -f 2 | sed 1d); \
            m=$(csvtk grep -t -p medium {} | cut -f 2 | sed 1d); \
            l=$(csvtk grep -t -p low    {} | cut -f 2 | sed 1d); \
            echo -e "{query}\t{tool}\t$n\t$h\t$m\t$l" ' \
        | csvtk add-header -t -n query,tool,sum,high-similarity,medium-similarity,low-similarity \
        | csvtk sort -t -k query -k tool \
        | tee atb.search_hits.tsv \
        | csvtk pretty -t -r 3-
    
    csvtk join -t -f 1,2 atb.search_hits.tsv atb.search_time.tsv \
        | csvtk comma -t -F -f sum,*simi* \
        | csvtk pretty -t -r 3-
                    
    query                     tool              sum   high-similarity   medium-similarity   low-similarity          time     memory
    -----------------------   --------   ----------   ---------------   -----------------   --------------   -----------   --------
    amr                       lexicmap   25,563,227         6,693,084           3,814,828       15,055,315   12h:17m:11s   17.68 GB
    amr                       phylign    11,742,865         5,796,412           2,871,952        3,074,501    2h:36m:08s   85.92 GB
    gene_E_coli_16S           lexicmap    1,857,974           496,867             556,951          804,156       25m:52s   14.82 GB
    gene_E_coli_16S           phylign     1,017,766           483,054             434,105          100,607    2h:10m:33s   77.02 GB
    gene_E_faecalis_SecY      lexicmap       38,062             7,935                  18           30,109        1m:52s    3.77 GB
    gene_E_faecalis_SecY      phylign         7,937             7,935                   1                1       30m:48s   77.62 GB
    plasmid_pCUVET18-1784.4   lexicmap      485,295                25               9,201          476,069       31m:31s   15.16 GB
    plasmid_pCUVET18-1784.4   phylign        46,822                27               9,832           36,963       47m:33s   82.56 GB

## MMseqs2

Indexing

    # The final index would produce a file > 4TB (not supported by our file system)
    
    # split the file list into n parts
    parts=5
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
    elapsed time: 1h:02m:40s
    peak rss: 8.22 GB


    ==> mmseqs2-01.log <==
    elapsed time: 1h:02m:31s
    peak rss: 8.21 GB


    ==> mmseqs2-02.log <==
    elapsed time: 1h:02m:54s
    peak rss: 8.21 GB


    ==> mmseqs2-03.log <==
    elapsed time: 1h:02m:57s
    peak rss: 8.21 GB


    ==> mmseqs2-04.log <==
    elapsed time: 1h:02m:55s
    peak rss: 8.21 GB

    
    
    # size
    ls -l  mmseqs2* | csvtk space2tab | csvtk summary -Ht -f 5:sum -w 0 | csvtk comma -Ht
    7,551,205,276,989

## Minimap2

Indexing

    parts=5
    for f in files.txt.n$parts-*; do
        echo $f;
        i=$(echo $f | cut -d "-" -f 2);
        memusg -t -s "rush -i $f 'zcat {}' -j 12 -k --eta \
            | minimap2 -t 48 -x map-ont -d minimap2-$i.mmi -" > minimap2-$i.mmi.log 2>&1;
    done

    # performance
    tail -n 3 minimap2*.mmi.log
    ==> minimap2-00.mmi.log <==
    elapsed time: 6h:27m:37s
    peak rss: 75.49 GB


    ==> minimap2-01.mmi.log <==
    elapsed time: 6h:27m:09s
    peak rss: 76.46 GB


    ==> minimap2-02.mmi.log <==
    elapsed time: 6h:28m:12s
    peak rss: 61.06 GB


    ==> minimap2-03.mmi.log <==
    elapsed time: 6h:28m:16s
    peak rss: 73.21 GB


    ==> minimap2-04.mmi.log <==
    elapsed time: 6h:25m:39s
    peak rss: 60.36 GB


    # size
    ls -l  minimap2-*.mmi | csvtk space2tab | csvtk summary -Ht -f 5:sum -w 0 | csvtk comma -Ht
    15,541,500,350,868
