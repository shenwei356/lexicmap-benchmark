# GTDB r214 complete genomes

## Data summary

    #.genomes: 402,538
    #.bases: 1,501,880,958,179
    #.gzip_size: 578 GB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -j 48 -S -X files.txt -O gtdb_complete.lmi --force" > gtdb_complete.lmi.log 2>&1

    elapsed time: 9h:38m:55s
    peak rss: 71.55 GB

    gtdb_complete.lmi: 972.10 GB (972,098,200,328)
     582.33 GB      seeds
     389.76 GB      genomes
      10.06 MB      genomes.map.bin
     160.03 kB      masks.bin
         616 B      info.toml
         168 B      genomes.chunks.bin

Searching

    db=gtdb_complete.lmi

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
        | tee lexicmap_search_time.tsv \
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
            csvtk freq -t -f type -k {}.type > {}.type.freq' --eta


## BLASTN

Indexing

    memusg -t -s "rush -i files.txt 'zcat {}' -j 12 -k --eta | makeblastdb -dbtype nucl -title blastdb -in - -out blastdb/blastdb" > blastdb.log 2>&1

    # EBI server
    time: 3 h 11m
     mem: 718 MB
    size: 386.73 GB
    
    # my server
    elapsed time: 4h:52m:58s
        peak rss: 670.1 MB


Searching (in different nodes, because blastn cache the index data in memory)

    blastdb=blastdb/blastdb

    q=b.amr.fasta
    q=b.gene_E_coli_16S.fasta
    q=b.gene_E_faecalis_SecY.fasta
    q=b.plasmid_pCUVET18-1784.4.fasta


    memusg -t -s "blastn -num_threads 48 -max_target_seqs 10000000 -db $blastdb -query $q \
            -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
        | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
        > $q.blastn.tsv" > $q.blastn.tsv.log 2>&1;

    # word size = 15
    memusg -t -s "blastn -num_threads 48 -max_target_seqs 10000000 -db $blastdb -query $q \
            -word_size 15 \
            -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
        | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
        > $q.blastn2.tsv" > $q.blastn2.tsv.log 2>&1;
        
    # --------------------------------------------------------------------------
    # resource

    ls b.*.{blastn,blastn2}.tsv.log \
        | rush -v 'query={@b\.(.+)\.fasta}' -v 'tool={@fasta\.(\w+)}' \
            'time=$(tail -n3 {} | head -n 1 | cut -d " " -f3-); \
              mem=$(tail -n2 {} | head -n 1 | cut -d " " -f3-); \
              echo -e "{query}\t{tool}\t$time\t$mem"; ' \
        | csvtk add-header -t -n query,tool,time,memory \
        | csvtk sort -t -k query -k tool \
        | tee lexicmap_search_time.tsv \
        | csvtk pretty -t -r 3-


Compute the number of genome hits

    # prepare sseqid2ass.tsv.gz
    cat files.txt \
        | rush --eta -k 'seqkit seq -ni {} | awk "{print \$1\"\t{%..}\"}"' \
        | gzip -c > sseqid2ass.tsv.gz

    # preprare a subseq of sseqid2ass.tsv.gz
    ls b.*.blastn.tsv \
        | rush --eta 'csvtk grep -t -P <(csvtk cut -tU -f sseqid {} | csvtk uniq -Ht) sseqid2ass.tsv.gz -o {}.sseqid2ass.tsv.gz'

    # filter and add species
    ls b.*.blastn.tsv \
        | rush --eta 'csvtk mutate -t -n sgenome -f sseqid {} \
                | csvtk replace -t -f sgenome -p "(.+)" -r "{kv}" -k {}.sseqid2ass.tsv.gz -o {}.with_sgenome.gz'

    # --------------------------------------------------------------------------

    qcov1gene=90
    qcov2gene=50
    qcov1plasmid=70
    qcov2plasmid=30
    pident1=90
    pident2=80

    # check similarity
    ls b.*.blastn.tsv.with_sgenome.gz \
        | rush -v qcov1gene=$qcov1gene       -v qcov2gene=$qcov2gene \
               -v qcov1plasmid=$qcov1plasmid -v qcov2plasmid=$qcov2plasmid \
               -v pident1=$pident1           -v pident2=$pident2 \
            'if [[ "{%}" =~ "plasmid" ]]; then qcov1={qcov1plasmid}; qcov2={qcov2plasmid}; \
                else qcov1={qcov1gene}; qcov2={qcov2gene}; fi; \
            csvtk uniq -t -f qseqid,sgenome {} \
                | csvtk mutate2 -t -n type -e "\$length/\$qlen*100>=$qcov1 && \$pident>={pident1} ? \
                    \"high\" : (\$length/\$qlen*100>=$qcov2 && \$pident>={pident2} ? \"medium\" : \"low\")" -o {}.type; \
            csvtk freq -t -f type -k {}.type > {}.type.freq' --eta

    
## MMseqs2

Indexing

    memusg -t -s "rush -i files.txt  'zcat {} ' -j 12 -k --eta | mmseqs createdb --shuffle 0 --dbtype 2 stdin mmseqs2"
    
    # EBI server
    elapsed time: 54m:46s
    peak rss: 7.59 GB    
    Size: 1.56TB (1,555,361,328,703)
    
    # my server
    elapsed time: 2h:42m:45s
    peak rss: 7.55 GB    

Searching

    ls b.*fasta | tac | while read q; do \
        echo $q; \
        memusg -t -s "mmseqs easy-search --threads 48 --max-seqs 500000 --max-seq-len 100000 --search-type 3 \
            --split 0 --split-mode 0 --split-memory-limit 0 \
            --format-mode 4 --format-output \
                query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,qcov \
            $q ../mmseqs2 $q.mmseqs2.tsv tmp" > $q.mmseqs2.tsv.log 2>&1;
    done

Compute the number of genome hits

    # preprare a subseq of sseqid2ass.tsv.gz
    ls b.*.mmseqs2.tsv \
        | rush --eta 'csvtk grep -t -P <(csvtk cut -tU -f target {} | csvtk uniq -Ht) sseqid2ass.tsv.gz -o {}.sseqid2ass.tsv.gz'

    # filter and add species
    ls b.*.mmseqs2.tsv \
        | rush --eta 'csvtk mutate -t -n sgenome -f target {} \
                | csvtk replace -t -f sgenome -p "(.+)" -r "{kv}" -k {}.sseqid2ass.tsv.gz -o {}.with_sgenome.gz'

    # --------------------------------------------------------------------------

    qcov1gene=90
    qcov2gene=50
    qcov1plasmid=70
    qcov2plasmid=30
    pident1=90
    pident2=80
    
    # check similarity
    ls b.*.mmseqs2.tsv.with_sgenome.gz \
        | rush -v qcov1gene=$qcov1gene       -v qcov2gene=$qcov2gene \
               -v qcov1plasmid=$qcov1plasmid -v qcov2plasmid=$qcov2plasmid \
               -v pident1=$pident1           -v pident2=$pident2 \
            'if [[ "{%}" =~ "plasmid" ]]; then qcov1={qcov1plasmid}; qcov2={qcov2plasmid}; \
                else qcov1={qcov1gene}; qcov2={qcov2gene}; fi; \
            csvtk uniq -t -f query,sgenome {} \
                | csvtk mutate2 -t -n type -e "\$qcov*100>=$qcov1 && \$pident>={pident1} ? \
                    \"high\" : (\$qcov*100>=$qcov2 && \$pident>={pident2} ? \"medium\" : \"low\")" -o {}.type; \
            csvtk freq -t -f type -k {}.type > {}.type.freq' --eta


## Minimap2

Indexing

    memusg -t -s "rush -i files.txt  'zcat {} ' -j 12 -k --eta | minimap2 -t 48 -x map-ont -d gtdb_complete.mmi - "
    
    # EBI
    elapsed time: 6h:51m:09s
    peak rss: 69.69 GB
    
    # my server
    elapsed time: 12h:25m:05s
    peak rss: 76.01 GB


    gtdb_complete.mmi: 3.41 TB (3,409,046,856,832)
      3.41 TB      gtdb_complete.mmi
    
Searching

    ls b.*fasta | tac | while read q; do \
        echo $q;
        memusg -t -s "minimap2 -t 48 -a -N 10000000 ../gtdb_complete.mmi $q > $q.minimap2.sam" > $q.minimap2.tsv.log 2>&1;
    done
    
Compute the number of genome hits

    # cat files.txt | rush --eta -k 'seqkit fx2tab -nil ../{}' > seqid2len.tsv

    ls b.*fasta | rush 'seqkit fx2tab -nil {} > {}.qlen'

    # preprare a subseq of sseqid2ass.tsv.gz
    ls b.*.minimap2.sam \
        | rush --eta 'csvtk grep -t -P <(grep -v ^@ $f | cut -f 3) sseqid2ass.tsv.gz -o {}.sseqid2ass.tsv.gz'

    # filter and add species
    ls b.*.minimap2.sam \
        | rush --eta 'cat {} \
                | sam2tsv.py \
                | csvtk mutate -t -n sgenome -f target \
                | csvtk replace -t -f sgenome -p "(.+)" -r "{kv}" -k {}.sseqid2ass.tsv.gz -o {}.with_sgenome.gz'

    # --------------------------------------------------------------------------

    qcov1gene=90
    qcov2gene=50
    qcov1plasmid=70
    qcov2plasmid=30
    pident1=90
    pident2=80
    
    # check similarity
    ls b.*.minimap2.sam.with_sgenome.gz \
        | rush -v qcov1gene=$qcov1gene       -v qcov2gene=$qcov2gene \
               -v qcov1plasmid=$qcov1plasmid -v qcov2plasmid=$qcov2plasmid \
               -v pident1=$pident1           -v pident2=$pident2 \
               -v 'qlen={@(.+).minimap2}.qlen' \
            'if [[ "{%}" =~ "plasmid" ]]; then qcov1={qcov1plasmid}; qcov2={qcov2plasmid}; \
                else qcov1={qcov1gene}; qcov2={qcov2gene}; fi; \
            csvtk uniq -t -f query,sgenome {} \
                | csvtk cut -t -f -query_length \
                | csvtk mutate -t -f query -n query_length --after target \
                | csvtk replace -t -f query_length -p "(.+)" -r "{kv}" -k {qlen} \
                | csvtk mutate2 -t -n type -e "\$alignment_length/\$query_length*100>=$qcov1 && \$identity*100>={pident1} ? \
                    \"high\" : (\$alignment_length/\$query_length*100>=$qcov2 && \$identity*100>={pident2} ? \"medium\" : \"low\")" -o {}.type; \
            csvtk freq -t -f type -k {}.type > {}.type.freq' --eta


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

    
## Summarise

    # time and memory
    ls b.*.{lexicmap,blastn,blastn2,mmseqs2,minimap2}.tsv.log \
        | rush -v 'query={@b\.(.+)\.fasta}' -v 'tool={@fasta\.(\w+)}' \
            'time=$(tail -n3 {} | head -n 1 | cut -d " " -f3-); \
              mem=$(tail -n2 {} | head -n 1 | cut -d " " -f3-); \
              echo -e "{query}\t{tool}\t$time\t$mem"; ' \
        | csvtk add-header -t -n query,tool,time,memory \
        | csvtk sort -t -k query -k tool \
        | tee gtdb_complete.search_time.tsv \
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
        | tee gtdb_complete.search_hits.tsv \
        | csvtk pretty -t -r 3-
    
    csvtk join -t -f 1,2 gtdb_complete.search_hits.tsv gtdb_complete.search_time.tsv \
        | csvtk comma -t -F -f sum,*simi* \
        | csvtk pretty -t -r 3-
                    
    query                     tool              sum   high-similarity   medium-similarity   low-similarity          time      memory
    -----------------------   --------   ----------   ---------------   -----------------   --------------   -----------   ---------
    amr                       blastn      5,357,772         1,150,407             772,858        3,434,507    1h:18m:06s   442.07 GB
    amr                       blastn2    10,877,544         1,150,410             840,464        8,886,670    1h:16m:01s   311.93 GB
    amr                       lexicmap    4,665,317         1,123,251             776,153        2,765,913    2h:23m:40s    10.69 GB
    amr                       minimap2    2,078,490           943,516             319,529          815,445    5h:07m:50s    20.22 GB
    amr                       mmseqs2    10,137,345         1,148,942             808,177        8,180,226   10h:34m:18s   406.87 GB
    gene_E_coli_16S           blastn        301,197            61,878             109,477          129,842       45m:60s   378.39 GB
    gene_E_coli_16S           blastn2       301,197            61,878             109,477          129,842       54m:51s   378.42 GB
    gene_E_coli_16S           lexicmap      306,064            60,999              69,293          175,772        5m:03s     5.19 GB
    gene_E_coli_16S           minimap2       17,656            15,998               1,652                6    4h:48m:33s    20.22 GB
    gene_E_coli_16S           mmseqs2       324,364            60,915              89,874          173,575    8h:39m:00s   400.65 GB
    gene_E_faecalis_SecY      blastn          7,121             2,311                  47            4,763       36m:11s   351.19 GB
    gene_E_faecalis_SecY      blastn2        57,741             2,311                  47           55,383       52m:51s   324.06 GB
    gene_E_faecalis_SecY      lexicmap        6,255             2,311                  46            3,898       29.645s     2.06 GB
    gene_E_faecalis_SecY      minimap2        2,312             2,312                                         4h:46m:48s    20.22 GB
    gene_E_faecalis_SecY      mmseqs2        67,537             2,304                  54           65,179    7h:16m:14s   400.66 GB
    plasmid_pCUVET18-1784.4   blastn         69,311                21               2,865           66,425       37m:42s   364.67 GB
    plasmid_pCUVET18-1784.4   blastn2        91,847                21               2,865           88,961       51m:22s   142.81 GB
    plasmid_pCUVET18-1784.4   lexicmap       65,029                21               2,808           62,200        6m:59s     6.77 GB
    plasmid_pCUVET18-1784.4   minimap2        3,033                35               1,873            1,125    5h:28m:35s    20.22 GB
    plasmid_pCUVET18-1784.4   mmseqs2        90,277                 7               1,650           88,620   12h:25m:10s   400.74 GB



