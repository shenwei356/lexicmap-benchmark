
## Intro

According to previous test, seed density is most important factor to decide lexicmap's tolerance to mutations.
Besides, the number of probles (masks) affects indexing and searching speed a lot.

Parameters to test.

    -m, --masks int                 ► Number of LexicHash masks. (default 40000)
    -D, --seed-max-desert int       ► Maximum length of sketching deserts, or maximum seed distance.
                                    Deserts with seed distance larger than this value will be filled by
                                    choosing k-mers roughly every --seed-in-desert-dist bases. (default 200)
    -d, --seed-in-desert-dist int   ► Distance of k-mers to fill deserts. (default 50)
  
Parameters combinations (6*3=18):

    (-D, -d): (200, 100), (200, 66), (150, 75), (150, 50), (100, 50), (100,33).

    -m: 10000, 20000, 40000
    
DB and queries: 10 genomes * 3 query lengths (250, 500, 1000) * 21 pidents (80-100) = 630

Runs: 630 * 18 = 11340

Result interpretation:

    # mapping rate
    
    for each number of probes:
        for each query length:
            x: pident
            y: average mapping_rate of 10 genomes
            color: seed density
            

## Indexing

Parameters

    $ cat parameters.tsv
    200     100     10000
    200     66      10000
    150     75      10000
    150     50      10000
    100     50      10000
    100     33      10000
    200     100     20000
    200     66      20000
    150     75      20000
    150     50      20000
    100     50      20000
    100     33      20000
    200     100     40000
    200     66      40000
    150     75      40000
    150     50      40000
    100     50      40000
    100     33      40000

Indexing

    for f in parameters/*.fna.gz; do
        echo $f;
        time cat parameters.tsv \
            | rush --eta -j 1 -v f=$f -v 'db={f}.lexicmap-D{1}-d{2}-m{3}' \
                'memusg -H -t -s "lexicmap index -j 48 {f} -O {db} -D {1} -d {2} -m {3} --force" > {db}.log 2>&1'
    done
    
    
Index size, indexing time and memory

    # parameters/GCF_022869705.1.fna.gz.lexicmap-D200-d66-m10000
    
    ls -d parameters/*.fna.gz.lexicmap* | grep -v log \
        | rush --eta -v 'ass={%@^(.+).fna}' \
                -v 'D={@\-D(\d+)\-}'   -v 'd={@\-d(\d+)\-}' -v 'm={@\-m(\d+)$}'\
                'size=$(du -sh -B 1 --apparent-size {} | cut -f 1); \
                 time=$(tail -n 3 {}.log | head -n 1 | cut -d " " -f 3); \
                 mem=$(tail -n 2 {}.log | head -n 1 | cut -d " " -f 3); \
                 echo -e "{ass}\t{D}\t{d}\t{m}\t$size\t$time\t$mem"' \
        | csvtk add-header -t -n assembly,window,density,probes,size,time,mem \
        | csvtk sort -t -k assembly -k window:n -k density:n -k probes:n \
        > pt.index_info.tsv

Ground-truth positions

    # sequence ID -> sequence length
    ls parameters/*.fna.gz | rush --eta 'seqkit fx2tab -nil {} -o {}.seqlen'
    
    i=1
    n=$(ls parameters/*.fastq.gz | wc -l)
    time ls parameters/*.fastq.gz \
        | while read f; do \
            echo "$i/$n";
            i=$((i+1));
            seqkit seq -n $f \
                | perl -ane '$F[1] =~ /^(.+),([+-])strand,(\d+)\-(\d+)/; print "$F[0]\t$1\t$2\t$3\t$4\n";' \
                | csvtk join -Ht -f "2;1" - $(echo $f | rush 'echo {....}').seqlen \
                | csvtk mutate2 -Ht -w 0 -e '$3 == "+" ? $4+1 : $6-$5+1' \
                | csvtk mutate2 -Ht -w 0 -e '$3 == "+" ? $5 : $6-$4' \
                | csvtk cut -Ht -f 1-3,7,8 \
                | csvtk add-header -t -n query,sseqid,ref_strand,ref_start,ref_end \
                > $f.pos.tsv;
        done

        
## Searching

Submitting 30 jobs, each with 378 runs.

    d=parameters
    
    for g in $d/*.fna.gz; do
        for len in 250 500 1000; do
            slurmzy run -t 48 -c 48 200 $f "for q in $g.i*.q$len.fastq.gz; do \
                ls -d $g.lexicmap* | grep -v log | while read d; do \
                    out=\$q.\$(basename \$d | cut -d . -f 5).tsv.gz; \
                    lexicmap search -d \$d -w \$q -o \$out --log \$out.log; \
                done; \
            done";
        done;
    done

            
Extracting results

    # GCF_000005845.2.fna.gz.i87.q500.fastq.gz.lexicmap-D100-d50-m20000.tsv.gz.log
    
    # match positions
    ls parameters/*.tsv.gz \
        | rush --eta 'cat {} \
            | csvtk cut  -t -f query,sseqid,sstart,send \
            | csvtk join -t -f "query,sseqid;query,sseqid" - {...}.pos.tsv \
            | csvtk filter2 -t -f "(\$sstart >= \$ref_start && \$sstart <= \$ref_end) \
                                || (\$send   >= \$ref_start &&   \$send <= \$ref_end)" \
            > {}.loc.tsv'

    ls parameters/*.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'D={@\-D(\d+)\-}'        -v 'd={@\-d(\d+)\-}'     -v 'm={@\-m(\d+)\.}'\
                -v 'qfile={@(.+).lexicmap}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f query {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t{D}\t{d}\t{m}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,window,density,probes,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly -k window:n -k density:n -k probes:n \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        > pt.mapping_rate.tsv
        
## Plot

    Rscript pt.plot.mapping_rate_index_size.R
    

    
