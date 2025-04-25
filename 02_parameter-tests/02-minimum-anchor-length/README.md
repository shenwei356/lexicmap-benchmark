
## LexicMap

Parameters to test

    -p, --seed-min-prefix int            ► Minimum (prefix) length of matched seeds. (default 15)
    -P, --seed-min-single-prefix int     ► Minimum (prefix) length of matched seeds if there's only one
                                         pair of seeds matched. (default 17)

    cat parameters2.tsv
    15      15
    15      17
    17      17
    17      19
    19      19
    19      21

Indexing

    ls parameters2/*.fna.gz \
        | rush --eta 'lexicmap index {} -O {}.lexicmap --quiet --force'

        
Ground-truth positions

    # sequence ID -> sequence length
    ls parameters2/*.fna.gz | rush --eta 'seqkit fx2tab -nil {} -o {}.seqlen'
    
    i=1
    n=$(ls parameters2/*.fastq.gz | wc -l)
    time ls parameters2/*.fastq.gz \
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
        
Searching

    d=parameters2
    
    for g in $d/*.fna.gz; do
        for len in 250 500 1000; do
            slurmzy run -t 48 -c 48 200 $g "cat parameters2.tsv \
                | rush -j 1 -v g=$g -v len=$len \
                    'for q in {g}.i*.q{len}.fastq.gz; do \
                        out=\$q.lexicmap-p{1}-P{2}.tsv.gz; \
                        lexicmap search -d {g}.lexicmap -w \$q -p {1} -P {2} -o \$out --log \$out.log; \
                    done' ";
        done;
    done
    
Extracting results
    
    # GCF_000005845.2.fna.gz.i100.q1000.fastq.gz.lexicmap-p15-P15.tsv.gz

    # match positions
    ls parameters2/*.lexicmap-*.tsv.gz \
        | rush --eta 'cat {} \
            | csvtk cut  -t -f query,sseqid,sstart,send \
            | csvtk join -t -f "query,sseqid;query,sseqid" - {...}.pos.tsv \
            | csvtk filter2 -t -f "(\$sstart >= \$ref_start && \$sstart <= \$ref_end) \
                                || (\$send   >= \$ref_start &&   \$send <= \$ref_end)" \
            > {}.loc.tsv'

    ls parameters2/*.lexicmap-*.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'p={@\-p(\d+)\-}'     -v 'P={@\-P(\d+)\.}'\
                -v 'qfile={@(.+).lexicmap}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f query {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t{p}-{P}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,wordsize,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly -k wordsize:N \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        > pt2.mapping_rate.lexicmap.tsv

## Blastn

Indexing

    ls parameters2/*.fna.gz \
        | rush --eta 'zcat {} | makeblastdb -dbtype nucl -in - -out {}.blastn/db -title {%..}'

Searching

    d=parameters2
    
    for g in $d/*.fna.gz; do
        for len in 250 500 1000; do
            slurmzy run -t 48 -c 48 200 $g "cat parameters2.blastn.tsv \
                | rush -j 1 -v g=$g -v len=$len \
                    'for q in {g}.i*.q{len}.fastq.gz; do \
                        echo \$q; \
                        out=\$q.blastn-w{1}.tsv.gz; \
                        zcat \$q | seqkit fq2fa \
                            | blastn -num_threads 48 -db {g}.blastn/db -word_size {1} -query - \
                                -outfmt \"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp\" \
                            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
                                -o \$out; \
                    done' ";
        done;
    done

Extracting results

    # parameters2/GCF_000005845.2.fna.gz.i100.q1000.fastq.gz.blastn-w15.tsv.gz

    # match positions
    ls parameters2/*.blastn-*.tsv.gz \
        | rush --eta 'cat {} \
            | csvtk cut  -t -f qseqid,sseqid,sstart,send \
            | csvtk join -t -f "qseqid,sseqid;query,sseqid" - {...}.pos.tsv \
            | csvtk filter2 -t -f "(\$sstart >= \$ref_start && \$sstart <= \$ref_end) \
                                || (\$send   >= \$ref_start &&   \$send <= \$ref_end)" \
            > {}.loc.tsv'
            
    ls parameters2/*.blastn-*.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'w={@\-w(\d+)\.}'\
                -v 'qfile={@(.+).blastn}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f qseqid {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t{w}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,wordsize,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly -k wordsize:N \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        > pt2.mapping_rate.blastn.tsv
        
## Plot

Only LexicMap

    Rscript pt2.plot.mapping_rate2.R
    
LexicMap + Blastn

    # csvtk pretty -t mapping_rate.tsv  > mapping_rate.tsv.pretty

