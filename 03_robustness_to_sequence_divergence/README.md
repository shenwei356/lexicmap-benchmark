## Summary

- Genomes: 10 bacteria
- Query length: 1 kb, 2 kb
- Query number: genome coverage 30X.
- Mutation rate: 0.01, 0.02, ..., 0.20
- Tools: LexicMap, Blastn, Minimap2, COBS.


## Data

Genomes of 10 most frequent bacteria.

|assembly       |taxid|species                   |strain                |num_seqs|genome_size|
|:--------------|----:|:-------------------------|:---------------------|-------:|----------:|
|GCF_000005845.2|562  |Escherichia coli          |K-12 substr. MG1655   |1       |4641652    |
|GCF_000240185.1|573  |Klebsiella pneumoniae     |HS11286               |7       |5682322    |
|GCF_000013425.1|1280 |Staphylococcus aureus     |NCTC 8325             |1       |2821361    |
|GCF_000006945.2|28901|Salmonella enterica       |LT2                   |2       |4951383    |
|GCF_001457635.1|1313 |Streptococcus pneumoniae  |NCTC7465              |1       |2110968    |
|GCF_000195955.2|1773 |Mycobacterium tuberculosis|H37Rv                 |1       |4411532    |
|GCF_000006765.1|287  |Pseudomonas aeruginosa    |PAO1                  |1       |6264404    |
|GCF_008632635.1|470  |Acinetobacter baumannii   |K09-14                |2       |3980230    |
|GCF_018885085.1|1496 |Clostridioides difficile  |S-0253                |2       |4095894    |
|GCF_022869705.1|1351 |Enterococcus faecalis     |PartL-Efaecalis-RM8376|1       |2866948    |


## Generating reads

Using https://github.com/rrwick/Badread

Steps:

    mkdir -p queries; cd queries;
    ls ../genomes/*.fna.gz | rush 'ln -s {}'
    cd ..

    qlen=1000
    # mutate
    for i in $(seq 80 100); do
        ls queries/*.fna.gz \
            | rush --eta -v i=$i -v qlen=$qlen \
                'badread simulate --seed 1 --reference {} \
                    --quantity 30x --length {qlen},0 \
                    --identity {i},{i},0 \
                    --error_model random --qscore_model ideal \
                    --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
                    --start_adapter_seq "" --end_adapter_seq "" \
                | seqkit seq -g -m {qlen} \
                | seqkit grep -i -s -v -p NNNNNNNNNNNNNNNNNNNN -o {}.i{i}.q{qlen}.fastq.gz'
    done

    
Ground-truth positions.

    # sequence ID -> sequence length
    ls queries/*.fna.gz | rush --eta 'seqkit fx2tab -nil {} -o {}.seqlen'
    
    i=1
    n=$(ls queries/*.fastq.gz | wc -l)
    time ls queries/*.fastq.gz \
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
    
    # not used.
    # convert to bed format
    i=1
    n=$(ls queries/*.fastq.gz | wc -l)
    ls queries/*.pos.tsv \
        | while read f; do \
            echo "$i/$n";
            i=$((i+1));
            awk '{print $2"\t"$4-1"\t"$5}' $f >  $(echo $f | rush 'echo {.}').bed;
        done

## LexicMap

Indexing

    # default
    ls queries/*.fna.gz \
        | rush --eta 'lexicmap index {} -O {}.lexicmap --quiet --force'

    # ls -d queries/*.fna.gz.lexicmap \
    #    | rush --eta 'cd {}; lexicmap utils seed-pos -d . -o sd.tsv -a --plot-dir hist --quiet'

Searching

    # default
    ls queries/*.fastq.gz \
        | rush -j 2 --eta -v "index={....}.lexicmap" \
            'lexicmap search -d {index} -w {} -o {}.lexicmap.tsv.gz --log {}.lexicmap.tsv.log'

Extracting results

    # match positions
    ls queries/*.lexicmap.tsv.gz \
        | rush --eta 'cat {} \
            | csvtk cut  -t -f query,sseqid,sstart,send \
            | csvtk join -t -f "query,sseqid;query,sseqid" - {...}.pos.tsv \
            | csvtk filter2 -t -f "(\$sstart >= \$ref_start && \$sstart <= \$ref_end) \
                                || (\$send   >= \$ref_start &&   \$send <= \$ref_end)" \
            > {}.loc.tsv'

    # stats
    ls queries/*.lexicmap.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).lexicmap}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f query {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"LexicMap"' \
        > mapping_rate.lexicmap.tsv

## Blastn

Indexing

    ls queries/*.fna.gz \
        | rush --eta 'zcat {} | makeblastdb -dbtype nucl -in - -out {}.blastn/db -title {%..}'

Searching

    # default word size
    threads=$(grep -c processor /proc/cpuinfo)
    ls queries/*.fastq.gz \
        | rush -j 2 --eta -v threads=$threads -v "index={....}.blastn/db" \
            'zcat {} \
                | seqkit fq2fa \
                | blastn -num_threads {threads} -db {index} -query - \
                    -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp" \
                | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
                -o {}.blastn.tsv.gz'

    # word size 15
    threads=$(grep -c processor /proc/cpuinfo)
    ls queries/*.fastq.gz \
        | rush -j 4 --eta -v threads=$threads -v "index={....}.blastn/db" \
            'zcat {} \
                | seqkit fq2fa \
                | blastn -num_threads {threads} -db {index} -query - \
                    -word_size 15 \
                    -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp" \
                | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
                -o {}.blastn2.tsv.gz'

Extracting results

    # match positions
    ls queries/*.blastn?(?).tsv.gz \
        | rush --eta 'cat {} \
            | csvtk cut  -t -f qseqid,sseqid,sstart,send \
            | csvtk join -t -f "qseqid,sseqid;query,sseqid" - {...}.pos.tsv \
            | csvtk filter2 -t -f "(\$sstart >= \$ref_start && \$sstart <= \$ref_end) \
                                || (\$send   >= \$ref_start &&   \$send <= \$ref_end)" \
            > {}.loc.tsv'

    # stats
    ls queries/*.blastn.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).blastn}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f qseqid {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"Blastn"' \
        > mapping_rate.blastn.tsv
        
    # word size 15
    # stats
    ls queries/*.blastn2.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).blastn2}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f qseqid {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"Blastn(ws=15)"' \
        > mapping_rate.blastn2.tsv

## Minimap2

Searching

    threads=$(grep -c processor /proc/cpuinfo)
    ls queries/*.fastq.gz \
        | rush -j 2 --eta -v threads=$threads -v "index={....}" \
            'minimap2 -t {threads} -ax map-ont {index} {} | samtools view -b > {}.minimap.bam'


Extracting results

    conda activate samtools

    # match positions
    ls queries/*.minimap.bam \
        | rush --eta 'samtools view {} \
            | sam2tsv.py \
            | csvtk cut  -t -f query,target,target_start,target_end \
            | csvtk join -t -f "query,target;query,sseqid" - {..}.pos.tsv \
            | csvtk filter2 -t -f "(\$target_start >= \$ref_start && \$target_start <= \$ref_end) \
                                || (\$target_end   >= \$ref_start &&   \$target_end <= \$ref_end)" \
            > {}.loc.tsv'

    # stats
    ls queries/*.minimap.bam \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).minimap}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f query {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"Minimap2"' \
        > mapping_rate.minimap.tsv

## MMseqs2

Indexing

    ls queries/*.fna.gz \
        | rush --eta 'mkdir -p {}.mmseqs; zcat {} | mmseqs createdb --shuffle 0 --dbtype 2 stdin {}.mmseqs/mmseqs'
        
Searching

    ls queries/*.fastq.gz \
        | rush -j 1 --eta -v "index={....}.mmseqs/mmseqs" \
            'slurmzy run -t 12 -c 12 100 {%} \
            "mkdir -p {}.tmp; mmseqs easy-search --threads 12 --max-seq-len 100000 --search-type 3 \
            --split 0 --split-mode 0 --split-memory-limit 100G \
            --format-mode 4 --format-output \
                query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,qcov \
            {} {index} {}.mmseqs2.tsv {}.tmp"'

Extracting results

    # match positions
    ls queries/*.mmseqs2.tsv \
        | rush --eta 'cat {} \
            | csvtk cut  -t -f query,target,tstart,tend \
            | csvtk join -t -f "query,target;query,sseqid" - {..}.pos.tsv \
            | csvtk filter2 -t -f "(\$tstart >= \$ref_start && \$tstart <= \$ref_end) \
                                || (\$tend   >= \$ref_start &&   \$tend <= \$ref_end)" \
            > {}.loc.tsv'
            
    # stats
    ls queries/*.mmseqs2.tsv \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).mmseqs2}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f query {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"MMseqs2"' \
        > mapping_rate.mmseqs2.tsv

## COBS

Indexing

    conda activate cobs

    ls queries/*.fna.gz \
        | rush --eta 'cobs compact-construct -k 31 {} {}.cobs_compact'

Searching

    conda activate cobs

    threshold=0.33
    ls queries/*.fastq.gz \
        | rush -j 2 --eta -v "index={....}.cobs_compact" -v t=$threshold \
            'cobs query --load-complete -t {t} -i {index} -f {} | gzip -c > {}.cobs.txt.gz'


Extracting results

    ls queries/*.cobs.txt.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).cobs}' \
            'queries=$(seqkit stats {qfile} -T | csvtk cut -tU -f num_seqs); \
             hits=$(zcat {} | grep "^*" | grep "1$" | wc -l); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"COBS"' \
        > mapping_rate.cobs.tsv
        
## Ropebwt3

Searching

    conda activate ropebwt3
    
    ls queries/*.fna.gz \
        | rush --eta -v 'p={}.ropebwt3/ropebwt3'\
                'mkdir -p {}.ropebwt3; \
                seqkit seq {} \
                    | ropebwt3 build -t 48 -bo {p}.fmr -; \
                ropebwt3 ssa -o {p}.fmd.ssa -t 48 {p}.fmr; \
                ropebwt3 build -i {p}.fmr -do {p}.fmd; \
                seqkit seq {} \
                    | seqtk comp \
                    | cut -f1,2 \
                    | gzip -c > {p}.fmd.len.gz;'

Searching

    conda activate ropebwt3
    
    ls queries/*.fastq.gz \
        | rush -j 12 --eta -v "index={....}.ropebwt3/ropebwt3.fmd" \
            'ropebwt3 sw -t 4 -N 25 -m 30 {index} {} > {}.ropebwt3.tsv'

Extracting results

    # match positions
    ls queries/*.ropebwt3.tsv \
        | rush --eta 'cat {} \
            | awk "{print \$1\"\t\"\$6\"\t\"(\$8+1)\"\t\"(\$9+1)}" \
            | csvtk add-header -t -n query,sseqid,sstart,send \
            | csvtk cut  -t -f query,sseqid,sstart,send \
            | csvtk join -t -f "query,sseqid;query,sseqid" - {..}.pos.tsv \
            | csvtk filter2 -t -f "(\$sstart >= \$ref_start && \$sstart <= \$ref_end) \
                                || (\$send   >= \$ref_start &&   \$send <= \$ref_end)" \
            > {}.loc.tsv'
            
    # stats
    ls queries/*.ropebwt3.tsv \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).ropebwt3}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -Ht -f 1 {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"Ropebwt3"' \
        > mapping_rate.ropebwt3.tsv


## Plot
    
    # mapping_rate.blastn2.tsv \
    #    | csvtk sort -t -k assembly -k ident:n -k qlen:n \
    
    # concatenate results
    csvtk concat -t \
        mapping_rate.blastn2.tsv \
        mapping_rate.ropebwt3.tsv \
        mapping_rate.mmseqs2.tsv \
        mapping_rate.minimap.tsv \
        mapping_rate.lexicmap.tsv \
        mapping_rate.blastn.tsv \
        mapping_rate.cobs.tsv \
        | csvtk join -t - <(csvtk cut -t -f assembly,genome_size genome_info.tsv) \
        | csvtk replace -t -f qlen -p '(.+)' -r 'Query length: $1 bp' \
        | csvtk sort -t -k ident:nr -k qlen:N -k tool -k assembly \
        > mapping_rate.tsv

    csvtk pretty -t mapping_rate.tsv  > mapping_rate.tsv.pretty
    
    Rscript plot.mapping_rate.R
