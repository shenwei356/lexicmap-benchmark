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


## LexicMap

Indexing

    # default
    ls queries/*.fna.gz \
        | rush --eta 'lexicmap index {} -O {}.lexicmap --quiet --force --save-seed-pos'

    ls -d queries/*.fna.gz.lexicmap \
        | rush --eta 'cd {}; lexicmap utils seed-pos -d . -o sd.tsv -a --plot-dir hist --quiet'

Searching

    # default
    ls queries/*.fastq.gz \
        | rush -j 2 --eta -v "index={....}.lexicmap" \
            'lexicmap search -d {index} -w {} -o {}.lexicmap.tsv.gz --log {}.lexicmap.tsv.log'

Extracting results

                    | csvtk filter2 -t -f "\$qcovGnm>=70" \

    # default
    ls queries/*.lexicmap.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' \-v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).lexicmap}' \
            'queries=$(seqkit stats {qfile} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk uniq -t -f query {} \
                    | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"LexicMap"' \
        > mapping_rate.lexicmap.tsv

## Blastn

Indexing

    ls queries/*.fna.gz \
        | rush --eta 'zcat {} | makeblastdb -dbtype nucl -in - -out {}.blastn/db -title {%..}'

Searching

    threads=$(grep -c processor /proc/cpuinfo)
    ls queries/*.fastq.gz \
        | rush -j 2 --eta -v threads=$threads -v "index={....}.blastn/db" \
            'zcat {} \
                | seqkit fq2fa \
                | blastn -num_threads {threads} -db {index} -query - \
                    -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp" \
                | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
                -o {}.blastn.tsv.gz'


Extracting results

                     | csvtk filter2 -t -f "\$qcovs>=70" \


    ls queries/*.blastn.tsv.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' \-v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).blastn}' \
            'queries=$(seqkit stats {qfile} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk uniq -t -f qseqid {} \
                     | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"Blastn"' \
        > mapping_rate.blastn.tsv

## Minimap2

Searching

    threads=$(grep -c processor /proc/cpuinfo)
    ls queries/*.fastq.gz \
        | rush -j 2 --eta -v threads=$threads -v "index={....}" \
            'minimap2 -t {threads} -ax map-ont {index} {} | samtools view -b > {}.minimap.bam'


Extracting results

    conda activate base

    ls queries/*.minimap.bam \
        | rush --eta -v 'ass={%@^(.+).fna}' \-v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).minimap}' \
            'queries=$(seqkit stats {qfile} -T | csvtk cut -tU -f num_seqs); \
             hits=$(samtools view -c -F 4 -F 256 -F 2048 {}); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"minimap2"' \
        > mapping_rate.minimap.tsv

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
        | rush --eta -v 'ass={%@^(.+).fna}' \-v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).cobs}' \
            'queries=$(seqkit stats {qfile} -T | csvtk cut -tU -f num_seqs); \
             hits=$(zcat {} | grep "^*" | grep "1$" | wc -l); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"COBS"' \
        > mapping_rate.cobs.tsv

## Plot

    # concatenate results
    csvtk concat -t \
        mapping_rate.blastn.tsv \
        mapping_rate.lexicmap.tsv \
        mapping_rate.minimap.tsv \
        mapping_rate.cobs.tsv \
        | csvtk join -t - <(csvtk cut -t -f assembly,genome_size genome_info.tsv) \
        | csvtk sort -t -k assembly -k ident:n -k tool -k qlen:n \
        | csvtk replace -t -f qlen -p '(.+)' -r 'Query length: $1 bp' \
        > mapping_rate.tsv

    csvtk pretty -t mapping_rate.tsv  > mapping_rate.tsv.pretty

    Rscript plot.mapping_rate.R
