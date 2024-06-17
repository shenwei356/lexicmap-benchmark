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

## Genome sequences

### LexicMap

Indexing

    rm genomes/*.qlen*
    lexicmap index -I genomes/ -O db.lmi --save-seed-pos --force

    lexicmap utils seed-pos -d db.lmi/ -a -o db.lmi/seed_pos.tsv --plot-dir db.lmi.dist --force

Searching

    # generate queries
    for len in 500 1000 1500 2000; do
        ls genomes/*.fna.gz \
            | while read f; do \
                seqkit sliding -s 200 -W $len $f | seqkit grep -i -s -v -p NNNNNNNNNNNNNNNNNNNN -o $f.qlen$len.fasta.gz; \
            done; \
    done

    for f in genomes/*.fasta.gz; do \
        lexicmap search -d db.lmi/ -w $f -o $f.lexicmap.tsv --log $f.lexicmap.tsv.log ; \
    done

Extracting results

    ls genomes/*.fasta.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'qlen={@qlen(\d+)}' \
            'queries=$(seqkit stats {} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk grep -t -f sgenome -p {ass} {}.lexicmap.tsv \
                     | csvtk mutate -t -p "^(.+)_sliding" -n qseqid \
                     | csvtk filter2 -t -f "\$qseqid==\$sseqid" \
                     | csvtk uniq -t -f query \
                     | csvtk nrow -t); \
             echo -e "{ass}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,qlen,queries,hits \
        | csvtk sort -t -k assembly -k qlen:n \
        | csvtk mutate2 -t -n recall -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --at 1 -e '"LexicMap"' \
        > aligned_fraction.lexicmap.tsv


### Blastn

Indexing

    zcat genomes/*.fna.gz | makeblastdb -dbtype nucl -in - -out blastdb/blastdb -title a

Searching

    for f in genomes/*.fasta.gz; do \
        zcat $f \
            | blastn -num_threads 16 -outfmt 6 -db blastdb/blastdb -query - \
            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore \
            > $f.blastn.tsv.gz; \
    done

Extracting results

    ls genomes/*.fasta.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'qlen={@qlen(\d+)}' \
            'queries=$(seqkit stats {} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk mutate -t -p "^(.+)_sliding" -n qseqid2 {}.blastn.tsv \
                     | csvtk filter2 -t -f "\$qseqid2==\$sseqid" \
                     | csvtk uniq -t -f qseqid \
                     | csvtk nrow -t); \
             echo -e "{ass}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,qlen,queries,hits \
        | csvtk sort -t -k assembly -k qlen:n \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"Blastn"' \
        > aligned_fraction.blastn.tsv

### Plot

    # concatenate results
    csvtk concat -t -o aligned_fraction.tsv \
        aligned_fraction.blastn.tsv \
        aligned_fraction.lexicmap.tsv

    Rscript plot.aligned_fraction.R


## Simulated long-reads

Simulated Oxford Nanopore R10.4.1 long-reads: simulated with [Badread](https://github.com/rrwick/Badread) from the 10 bacterial genomes.

    mkdir -p long-reads; cd long-reads;
    ls ../genomes/*.fna.gz | rush 'ln -s {}'
    cd ..

    # simulate, rename, and filter
    ls long-reads/*.fna.gz | rush --eta 'badread simulate --reference {} --quantity 50x | seqkit seq -m 500 | seqkit replace -p ".+" -r "{%..}_r{nr}" -o {}.fastq.gz'

    $ seqkit stats -a *.fastq.gz -j 16
    file                                        format  type  num_seqs      sum_len  min_len   avg_len  max_len       Q1        Q2        Q3  sum_gap     N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)
    long-reads/GCF_000005845.2.fna.gz.fastq.gz  FASTQ   DNA     15,208  232,031,399      500  15,257.2  125,628  5,742.5  11,646.5    20,894        0  22,784    3,093   69.67   50.89    14.13  50.29
    long-reads/GCF_000006765.1.fna.gz.fastq.gz  FASTQ   DNA     20,695  313,165,917      500  15,132.4  152,538  5,740.5    11,699  20,853.5        0  22,535    4,120    69.6   50.85     14.1     65
    long-reads/GCF_000006945.2.fna.gz.fastq.gz  FASTQ   DNA     16,259  247,517,766      501  15,223.4  123,040  5,787.5    11,688  20,825.5        0  22,624    3,311    69.6   50.82    14.11   51.7
    long-reads/GCF_000013425.1.fna.gz.fastq.gz  FASTQ   DNA      9,316  141,033,333      508  15,138.8  159,977  5,850.5  11,677.5  20,594.5        0  22,283    1,952   69.31   50.54    14.06  33.67
    long-reads/GCF_000195955.2.fna.gz.fastq.gz  FASTQ   DNA     14,526  220,532,517      500  15,181.9  126,681    5,707    11,692    20,852        0  22,509    2,962   69.57   50.84    14.09  64.13
    long-reads/GCF_000240185.1.fna.gz.fastq.gz  FASTQ   DNA     18,876  284,058,647      504  15,048.7  138,364  5,567.5    11,354    20,693        0  22,703    3,729   69.74   50.97    14.13  56.29
    long-reads/GCF_001457635.1.fna.gz.fastq.gz  FASTQ   DNA      6,966  105,545,232      501  15,151.5   94,172    5,847  11,672.5    20,899        0  22,473    1,480   69.57   50.77    14.11  40.05
    long-reads/GCF_008632635.1.fna.gz.fastq.gz  FASTQ   DNA     13,145  198,974,545      500  15,136.9  106,406    5,757    11,507    20,743        0  22,717    2,672   69.67   50.88    14.13  39.36
    long-reads/GCF_018885085.1.fna.gz.fastq.gz  FASTQ   DNA     13,477  204,770,825      501  15,194.1  120,610    5,819    11,502    20,717        0  22,654    2,732   69.87   51.04    14.16  29.57
    long-reads/GCF_022869705.1.fna.gz.fastq.gz  FASTQ   DNA      9,258  143,349,138      503  15,483.8  114,629    5,799    11,880    21,188        0  23,215    1,932   69.96   51.13    14.18  38.03


### LexicMap

Searching

    for f in long-reads/*.fastq.gz; do \
        echo $f;
        memusg -t -s "lexicmap search -d db.lmi/ -w $f -o $f.lexicmap.tsv.gz" > $f.lexicmap.tsv.gz.log 2>&1 ; \
    done

Extracting results

    ls long-reads/*.fna.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'query={}.fastq.gz' \
            'queries=$(seqkit stats {query} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk grep -t -f sgenome -p {ass} {query}.lexicmap.tsv.gz \
                     | csvtk mutate -t -p "^(.+)_r\d+" -n qseqid \
                     | csvtk filter2 -t -f "\$qseqid==\$sgenome" \
                     | csvtk uniq -t -f query \
                     | csvtk nrow -t); \
             echo -e "{ass}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,queries,hits \
        | csvtk sort -t -k assembly \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"LexicMap"' \
        > long_reads.aligned_fraction.lexicmap.tsv

    csvtk pretty -t long_reads.aligned_fraction.lexicmap.tsv
    tool       assembly          queries   hits    recall
    --------   ---------------   -------   -----   ------
    LexicMap   GCF_000005845.2   15208     14863   97.73
    LexicMap   GCF_000006765.1   20695     20134   97.29
    LexicMap   GCF_000006945.2   16259     15874   97.63
    LexicMap   GCF_000013425.1   9316      9117    97.86
    LexicMap   GCF_000195955.2   14526     14178   97.60
    LexicMap   GCF_000240185.1   18876     18373   97.34
    LexicMap   GCF_001457635.1   6966      6824    97.96
    LexicMap   GCF_008632635.1   13145     12868   97.89
    LexicMap   GCF_018885085.1   13477     13141   97.51
    LexicMap   GCF_022869705.1   9258      9089    98.17

Filter by query coverage (70%)

    ls long-reads/*.fna.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'query={}.fastq.gz' \
            'queries=$(seqkit stats {query} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk grep -t -f sgenome -p {ass} {query}.lexicmap.tsv.gz \
                     | csvtk mutate -t -p "^(.+)_r\d+" -n qseqid \
                     | csvtk filter2 -t -f "\$qseqid==\$sgenome" \
                     | csvtk uniq -t -f query \
                     | csvtk filter2 -t -f "\$qcovGnm>=70" \
                     | csvtk nrow -t); \
             echo -e "{ass}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,queries,hits \
        | csvtk sort -t -k assembly \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"lexicmap"' \
        > long_reads.aligned_fraction_qcov_ge70.lexicmap.tsv

    csvtk pretty -t long_reads.aligned_fraction_qcov_ge70.lexicmap.tsv
    tool       assembly          queries   hits    recall
    --------   ---------------   -------   -----   ------
    lexicmap   GCF_000005845.2   15208     14848   97.63
    lexicmap   GCF_000006765.1   20695     20077   97.01
    lexicmap   GCF_000006945.2   16259     15856   97.52
    lexicmap   GCF_000013425.1   9316      9111    97.80
    lexicmap   GCF_000195955.2   14526     14161   97.49
    lexicmap   GCF_000240185.1   18876     18342   97.17
    lexicmap   GCF_001457635.1   6966      6821    97.92
    lexicmap   GCF_008632635.1   13145     12850   97.76
    lexicmap   GCF_018885085.1   13477     13127   97.40
    lexicmap   GCF_022869705.1   9258      9084    98.12

### Blastn

Searching

    for f in long-reads/*.fastq.gz; do \
        echo $f;
        zcat $f \
            | seqkit fq2fa \
            | blastn -num_threads 16 -db blastdb/blastdb -query - \
                 -outfmt '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
            | gzip -c > $f.blastn.tsv.gz; \
    done

Extracting results

    # mmaping sequence id to assembly
    ls long-reads/*.fna.gz \
        | rush --eta 'seqkit seq -ni {} | sed "s/$/\t{%..}/"' \
        | pigz -c > sseqid2ass.tsv.gz

    ls long-reads/*.fna.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'query={}.fastq.gz' \
            'queries=$(seqkit stats {query} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk mutate -t -p "^(.+)_r\d+" -n sgenome {query}.blastn.tsv.gz \
                     | csvtk mutate -t -f sseqid -n sgenome2 \
                     | csvtk replace -t -f sgenome2 -p "(.+)" -r "{kv}" -k sseqid2ass.tsv.gz \
                     | csvtk filter2 -t -f "\$sgenome2==\$sgenome" \
                     | csvtk uniq -t -f qseqid \
                     | csvtk nrow -t); \
             echo -e "{ass}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,queries,hits \
        | csvtk sort -t -k assembly  \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"Blastn"' \
        > long_reads.aligned_fraction.blastn.tsv

    tool     assembly          queries   hits    recall
    ------   ---------------   -------   -----   ------
    blastn   GCF_000005845.2   15208     14932   98.19
    blastn   GCF_000006765.1   20695     20299   98.09
    blastn   GCF_000006945.2   16259     15945   98.07
    blastn   GCF_000013425.1   9316      9132    98.02
    blastn   GCF_000195955.2   14526     14235   98.00
    blastn   GCF_000240185.1   18876     18489   97.95
    blastn   GCF_001457635.1   6966      6827    98.00
    blastn   GCF_008632635.1   13145     12909   98.20
    blastn   GCF_018885085.1   13477     13199   97.94
    blastn   GCF_022869705.1   9258      9099    98.28

Filter by query coverage (70%)

    ls long-reads/*.fna.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'query={}.fastq.gz' \
            'queries=$(seqkit stats {query} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk mutate -t -p "^(.+)_r\d+" -n sgenome {query}.blastn.tsv.gz \
                     | csvtk mutate -t -f sseqid -n sgenome2 \
                     | csvtk replace -t -f sgenome2 -p "(.+)" -r "{kv}" -k sseqid2ass.tsv.gz \
                     | csvtk filter2 -t -f "\$sgenome2==\$sgenome" \
                     | csvtk uniq -t -f qseqid \
                     | csvtk filter2 -t -f "\$qcovs>=70" \
                     | csvtk nrow -t); \
             echo -e "{ass}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,queries,hits \
        | csvtk sort -t -k assembly  \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"Blastn"' \
        > long_reads.aligned_fraction_qov_ge70.blastn.tsv

    csvtk pretty -t long_reads.aligned_fraction_qov_ge70.blastn.tsv
    tool     assembly          queries   hits    recall
    ------   ---------------   -------   -----   ------
    blastn   GCF_000005845.2   15208     14931   98.18
    blastn   GCF_000006765.1   20695     20294   98.06
    blastn   GCF_000006945.2   16259     15941   98.04
    blastn   GCF_000013425.1   9316      9130    98.00
    blastn   GCF_000195955.2   14526     14232   97.98
    blastn   GCF_000240185.1   18876     18479   97.90
    blastn   GCF_001457635.1   6966      6826    97.99
    blastn   GCF_008632635.1   13145     12904   98.17
    blastn   GCF_018885085.1   13477     13191   97.88
    blastn   GCF_022869705.1   9258      9096    98.25

### Minimap2

Searching

    for f in long-reads/*.fna.gz; do \
        echo $f;
        minimap2 -ax map-ont $f  $f.fastq.gz | samtools view -b > $f.fastq.gz.minimap.bam; \
    done

Extracting results

    ls long-reads/*.fna.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'query={}.fastq.gz' \
            'queries=$(seqkit stats {query} -T | csvtk cut -tU -f num_seqs); \
             hits=$(samtools view -c -F 4 -F 256 -F 2048 {query}.minimap.bam); \
             echo -e "{ass}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,queries,hits \
        | csvtk sort -t -k assembly  \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"minimap2"' \
        > long_reads.aligned_fraction.minimap2.tsv

    tool       assembly          queries   hits    recall
    --------   ---------------   -------   -----   ------
    minimap2   GCF_000005845.2   15208     14932   98.19
    minimap2   GCF_000006765.1   20695     20298   98.08
    minimap2   GCF_000006945.2   16259     15945   98.07
    minimap2   GCF_000013425.1   9316      9132    98.02
    minimap2   GCF_000195955.2   14526     14237   98.01
    minimap2   GCF_000240185.1   18876     18490   97.96
    minimap2   GCF_001457635.1   6966      6829    98.03
    minimap2   GCF_008632635.1   13145     12909   98.20
    minimap2   GCF_018885085.1   13477     13199   97.94
    minimap2   GCF_022869705.1   9258      9101    98.30

### Plot

Aligned fraction

    # concatenate results
    csvtk concat -t -o long_reads.aligned_fraction.tsv \
        long_reads.aligned_fraction.blastn.tsv \
        long_reads.aligned_fraction.lexicmap.tsv \
        long_reads.aligned_fraction.minimap2.tsv

    Rscript plot.long_reads.aligned_fraction.R

Aligned fraction with qcov >= 70%

    # concatenate results
    csvtk concat -t -o long_reads.aligned_fraction_qcov_ge70.tsv \
        long_reads.aligned_fraction_qcov_ge70.blastn.tsv \
        long_reads.aligned_fraction_qcov_ge70.lexicmap.tsv

    Rscript plot.long_reads.aligned_fraction_qcov_ge70.R

Qcov and pident between Blastn and LexicMap for qcovHSP > 70%


    for f in long-reads/*.fna.gz.fastq.gz; do
        csvtk join -t -f '1,2,3' \
            <(csvtk uniq -t -f qseqid $f.blastn.tsv.gz   | csvtk filter2 -t -f "\$qcovhsp>=70" | csvtk cut -t -f qseqid,sseqid,qlen,qcovhsp,pident) \
            <(csvtk uniq -t -f query  $f.lexicmap.tsv.gz | csvtk filter2 -t -f "\$qcovHSP>=70" | csvtk cut -t -f query,sseqid,qlen,qcovHSP,pident) \
            | csvtk rename -t -f 4-7 -n qcov_blastn,pident_blastn,qcov_lexicmap,pident_lexicmap \
            > $f.tsv
    done

    csvtk concat long-reads/*.fastq.gz.tsv -o qcov_pident.tsv

    Rscript plot.long_reads.qcov_pident.R
