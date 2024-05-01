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

## LexicMap

Indexing

    lexicmap index -I genomes/ -O db.lmi --save-seed-pos --force

    lexicmap utils seed-pos -d db.lmi/ -a -o db.lmi/seed_pos.tsv --plot-dir db.lmi.dist

Searching

    # generate queries
    for len in 500 1000 1500 2000; do
        ls genomes/*.fna.gz \
            | while read f; do \
                seqkit sliding -s 200 -W $len $f -o $f.qlen$len.fasta.gz; \
            done; \
    done

    for f in genomes/*.fasta.gz; do \
        lexicmap search -d db.lmi/ -w $f -o $f.lexicmap.tsv.gz --log $f.lexicmap.tsv.gz.log ; \
    done

Extracting results

    ls genomes/*.fasta.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'qlen={@qlen(\d+)}' \
            'queries=$(seqkit stats {} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk grep -t -f sgenome -p {ass} {}.lexicmap.tsv.gz \
                     | csvtk mutate -t -p "^(.+)_sliding" -n qseqid \
                     | csvtk filter2 -t -f "\$qseqid==\$sseqid" \
                     | csvtk uniq -t -f query \
                     | csvtk nrow -t); \
             echo -e "{ass}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,qlen,queries,hits \
        | csvtk sort -t -k assembly -k qlen:n \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"lexicmap"' \
        > aligned_fraction.lexicmap.tsv


## Blastn

Indexing

    zcat genomes/*.fna.gz | makeblastdb -dbtype nucl -in - -out blastdb/blastdb -title a

Searching

    for f in genomes/*.fasta.gz; do \
        zcat $f \
            | blastn -num_threads 16 -outfmt 6 -db blastdb/blastdb -query - \
            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore \
            | gzip -c > $f.blastn.tsv.gz; \
    done

Extracting results

    ls genomes/*.fasta.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'qlen={@qlen(\d+)}' \
            'queries=$(seqkit stats {} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk mutate -t -p "^(.+)_sliding" -n qseqid2 {}.blastn.tsv.gz \
                     | csvtk filter2 -t -f "\$qseqid2==\$sseqid" \
                     | csvtk uniq -t -f qseqid \
                     | csvtk nrow -t); \
             echo -e "{ass}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,qlen,queries,hits \
        | csvtk sort -t -k assembly -k qlen:n \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"blastn"' \
        > aligned_fraction.blastn.tsv


## Simulated long-reads

Simulated Oxford Nanopore R10.4.1 long-reads: simulated with [Badread](https://github.com/rrwick/Badread) from the 10 bacterial genomes.

    mkdir -p long-reads; cd long-reads;
    ls ../genomes/*.fna.gz | rush 'ln -s {}'
    cd ..

    # simulate, rename, and filter
    ls long-reads/*.fna.gz | rush --eta 'badread simulate --reference {} --quantity 50x | seqkit seq -m 500 | seqkit replace -p ".+" -r "{%..}_r{nr}" -o {}.fastq.gz'

    $ seqkit stats *.fastq.gz -j 16
    file                             format  type  num_seqs      sum_len  min_len   avg_len  max_len
    GCF_000005845.2.fna.gz.fastq.gz  FASTQ   DNA     15,208  232,031,399      500  15,257.2  125,628
    GCF_000006765.1.fna.gz.fastq.gz  FASTQ   DNA     20,695  313,165,917      500  15,132.4  152,538
    GCF_000006945.2.fna.gz.fastq.gz  FASTQ   DNA     16,259  247,517,766      501  15,223.4  123,040
    GCF_000013425.1.fna.gz.fastq.gz  FASTQ   DNA      9,316  141,033,333      508  15,138.8  159,977
    GCF_000195955.2.fna.gz.fastq.gz  FASTQ   DNA     14,526  220,532,517      500  15,181.9  126,681
    GCF_000240185.1.fna.gz.fastq.gz  FASTQ   DNA     18,876  284,058,647      504  15,048.7  138,364
    GCF_001457635.1.fna.gz.fastq.gz  FASTQ   DNA      6,966  105,545,232      501  15,151.5   94,172
    GCF_008632635.1.fna.gz.fastq.gz  FASTQ   DNA     13,145  198,974,545      500  15,136.9  106,406
    GCF_018885085.1.fna.gz.fastq.gz  FASTQ   DNA     13,477  204,770,825      501  15,194.1  120,610
    GCF_022869705.1.fna.gz.fastq.gz  FASTQ   DNA      9,258  143,349,138      503  15,483.8  114,629

### LexicMap

Searching

    for f in long-reads/*.fastq.gz; do \
        lexicmap search -d db.lmi/ -w $f -o $f.lexicmap.tsv.gz --log $f.lexicmap.tsv.gz.log ; \
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
        | csvtk mutate2 -t -n tool --at 1 -e '"lexicmap"' \
        > long_reads.aligned_fraction.lexicmap.tsv

    csvtk pretty -t long_reads.aligned_fraction.lexicmap.tsv
    tool       assembly          queries   hits    recall
    --------   ---------------   -------   -----   ------
    lexicmap   GCF_000005845.2   15208     14860   97.71
    lexicmap   GCF_000006765.1   20695     20112   97.18
    lexicmap   GCF_000006945.2   16259     15873   97.63
    lexicmap   GCF_000013425.1   9316      9117    97.86
    lexicmap   GCF_000195955.2   14526     14176   97.59
    lexicmap   GCF_000240185.1   18876     18367   97.30
    lexicmap   GCF_001457635.1   6966      6826    97.99
    lexicmap   GCF_008632635.1   13145     12863   97.85
    lexicmap   GCF_018885085.1   13477     13136   97.47
    lexicmap   GCF_022869705.1   9258      9087    98.15

Filter by query coverage (70%)

    ls long-reads/*.fna.gz \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'query={}.fastq.gz' \
            'queries=$(seqkit stats {query} -T | csvtk cut -tU -f num_seqs); \
             hits=$(csvtk grep -t -f sgenome -p {ass} {query}.lexicmap.tsv.gz \
                     | csvtk mutate -t -p "^(.+)_r\d+" -n qseqid \
                     | csvtk filter2 -t -f "\$qseqid==\$sgenome" \
                     | csvtk uniq -t -f query \
                     | csvtk filter2 -t -f "\$qcovHSP>=70" \
                     | csvtk nrow -t); \
             echo -e "{ass}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,queries,hits \
        | csvtk sort -t -k assembly \
        | csvtk mutate2 -t -n recall -e '$hits/$queries*100' \
        | csvtk mutate2 -t -n tool --at 1 -e '"lexicmap"' \
        > long_reads.aligned_fraction_qov_ge70.lexicmap.tsv

    csvtk pretty -t long_reads.aligned_fraction_qov_ge70.lexicmap.tsv
    tool       assembly          queries   hits    recall
    --------   ---------------   -------   -----   ------
    lexicmap   GCF_000005845.2   15208     14666   96.44
    lexicmap   GCF_000006765.1   20695     19706   95.22
    lexicmap   GCF_000006945.2   16259     15668   96.37
    lexicmap   GCF_000013425.1   9316      9018    96.80
    lexicmap   GCF_000195955.2   14526     13953   96.06
    lexicmap   GCF_000240185.1   18876     18072   95.74
    lexicmap   GCF_001457635.1   6966      6740    96.76
    lexicmap   GCF_008632635.1   13145     12683   96.49
    lexicmap   GCF_018885085.1   13477     12995   96.42
    lexicmap   GCF_022869705.1   9258      8991    97.12

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
        | csvtk mutate2 -t -n tool --at 1 -e '"blastn"' \
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
        | csvtk mutate2 -t -n tool --at 1 -e '"blastn"' \
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

