## Why?

> LexicMap might be more sensitive when searching within frequently sequenced genomes,
> whereas fall short in rare ones because the pool of probes might be biased toward the former. 

We will prove this is not true by simulating queries from frequently sequenced genomes and other genomes and aligning to databases like GTDB or ATB.


## Data

Count species in GTDB genomes

    $ csvtk freq -Ht -f 2 -nr ref2species.tsv > ref2species.tsv.freq
    
    $ head -n 5 ref2species.tsv.freq | csvtk pretty -Ht
    Escherichia coli           33849
    Klebsiella pneumoniae      14975
    Staphylococcus aureus      14959
    Salmonella enterica        13832
    Streptococcus pneumoniae   8895
    
    $ tail -n 10 ref2species.tsv.freq | csvtk pretty -Ht
    Zwartia sp903929925     1
    Zwartia sp903932425     1
    Zymogenus saltonus      1
    Zymomonas mobilis_B     1
    Zymomonas sp020216885   1
    
    $ filter species of different genus
    $ tac ref2species.tsv.freq \
        | csvtk mutate -Ht -p '^(\w)' \
        | csvtk uniq -Ht -f 3 \
        | csvtk cut -Ht -f -3 \
        | head -n 15 \
        | tee rare_species.txt \
        | csvtk pretty -Ht
    Zymomonas sp020216885       1
    Youngiibacter multivorans   1
    Xylophilus sp022702125      1
    Ww85 sp024643445            1
    Vulgatibacter incomptus     1
    Uzinura diaspidicola        1
    Tyzzerella sp944384215      1
    Szabonella alba             1
    Ruthia sp021730205          1
    Quinella sp905236665        1
    Pyruvatibacter ectocarpi    1
    Ozemobacter sp012517355     1
    Nriv7 sp018831625           1
    Myxosarcina sp000756305     1
    Lysobacter_G lycopersici    1

Choose genomes

- Frequently sequenced species

        $ csvtk grep -Ht -f 2 -p 'Escherichia coli' ref2species.tsv \
            | csvtk grep -Ht -f 1 -r -p ^GCF \
            | csvtk sort -Ht -k 1 | head -n 1
        GCF_000005845.2 Escherichia coli

- Less sequenced species

        $ csvtk grep -Ht -f 2 -P <(cut -f 1 rare_species.txt) ref2species.tsv \
            | tee rare_species.txt.ass
        GCF_018831625.1 Nriv7 sp018831625
        GCF_000756305.1 Myxosarcina sp000756305
        GCF_007556775.1 Lysobacter_G lycopersici
        GCF_001263175.1 Vulgatibacter incomptus
        GCF_000689395.1 Pyruvatibacter ectocarpi
        GCA_905236665.1 Quinella sp905236665
        GCA_022702125.1 Xylophilus sp022702125
        GCA_021730205.1 Ruthia sp021730205
        GCA_020216885.1 Zymomonas sp020216885
        GCA_024643445.1 Ww85 sp024643445
        GCA_012517355.1 Ozemobacter sp012517355
        GCF_017873395.1 Youngiibacter multivorans
        GCF_016765775.1 Szabonella alba
        GCA_000331975.1 Uzinura diaspidicola
        GCA_944384215.1 Tyzzerella sp944384215
        
        # genomes without gaps
        $ grep -w -f <(cut -f 1 rare_species.txt.ass) files.txt \
            | seqkit stats -X - -T -a \
            | csvtk filter2 -t -f '$sum_n == 0' \
            | csvtk head -t -n 10 \
            | tee rare_species.txt.ass_no_gap \
            | csvtk comma -t -f 4- \
            | csvtk pretty -t -r 4-
        file                         format  type  num_seqs    sum_len    min_len    avg_len    max_len         Q1         Q2         Q3  sum_gap        N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n
        gtdb/GCA_024643445.1.fna.gz  FASTA   DNA        682  4,066,996      2,003    5,963.3     32,248      3,115      4,820      7,493        0      7,393      176       0       0        0  52.58    392
        gtdb/GCA_905236665.1.fna.gz  FASTA   DNA        308  1,639,765      2,005    5,323.9     31,958      2,790    4,078.5    6,703.5        0      6,522       80       0       0        0  44.01      0
        gtdb/GCF_001263175.1.fna.gz  FASTA   DNA          1  4,350,553  4,350,553  4,350,553  4,350,553  4,350,553  4,350,553  4,350,553        0  4,350,553        1       0       0        0  68.94      0
        gtdb/GCF_017873395.1.fna.gz  FASTA   DNA         77  3,670,578      1,050   47,669.8    312,343      5,017     21,792     77,306        0     96,510       13       0       0        0  44.84      0
        gtdb/GCA_021730205.1.fna.gz  FASTA   DNA          1  1,184,132  1,184,132  1,184,132  1,184,132  1,184,132  1,184,132  1,184,132        0  1,184,132        1       0       0        0  36.58      0
        gtdb/GCF_016765775.1.fna.gz  FASTA   DNA         52  4,219,778        644   81,149.6    961,539      1,474      5,604     57,926        0    531,064        3       0       0        0   64.3      0
        gtdb/GCA_020216885.1.fna.gz  FASTA   DNA        184  1,921,851      2,704   10,444.8     19,915      8,871     10,000     12,713        0     10,000       66       0       0        0  45.89      0
        gtdb/GCA_022702125.1.fna.gz  FASTA   DNA        347  4,681,335      2,501   13,490.9    124,187    4,516.5      8,755     17,850        0     21,858       65       0       0        0   70.6      0
        gtdb/GCA_000331975.1.fna.gz  FASTA   DNA          1    263,431    263,431    263,431    263,431    263,431    263,431    263,431        0    263,431        1       0       0        0   30.2      0
        gtdb/GCA_944384215.1.fna.gz  FASTA   DNA        157  1,915,325      2,604   12,199.5     63,848      4,209      7,205     16,616        0     20,937       32       0       0        0  26.31      0
        
        # copy genomes without gaps
        $ mkdir -p genome_bias/
        $ grep -w -f <(cut -f 1 rare_species.txt.ass_no_gap | sed 1d) files.txt \
            | grep -v GCA_024643445.1 \
            | rush 'cp {} genome_bias/'

## Analysis

Generate reads of coverage 10x.

    cd genome_bias/

    qlen=500
    # mutate
    for i in 100; do
        ls *.fna.gz \
            | rush --eta -v i=$i -v qlen=$qlen \
                'badread simulate --seed 1 --reference {} \
                    --quantity 10x --length {qlen},0 \
                    --identity {i},{i},0 \
                    --error_model random --qscore_model ideal \
                    --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
                    --start_adapter_seq "" --end_adapter_seq "" \
                | seqkit seq -g -m {qlen} \
                | seqkit grep -i -s -v -p NNNNNNNNNNNNNNNNNNNN -o {}.i{i}.q{qlen}.fastq.gz'
    done
    
    cd ..
    
Ground-truth positions

    # sequence ID -> sequence length
    ls genome_bias/*.fna.gz | rush --eta 'seqkit fx2tab -nil {} -o {}.seqlen'
    
    i=1
    n=$(ls genome_bias/*.fastq.gz | wc -l)
    time ls genome_bias/*.fastq.gz \
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

Searching with lexicmap
    
    cd genome_bias/
    
    ln -s ../gtdb_complete.lmi

    ls *.fastq.gz \
        | rush -j 1 --eta -c -C genome_bias.rush \
            'seqkit fq2fa {} \
                | lexicmap search -d gtdb_complete.lmi/ -n 1000 \
                    -o {}.lexicmap.tsv --log {}.lexicmap.tsv.log'
    
    cd ..
    
Checking mapping rate.

    # match positions
    ls genome_bias/*.lexicmap.tsv \
        | rush --eta 'cat {} \
            | csvtk cut  -t -f query,sseqid,sstart,send \
            | csvtk join -t -f "query,sseqid;query,sseqid" - {..}.pos.tsv \
            | awk "(\$sstart >= \$ref_start && \$sstart <= \$ref_end) \
                || (\$send   >= \$ref_start &&   \$send <= \$ref_end)" \
            > {}.loc.tsv'

    # stats
    ls genome_bias/*.lexicmap.tsv \
        | rush --eta -v 'ass={%@^(.+).fna}' -v 'ident={@\.i(\d+)\.}' -v 'qlen={@\.q(\d+)\.}' \
                -v 'qfile={@(.+).lexicmap}' \
            'queries=$(csvtk nrow -t {qfile}.pos.tsv); \
             hits=$(csvtk uniq -t -f query {}.loc.tsv | csvtk nrow -t); \
             echo -e "{ass}\t{ident}\t{qlen}\t$queries\t$hits"' \
        | csvtk add-header -t -n assembly,ident,qlen,queries,hits \
        | csvtk sort -t -k ident:n -k qlen:n -k assembly \
        | csvtk mutate2 -t -n recall -w 4 -e '$hits*100/$queries' \
        | csvtk mutate2 -t -n tool --after queries -e '"LexicMap"' \
        | tee  mapping_rate.lexicmap.tsv \
        | csvtk pretty -t

    assembly          ident   qlen   queries   tool       hits     recall  
    ---------------   -----   ----   -------   --------   ------   --------
    GCA_000331975.1   100     500    5263      LexicMap   5263     100.0000
    GCA_012517355.1   100     500    141361    LexicMap   141361   100.0000
    GCA_020216885.1   100     500    37463     LexicMap   37463    100.0000
    GCA_021730205.1   100     500    23678     LexicMap   23678    100.0000
    GCA_022702125.1   100     500    91873     LexicMap   91873    100.0000
    GCA_905236665.1   100     500    31240     LexicMap   31240    100.0000
    GCF_000756305.1   100     500    141007    LexicMap   141007   100.0000
    GCF_001263175.1   100     500    87009     LexicMap   87009    100.0000
    GCF_016765775.1   100     500    84144     LexicMap   84144    100.0000
    GCF_017873395.1   100     500    73024     LexicMap   73024    100.0000
