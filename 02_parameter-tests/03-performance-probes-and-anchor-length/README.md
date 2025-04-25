## Searching performance with different number of probes

- Dataset: 402K GTDB complete genomes
- Probes: 10,000, 20,000, 40,000
- Query: 3 queries

Indexing

    probes index-size indexing-time indexing-mem
    10000   955GB     10h:11m       69.77GB
    20000   972GB     12h:11m       77.03GB
    40000  1029GB     19h:50m       86.59GB

Searching

    # rare gene    
    probes time  mem     hits     hits(qcov>=70)   hits(qcov>=70&&pident>=80)
    10000  67s   1.77GB  6100     6068             2355
    20000  70s   2.42GB  6230     6197*            2356*
    40000  64s   3.96GB  5882     5848             2356*
    
    # 16S gene
    probes time  mem     hits     hits(qcov>=70)   hits(qcov>=70&&pident>=80)
    10000  749s  3.35GB  303037   227095           117104
    20000  723s  4.31GB  303664   227275*          117105
    40000  714s  6.14GB  305360   227264           117107*
    
    # plasmid
    probes time  mem     hits     hits(qcov>=50)   hits(qcov>=50&&pident>=80)
    10000  483s  4.31GB  60800     755              755
    20000  627s  5.86GB  63517    1191*             1191*
    40000  732s  7.52GB  64970    1191*             1191*
    
## Searching performance with different anchor length

    # rare gene
    anchor_lengths time  mem     hits     hits(qcov>=70)   hits(qcov>=70&&pident>=80)
    15-15          330s  2.19GB  6505*    6269*            2356*
    15-17           67s  2.11GB  6230     6197             2356*
    17-17           66s  1.97GB  6211     6178             2356*
    17-19           34s  1.88GB  6111     6081             2356*
    19-19           34s  1.93GB  6084     6055             2355
    19-21           29s  1.90GB  5510     5487             2352

    # 16S gene
    anchor_lengths time  mem       hits     hits(qcov>=70)   hits(qcov>=70&&pident>=80)
    15-15          821s  4.19GB    304867*  227285*          117105*
    15-17          723s  4.29GB    303664   227275           117105*
    17-17          721s  4.18GB    303353   227212           117105*
    17-19          710s  4.11GB    302071   227137           117105*
    19-19          708s  4.11GB    301867   227053           117104
    19-21          700s  4.04GB    300413   226953           117104

    # plasmid
    anchor_lengths time   mem     hits     hits(qcov>=50)   hits(qcov>=50&&pident>=80)
    15-15          1191s  6.06GB  64328    1191             1191
    15-17           627s  6.08GB  63517    1191             1191
    17-17           622s  5.74GB  63489    1191             1191
    17-19           362s  5.85GB  62002    1191             1191
    19-19           363s  6.24GB  61938    1191             1191
    19-21           327s  5.67GB  60756    1191             1191
