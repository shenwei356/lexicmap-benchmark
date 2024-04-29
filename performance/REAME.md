

## Indexing

csvtk csv2md -t performance.indexing.tsv -a l,r,r,,l,r,r,r
|dataset          |genomes  |gzip_size|tool    |db_size|time     |RAM   |
|:----------------|--------:|--------:|:-------|------:|--------:|-----:|
|GTDB complete    |402,538  |578 GB   |LexicMap|510 GB |3 h 26 m |35 GB |
|                 |         |         |Blastn  |360 GB |3 h 11 m |718 MB|
|Genbank+RefSeq   |2,340,672|3.5 TB   |LexicMap|2.91 TB|16 h 40 m|79 GB |
|                 |         |         |Blastn  |2.15 TB|14 h 4 m |4.3 GB|
|AllTheBacteria HQ|1,858,610|3.1 TB   |LexicMap|2.32 TB|10 h 48 m|43 GB |
|                 |         |         |Blastn  |1.76 TB|14 h 3 m |2.9 GB|
|                 |         |         |Phylign |248 GB |/        |/     |

# Searching

Queries

    file                             format  type  num_seqs  sum_len  min_len  avg_len  max_len
    b.gene_E_coli_16S.fasta          FASTA   DNA          1    1,542    1,542    1,542    1,542
    b.gene_E_faecalis_SecY.fasta     FASTA   DNA          1    1,299    1,299    1,299    1,299
    b.plasmid_pCUVET18-1784.4.fasta  FASTA   DNA          1   52,830   52,830   52,830   52,830

Results

|dataset          |genomes  |query          |query_len|tool           |genome_hits|time     |RAM    |
|:----------------|--------:|:--------------|--------:|:--------------|----------:|--------:|------:|
|GTDB complete    |402,538  |a marker gene  |1,299 bp |LexicMap       |2313       |0.9 s    |577 MB |
|                 |         |               |         |Blastn         |7121       |36 m 11 s|351 GB |
|                 |         |a 16S rRNA gene|1,542 bp |LexicMap       |107557     |3 m 05 s |2.7 GB |
|                 |         |               |         |Blastn         |301197     |39 m 13 s|378 GB |
|                 |         |a plasmid      |52,830 bp|LexicMap       |3217       |1 m 10 s |3.2 GB |
|                 |         |               |         |Blastn         |69311      |37 m 42 s|365 GB |
|Genbank+RefSeq   |2,340,672|a marker gene  |1,299 bp |LexicMap       |817        |6.0 s    |1.4 GB |
|                 |         |a 16S rRNA gene|1,542 bp |LexicMap       |832161     |18 m 58 s|8.3 GB |
|                 |         |a plasmid      |52,830 bp|LexicMap       |17822      |5 m 02 s |13.7 GB|
|AllTheBacteria HQ|1,858,610|a marker gene  |1,299 bp |LexicMap       |7936       |16.0 s   |1.1 GB |
|                 |         |               |         |Phylign_local  |7937       |1 h 44 m |27.1 GB|
|                 |         |               |         |Phylign_cluster|7937       |32 m 52 s|/      |
|                 |         |a 16S rRNA gene|1,542 bp |LexicMap       |1031705    |21 m 25 s|8.2 GB |
|                 |         |               |         |Phylign_local  |1032948    |8 h 06 m |28.1 GB|
|                 |         |               |         |Phylign_cluster|1032948    |84 m 30 s|/      |
|                 |         |a plasmid      |52,830 bp|LexicMap       |10875      |4 m 28 s |9.7 GB |
|                 |         |               |         |Phylign_local  |1007       |2 h 50 m |20.6 GB|
|                 |         |               |         |Phylign_cluster|1007       |32 m 23 s|/      |