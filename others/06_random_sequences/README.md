

## References

Random DNA sequences are generated from https://faculty.ucr.edu/~mmaduro/random.htm

    $ seqkit  stats refs/*.fa.gz -a
    file            format  type  num_seqs    sum_len    min_len    avg_len    max_len         Q1         Q2         Q3  sum_gap        N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n
    refs/s01.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  49.99      0
    refs/s02.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  49.99      0
    refs/s03.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  49.99      0
    refs/s04.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  49.97      0
    refs/s05.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  49.96      0
    refs/s06.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  50.01      0
    refs/s07.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  50.03      0
    refs/s08.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  50.04      0
    refs/s09.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  50.02      0
    refs/s10.fa.gz  FASTA   DNA          1  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000  5,000,000        0  5,000,000        1       0       0        0  49.97      0

## Steps

Indexing

    ls refs/*.fa.gz | rush 'lexicmap index {} -O {}.lmi'

Alignment

    ls refs/*.fa.gz \
        | rush -j 1 --eta \
            'seqkit sliding -s 200 -W 500 {} \
                | lexicmap search -d {}.lmi -o {}.lexicmap.tsv.gz --log {}.lexicmap.tsv.gz.log'

Stats

    cat refs/*.log | grep 'queries matched'
    
    13:07:30.838 [INFO] 100.0000% (24998/24998) queries matched
    13:07:40.483 [INFO] 100.0000% (24998/24998) queries matched
    13:07:50.300 [INFO] 100.0000% (24998/24998) queries matched
    13:07:59.915 [INFO] 100.0000% (24998/24998) queries matched
    13:08:09.694 [INFO] 100.0000% (24998/24998) queries matched
    13:08:19.529 [INFO] 100.0000% (24998/24998) queries matched
    13:08:29.156 [INFO] 100.0000% (24998/24998) queries matched
    13:08:38.736 [INFO] 100.0000% (24998/24998) queries matched
    13:08:48.414 [INFO] 100.0000% (24998/24998) queries matched
    13:08:57.955 [INFO] 100.0000% (24998/24998) queries matched

