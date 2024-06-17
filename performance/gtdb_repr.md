# GTDB r214 representative genomes

## Data

    #.genomes: 85205
    #.bases: 273,864,531,576
    #.gzip_size: 75 GB

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -O gtdb_repr.lmi --force" > gtdb_repr.lmi.log 2>&1

    elapsed time: 1h:04m:27s
    peak rss: 48.05 GB

    gtdb_repr.lmi: 110.33 GB
      66.78 GB      genomes
      43.55 GB      seeds
       2.03 MB      genomes.map.bin
     312.53 KB      masks.bin
      266.00 B      info.toml

Searching

    db=gtdb_repr.lmi

    ls b.* | grep -v lexicmap | while read q; do echo $q; memusg -t -s "lexicmap search -d $db $q -n 0 -o $q.lexicmap.tsv" > $q.lexicmap.tsv.log 2>&1; done

    # hits
    ls b.*.tsv | rush -k 'echo -ne "{}\t" ; csvtk head -n 1 -t {} | csvtk cut -t -f hits -U;'
    b.16S.fasta.lexicmap.tsv        4374
    b.plasmid.fasta.lexicmap.tsv
    b.sm.MutL.fasta.lexicmap.tsv    2

    # resource
    ls b.*.tsv.log | rush -k 'echo {} ; tail -n 3 {};'
    b.16S.fasta.lexicmap.tsv.log
    elapsed time: 38.467s
    peak rss: 747.42 MB

    b.plasmid.fasta.lexicmap.tsv.log
    elapsed time: 13.021s
    peak rss: 768.07 MB

    b.sm.MutL.fasta.lexicmap.tsv.log
    elapsed time: 3.223s
    peak rss: 436.21 MB


## BLASTN

Indexing

    memusg -t -s "rush -i files.txt 'zcat {}' -j 12 -k --eta | makeblastdb -dbtype nucl -title blastdb -in - -out blastdb/blastdb" > blastdb.log 2>&1

    time: 1 h
     mem: 215 MB
    size: 67 GB

Searching

### MMseqs2

Indexing

    memusg -t -s "rush -i testall.txt  'zcat {} ' -j 12 -k --eta | mmseqs createdb --dbtype 2 stdin mmseqs2"

    elapsed time: 23m:28s
    peak rss: 2.15 GB

Searching

    memusg -t -s " mmseqs easy-search --search-type 3 t.sm.MutL.fasta mmseqs2 t.sm.MutL.fasta.mmseqs2.tsv tmp/"

    elapsed time: 1h:11m:58s
    peak rss: 888.85 GB
