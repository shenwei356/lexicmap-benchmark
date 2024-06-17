# GTDB r214 complete genomes


Indexing

    lexicmap index -S -X files.txt -O gtdb_complete.lmi.wsp --save-seed-pos --force

    gtdb_complete.lmi.wsp/: 536.84 GB
     388.98 GB      genomes
     147.85 GB      seeds
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      269.00 B      info.toml

Retreive seed positions >= 1000 bp.

    cd gtdb_complete.lmi.wsp

    time lexicmap utils seed-pos -d . -a -v -D 1000 -o sd.tsv

    $ cat sd.tsv | csvtk grep -t -f 7 -r -p AAAAAAAAAAAAAAAAA -v | wc -l
    5230
