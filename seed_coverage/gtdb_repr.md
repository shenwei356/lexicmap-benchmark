# GTDB r214 complete genomes


Indexing

    memusg -t -s "lexicmap index -S -X files.txt -O gtdb_repr.lmi.wsp --save-seed-pos --force" > gtdb_repr.lmi.log 2>&1


Retreive seed positions >= 1000 bp.

    cd gtdb_complete.lmi.wsp

    time lexicmap utils seed-pos -d . -a -v -D 1000 -o sd.tsv

    $ cat sd.tsv | csvtk grep -t -f 7 -r -p AAAAAAAAAAAAAAAAA -v | wc -l
    5230


# Numbers of seeds in sliding windows
