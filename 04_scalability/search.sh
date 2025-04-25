#!/bin/sh

q=scalability.query.fasta
out=scalability.search
mkdir -p $out $out/tmp

threads=48


# blastn

ls -d -r scalability/ref_n*.blastn \
    | rush -j 1 --eta -v 'db={}/blastdb' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
        'memusg -t -H -s "blastn -num_threads {threads} -max_target_seqs 1000000 -db {db} -query {q} \
                -outfmt \"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp\" \
            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
            > {out}/{n}.blastn.tsv" > {out}/{n}.blastn.tsv.log 2>&1'


# mmseqs2
conda activate base
ls -d -r scalability/ref_n*.mmseqs2 \
    | rush -j 1 --eta -v 'db={}/mmseqs2' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
        'memusg -t -H -s "mmseqs easy-search --threads {threads} --max-seqs 1000000 --max-seq-len 100000 --search-type 3 \
        --split 0 --split-mode 0 --split-memory-limit 0 \
        --format-mode 4 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,qcov \
        {q} {db} {out}/{n}.mmseqs2.tsv {out}/tmp" > {out}/{n}.mmseqs2.tsv.log 2>&1'


# minimap2

conda activate base
ls -d -r scalability/ref_n*.minimap2 \
    | grep -v ref_n1000000 \
    | rush -j 1 --eta -v 'db={}' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
        'memusg -t -H -s "minimap2 -t {threads} -a -N 1000000 {db} {q} > {out}/{n}.minimap2.sam" \
            > {out}/{n}.minimap2.sam.log 2>&1'

memusg -t -H -s "ls scalability/ref_n1000000.minimap2-* | grep -v log | while read db; do \
    s=\$(echo \$db | cut -d. -f 2); \
    minimap2 -t {threads} -a -N 1000000 \$db $q > $out/1000000.\$s.sam; done" \
    > $out/1000000.minimap2.sam.log 2>&1


# cobs

conda activate cobs
ls -d -r scalability/ref_n*.cobs_compact \
    | rush -j 1 --eta -v 'db={}' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
        'memusg -t -H -s "cobs query -T {threads} -t 0.33 -i {db} -f {q} > {out}/{n}.cobs.txt" \
            > {out}/{n}.cobs.txt.log 2>&1'


# ropebwt3

conda activate ropebwt3
ls -d -r scalability/ref_n*.ropebwt3 \
    | rush -j 1 --eta -v 'db={}/ropebwt3.fmd' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
        'memusg -t -H -s "ropebwt3 sw -t {threads} -e -N 100000 -m 10 {db} {q} > {out}/{n}.ropebwt3.tsv" \
            > {out}/{n}.ropebwt3.tsv.log 2>&1'


# lexicmap

ls -d -r scalability/ref_n*.lexicmap \
    | rush -j 1 --eta -v 'db={}' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
        'memusg -t -H -s "lexicmap search -j {threads} -d {db} {q} --debug -o {out}/{n}.lexicmap.tsv" \
            > {out}/{n}.lexicmap.tsv.log 2>&1'
