# GTDB r214 representative genomes


Indexing

    memusg -t -s "lexicmap index -S -X files.txt -O gtdb_repr.lmi.wsp --save-seed-pos --force" > gtdb_repr.lmi.wsp.log 2>&1

    elapsed time: 1h:45m:38s
    peak rss: 116.66 GB

    gtdb_repr.lmi.wsp: 220.14 GB
     145.79 GB      seeds
      74.34 GB      genomes
       2.03 MB      genomes.map.bin
     312.53 KB      masks.bin
      328.00 B      info.toml

## Seed distances

All seed positions.

    cd gtdb_repr.lmi.wsp

    time lexicmap utils seed-pos -d . -a -o sd.tsv.gz
    # 38 min

Seed positions >= 500 bp.

    cd gtdb_repr.lmi.wsp

    time lexicmap utils seed-pos -d . -a -v -D 500 -o sd.ge500.tsv

    cat sd.ge500.tsv | sed 1d | wc -l
    39739

    # proportion of A's > 0.5
    cat sd.ge500.tsv | sed 1d | awk '$7/$6 > 0.5' | wc -l
    39567

    # left 172 sequences
    cat sd.ge500.tsv | sed 1d | awk '$7/$6 <= 0.5' \
      | awk '{print ">"$1":"$2":"$4":"$5":"$6"\n"$8}' \
      > sd.ge500.tsv.fasta

    # check with dustmasker
    cat sd.ge500.tsv.fasta | dustmasker -outfmt acclist | sed 's/^>//' \
      | csvtk add-header -t -n seq,start,end \
      | csvtk mutate -t -n seqlen -f 1 -p ':(\d+)$' \
      | csvtk mutate2 -t -n lclen -e '$3-$2+1' -w 0 \
      | csvtk mutate2 -t -n lcfrac -e '$lclen/$seqlen * 100' \
      > sd.ge500.tsv.fasta.mask.tsv

    # summarize by sequence
    cat sd.ge500.tsv.fasta.mask.tsv \
      | csvtk summary -t -g seq -f lcfrac:sum \
      | csvtk filter2 -t -f '$lcfrac:sum >= 80' \
      | csvtk nrow -t
    172

## Number of seeds for all genomes

    cat  files.txt.stats.tsv  \
        | csvtk replace -t -p '^.+\/|\.fna\.gz$' \
        | csvtk cut -t -f file,sum_len \
        | csvtk rename -t -f 1,2 -n assembly,size \
        > ass2size.tsv

    memusg -t -s "zcat sd.tsv.gz | csvtk freq -t -nr | csvtk rename -t -f 1,2 -n assembly,seeds > nseeds.tsv"
    # 53 min

    csvtk join -t nseeds.tsv ../ass2size.tsv > ass2seeds2gsize.tsv

    Rscript seeds2gsize.R

## Numbers of seeds in sliding windows


## Mapping rate all non-low-complexicty sliding windows

Test

    # filter out low-complexity sequences
    cat test.fasta | dustmasker -hard_masking -outfmt fasta | seqkit replace -s -p 'N+' | seqkit seq -m 500

    memusg -t -s "head -n 1 files.txt | seqkit sliding -s 200 -W 500 --infile-list - | dustmasker -hard_masking -outfmt fasta | seqkit replace -s -p 'N+' | seqkit seq -m 500 | lexicmap search -d gtdb_repr.lmi.wsp -w -o t.txt"

Splitting genome file lists into 100 chunks

    # filter out not indexed genomes
    lexicmap utils genomes -d gtdb_repr.lmi.wsp/ > genomes.txt

    cp files.txt files.txt0
    cat files.txt0 \
      | csvtk mutate -Ht -p '([^/]+).fna.gz$' \
      | csvtk grep -Ht -f 2 -P genomes.txt  -v \
      > files.txt

    n=100
    split -n r/$n -d files.txt files.n-

Run

    ls files.n-* \
      | grep -v lexicmap \
      | grep -v  -E '\.[eo]$' \
      | while read f; do \
        echo $f; \
        slurmzy run -t 48 -c 48 200 $f "cat $f | seqkit sliding -s 200 -W 500 --infile-list - | dustmasker -hard_masking -outfmt fasta | seqkit replace -s -p 'N+' | seqkit seq -g -m 500 | lexicmap search -w -d gtdb_repr.lmi.wsp -o $f.lexicmap.tsv.gz --quiet --log $f.lexicmap.tsv.gz.log"; \
      done;

Stats

    grep matched *.log > stats.txt

    cat stats.txt \
      | csvtk mutate -Ht -p '(\d+)/' | csvtk mutate -Ht -p '/(\d+)' \
      | csvtk summary -Ht -f 2:sum -f 3:sum -w 0 | csvtk mutate2 -Ht -e '$1/$2*100' -w 6 \
      | csvtk add-header -t -n hists,queries,rate \
      | csvtk pretty -Ht

    hists        queries      rate
    1076276794   1076279675   99.999732


## Length of shared prefix between probes and captured k-mers

    cd gtdb_repr.lmi.wsp

    time lexicmap utils kmers -d . --mask 0 -f -o kmers.tsv.gz
    # 5h10m49.570262211s

    memusg -t -s 'zcat kmers.tsv.gz | sed 1d | cut -f 2,3 | uniq | cut -f 2 | csvtk freq -Ht | csvtk add-header -t -n prefix,freq | csvtk sort -t -k 1:n -o prefix.freq.tsv'
    # elapsed time: 3h:31m:24s
    # peak rss: 57.16 MB

    Rscript prefix.freq.R


## Length of shared prefix between probes and captured k-mers without desert filling

Index

    lexicmap index -S -X files.txt -O gtdb_repr.lmi.wsp-no-df --save-seed-pos --no-desert-filling  --force

Count

    cd gtdb_repr.lmi.wsp-no-df
    seq 1 $(lexicmap utils masks --quiet -d . | wc -l) \
      | rush --eta -k \
          'lexicmap utils kmers --quiet -d . -f -m {} \
              | sed 1d | cut -f 1,2,3 | uniq | cut -f 1,3 | csvtk freq -Ht -f 1,2 \
               | csvtk sort -Ht -k 2:n' \
      | csvtk add-header -t -n mask,prefix,freq \
      > mask_prefix.freq.tsv

Frequency of genomes of each probe.

    csvtk summary -t -g mask -f freq:sum mask_prefix.freq.tsv -w 0 \
        | csvtk rename -t -f 2 -n genomes \
        | csvtk sort -t -k genomes:n \
        > probe2genomes.tsv
    csvtk plot hist -t -f genomes -o probe2genomes.png probe2genomes.tsv \
        --xlab "Number of genomes of a probe" --ylab "Frequency"

## Number of k-mers

kmc requires single-line fasta format ...

    name=atb

    # count k-mers
    memusg -t -s "mkdir -p $name.tmp; kmc -k31 -m300 -fm -ci1 @files.txt $name $name.tmp"
    # elapsed time: 1.0days 12h:06m:41s
    # peak rss: 292.0 G

    # number of k-mers
    time kmc_tools transform $name dump /dev/stdout | wc -l

gtdb_repr     : 242,928,092,491 (   85,205 genomes)
gtdb_complete : 292,087,847,711 (  402,538 genomes)
atb_hq        :  59,449,535,089 (1,858,610 genomes)

unikmer

    name=atb

    mkdir -p $name.tmp.unikmer

    cat files.txt | rush --eta -v 'dir={%@^(.+?)...\.}' -v "outdir=$name.tmp.unikmer" \
        'unikmer count -k 31 -K -s -u {} -o {outdir}/{dir}/{%}' \
        -c -C count.rush

    fd .unik$ $name.tmp.unikmer/ | unikmer merge --skip-flag-check --skip-file-check -i - -u -o $name.unik --verbose
