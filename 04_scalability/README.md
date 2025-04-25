## Data

Note: on EBI server, 48 threads are used, while on my server, 80 was used.

### Genomes for indexing

From Genbank+RefSeq

    #.genomes: 2,340,672
    #.bases: 9,192,667,695,899
    #.gzip_size: 3.5 TB

Sampling

    # shuffle input file list
    shuf files.txt > files.txt.shuf
    
    mkdir -p scalability
    
    cat files.txt.shuf | awk 'NR==1'                      > scalability/ref_n1.txt
    cat files.txt.shuf | awk 'NR>=11      && NR<=20'      > scalability/ref_n10.txt
    cat files.txt.shuf | awk 'NR>=101     && NR<=200'     > scalability/ref_n100.txt
    cat files.txt.shuf | awk 'NR>=1001    && NR<=2000'    > scalability/ref_n1000.txt
    cat files.txt.shuf | awk 'NR>=10001   && NR<=20000'   > scalability/ref_n10000.txt
    cat files.txt.shuf | awk 'NR>=100001  && NR<=200000'  > scalability/ref_n100000.txt
    cat files.txt.shuf | awk 'NR>=1000001 && NR<=2000000' > scalability/ref_n1000000.txt
    
### Queries

Combine the rare gene and 16S rRNA gene

    seqkit seq bench/b.gene_E_faecalis_SecY.fasta bench/b.gene_E_coli_16S.fasta > scalability.query.fasta

    seqkit stats scalability.query.fasta 
    file                     format  type  num_seqs  sum_len  min_len  avg_len  max_len
    scalability.query.fasta  FASTA   DNA          2    2,841    1,299  1,420.5    1,542
    
## LexicMap

Indexing

    threads=48
    ls -r scalability/ref_n*.txt \
        | rush -j 1 --eta -v threads=$threads \
            'memusg -t -s "lexicmap index -j {threads} -S -X {} -O {.}.lexicmap --force" > {.}.lexicmap.log 2>&1'

Searching
    
    q=scalability.query.fasta
    out=scalability.search
    mkdir -p $out
    
    threads=48
    ls -d -r scalability/ref_n*.lexicmap \
        | rush -j 1 --eta -v 'db={}' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
            'memusg -t -H -s "lexicmap search -j {threads} -d {db} {q} --debug -o {out}/{n}.lexicmap.tsv" \
                > {out}/{n}.lexicmap.tsv.log 2>&1'
            
## Blastn

Indexing

    ls -r scalability/ref_n*.txt \
        | rush -j 1 --eta \
            'memusg -t -s "rush -i {} \"zcat {{}}\" -j 12 -k \
                | makeblastdb -dbtype nucl -title blastdb -in - -out {.}.blastn/blastdb" > {.}.blastn.log 2>&1'

Searching
    
    q=scalability.query.fasta
    out=scalability.search
    mkdir -p $out
    
    threads=48
    ls -d -r scalability/ref_n*.blastn \
        | rush -j 1 --eta -v 'db={}/blastdb' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
            'memusg -t -H -s "blastn -num_threads {threads} -max_target_seqs 1000000 -db {db} -query {q} \
                    -outfmt \"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp\" \
                | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
                > {out}/{n}.blastn.tsv" > {out}/{n}.blastn.tsv.log 2>&1'
                
    # word size = 15
    threads=48
    ls -d -r scalability/ref_n*.blastn \
        | rush -j 1 --eta -v 'db={}/blastdb' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
            'memusg -t -H -s "blastn -num_threads {threads} -max_target_seqs 1000000 -db {db} -query {q} \
                    -word_size 15 \
                    -outfmt \"6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp\" \
                | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
                > {out}/{n}.blastn2.tsv" > {out}/{n}.blastn2.tsv.log 2>&1'
                
## MMseqs2

Indexing

    conda activate base
    ls -r scalability/ref_n*.txt \
        | rush -j 1 --eta \
            'mkdir -p {.}.mmseqs2; \
             memusg -t -s "rush -i {} \"zcat {{}}\" -j 12 -k \
                | mmseqs createdb --shuffle 0 --dbtype 2 stdin {.}.mmseqs2/mmseqs2" > {.}.mmseqs2.log 2>&1'
                
Searching
    
    q=scalability.query.fasta
    out=scalability.search
    mkdir -p $out $out/tmp
    
    conda activate base
    threads=48
    ls -d -r scalability/ref_n*.mmseqs2 \
        | rush -j 1 --eta -v 'db={}/mmseqs2' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
            'memusg -t -H -s "mmseqs easy-search --threads {threads} --max-seqs 1000000 --max-seq-len 100000 --search-type 3 \
            --split 0 --split-mode 0 --split-memory-limit 0 \
            --format-mode 4 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,qcov \
            {q} {db} {out}/{n}.mmseqs2.tsv {out}/tmp" > {out}/{n}.mmseqs2.tsv.log 2>&1'
            

## Minimap2

Indexing

    conda activate base
    threads=48
    ls -r scalability/ref_n*.txt \
        | grep -v ref_n1000000 \
        | rush -j 1 --eta -v threads=$threads  \
            'memusg -t -s "rush -i {} \"zcat {{}}\" -j 12 -k \
                | minimap2 -t {threads} -x map-ont -d {.}.minimap2 -" > {.}.minimap2.log 2>&1'

    # for ref_n1000000, which would produce a file > 4TB (not supported by our file system)
    split -n r/3 -d scalability/ref_n1000000.txt scalability/ref_n1000000.txt.n3-
    ls -r scalability/ref_n1000000.txt.n3-* \
        | rush -j 1 --eta -v threads=$threads \
            'memusg -t -s "rush -i {} \"zcat {{}}\" -j 12 -k \
                | minimap2 -t {threads} -x map-ont -d {..}.minimap2-{@-(\d+)$} -" > {..}.minimap2-{@-(\d+)$}.log 2>&1'
    
Searching
    
    q=scalability.query.fasta
    out=scalability.search
    mkdir -p $out
    
    conda activate base
    threads=48
    ls -d -r scalability/ref_n*.minimap2 \
        | grep -v ref_n1000000 \
        | rush -j 1 --eta -v 'db={}' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
            'memusg -t -H -s "minimap2 -t {threads} -a -N 1000000 {db} {q} > {out}/{n}.minimap2.sam" \
                > {out}/{n}.minimap2.sam.log 2>&1'
                
    memusg -t -H -s "ls scalability/ref_n1000000.minimap2-* | grep -v log | while read db; do \
        s=\$(echo \$db | cut -d. -f 2); \
        minimap2 -t $threads -a -N 1000000 \$db $q > $out/1000000.\$s.sam; done" \
        > $out/1000000.minimap2.sam.log 2>&1
    
## COBS

Indexing

    # copy file list in current directory
    cp scalability/ref_n*.txt .

    conda activate cobs
    threads=48
    ls -r ref_n*.txt \
        | rush -j 1 --eta -v threads=$threads \
            'memusg -t -H -s "cobs compact-construct -T {threads} --file-type list -k 31 -C {} scalability/{.}.cobs_compact" \
                > scalability/{.}.cobs_compact.log 2>&1'

Searching
    
    q=scalability.query.fasta
    out=scalability.search
    mkdir -p $out

    conda activate cobs
    threads=48
    ls -d -r scalability/ref_n*.cobs_compact \
        | rush -j 1 --eta -v 'db={}' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
            'memusg -t -H -s "cobs query -T {threads} -t 0.33 -i {db} -f {q} > {out}/{n}.cobs.txt" \
                > {out}/{n}.cobs.txt.log 2>&1'
              

              
## Ropebwt3

Indexing
   
    conda activate ropebwt3
    
    # rm {p}.fmr
    
    threads=48
    ls -r scalability/ref_n*.txt \
        | rush -j 1 --eta -v 'p={.}.ropebwt3/ropebwt3' -v threads=$threads \
            'memusg -t -H -s "\
                mkdir -p {.}.ropebwt3; \
                seqkit seq --skip-file-check -X {} \
                    | ropebwt3 build -t {threads} -bo {p}.fmr -; \
                ropebwt3 ssa -o {p}.fmd.ssa -t {threads} {p}.fmr; \
                ropebwt3 build -i {p}.fmr -do {p}.fmd; \
                seqkit seq --skip-file-check -X {} \
                    | seqtk comp \
                    | cut -f1,2 \
                    | gzip -c > {p}.fmd.len.gz; \
                " \
            > {.}.ropebwt3.log 2>&1'

Searching
    
    q=scalability.query.fasta
    out=scalability.search
    mkdir -p $out

    threads=48
    conda activate ropebwt3
    ls -d -r scalability/ref_n*.ropebwt3 \
        | rush -j 1 --eta -v 'db={}/ropebwt3.fmd' -v q=$q -v out=$out -v 'n={@ref_n(\d+)}' -v threads=$threads \
            'memusg -t -H -s "ropebwt3 sw -t {threads} -e -N 100000 -m 10 {db} {q} > {out}/{n}.ropebwt3.tsv" \
                > {out}/{n}.ropebwt3.tsv.log 2>&1'
 
## Result

indexing: manually type in

searching:

    r=1
    ls -r scalability.search.round$r/*.log \
        | rush 'second=$(tail -n 3 {} | head -n 1 | cut -d " " -f 3); \
                kb=$(tail -n 2 {} | head -n 1 | cut -d " " -f 3); \
                echo -e {%:}"\t"{@(\w+)\.\w+\.log$}"\t$second\t$kb"' \
        | csvtk add-header -t -n genomes,tool,seconds,kb \
        | csvtk sort -t -k tool -k genomes:n  \
        | csvtk replace -t -f tool -p '(.+)' -r '{kv}' -k tools_name.tsv \
        > searching_r$r.tsv
    
    csvtk concat -t \
            <(csvtk mutate2 -t -n round -e '"r1"' searching_r1.tsv) \
            <(csvtk mutate2 -t -n round -e '"r2"' searching_r2.tsv) \
            <(csvtk mutate2 -t -n round -e '"r3"' searching_r3.tsv) \
            <(csvtk mutate2 -t -n round -e '"r4"' searching_r4.tsv) \
        | csvtk sort -t -k tool -k genomes:n -k round:N \
        > searching.tsv

    Rscript scale.plot.R
    
