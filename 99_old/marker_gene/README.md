## Taxonomy

    wget https://github.com/shenwei356/gtdb-taxdump/releases/download/v0.4.0/gtdb-taxdump.tar.gz
    tar -zxvf gtdb-taxdump.tar.gz
    cp -r gtdb-taxdump/R214/ taxdump

    cat taxdump/taxid.map | taxonkit reformat --data-dir taxdump/ -I 2 -f '{s}' > ass2taxid2species.tsv

Top 20 Bacteria species

    echo Bacteria | taxonkit name2taxid
    Bacteria        81602897

    cat ass2taxid2species.tsv \
        | csvtk grep -Ht -f 2 -P <(taxonkit list --data-dir taxdump/ -I "" --ids 81602897) \
        | csvtk freq -Ht -f 3 -nr \
        | csvtk head -n 20 \
        | tee top20.bac.tsv | csvtk pretty -Ht
    Escherichia coli             33849
    Klebsiella pneumoniae        14975
    Staphylococcus aureus        14959
    Salmonella enterica          13832
    Streptococcus pneumoniae     8895
    Mycobacterium tuberculosis   7132
    Pseudomonas aeruginosa       7037
    Acinetobacter baumannii      6912
    Clostridioides difficile     2701
    Enterococcus_B faecium       2657
    Enterobacter hormaechei_A    2605
    Campylobacter_D jejuni       2442
    Enterococcus faecalis        2314
    Listeria monocytogenes       2307
    Neisseria meningitidis       2243
    Streptococcus pyogenes       2220
    Listeria monocytogenes_B     1985
    Mycobacterium abscessus      1903
    Vibrio parahaemolyticus      1892
    Burkholderia mallei          1824
    Vibrio cholerae              1663


## Marker genes

    # complete
    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.1/genomic_files_all/ar53_marker_genes_all_r214.tar.gz
    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/214.1/genomic_files_all/bac120_marker_genes_reps_r214.tar.gz

    for f in *.tar.gz; do tar -zxvf $f; done

    # remove faa
    rm -rf */faa/

    # compress theme to save space
    fd .fna  | rush --eta 'gzip {}'

Common genes

    fd .fna.gz | rush 'echo {%..}' | csvtk freq -Ht -nr | csvtk filter2 -Ht -f '$2 > 1' | csvtk cut -Ht -f 1 > common.txt
    PF00410.20
    PF00466.21
    TIGR00064
    TIGR00967
    TIGR01171

    $ fd .fna | grep -f common.txt | seqkit stats -X -
    file                                                format  type  num_seqs      sum_len  min_len  avg_len  max_len
    ar53_marker_genes_all_r214/fna/PF00410.20.fna.gz    FASTA   DNA      6,988    2,725,422       84      390    1,140
    ar53_marker_genes_all_r214/fna/PF00466.21.fna.gz    FASTA   DNA      7,026    6,653,484      135      947    2,157
    ar53_marker_genes_all_r214/fna/TIGR00064.fna.gz     FASTA   DNA      6,839    7,802,559      546  1,140.9    2,958
    ar53_marker_genes_all_r214/fna/TIGR00967.fna.gz     FASTA   DNA      6,082    9,106,488    1,077  1,497.3    2,055
    ar53_marker_genes_all_r214/fna/TIGR01171.fna.gz     FASTA   DNA        950      687,681      435    723.9      828
    bac120_marker_genes_all_r214/fna/PF00410.20.fna.gz  FASTA   DNA    368,888  146,169,498       63    396.2    2,703
    bac120_marker_genes_all_r214/fna/PF00466.21.fna.gz  FASTA   DNA    370,108  190,163,172       96    513.8    3,417
    bac120_marker_genes_all_r214/fna/TIGR00064.fna.gz   FASTA   DNA    371,270  451,371,021      471  1,215.7    4,284
    bac120_marker_genes_all_r214/fna/TIGR00967.fna.gz   FASTA   DNA    368,369  485,318,646      597  1,317.5    2,721
    bac120_marker_genes_all_r214/fna/TIGR01171.fna.gz   FASTA   DNA    367,287  303,915,555      300    827.5    1,845

So we choose TIGR00967.

## Query sequence

    # using the same genome_info.tsv file

    gene=TIGR00967
    outdir=gene_$gene

    mkdir -p $outdir
    csvtk cut -t -f 1 -U genome_info.tsv \
        | rush --eta -v gene=$gene -v outdir=$outdir \
            'seqkit grep -r -p {} bac120_marker_genes_all_r214/fna/{gene}.fna.gz \
                | seqkit replace -p "^...(.+)$" -r "\$1" -o {outdir}/{}.fasta'

    seqkit stats $outdir/*.fasta
    file                                  format  type  num_seqs  sum_len  min_len  avg_len  max_len
    gene_TIGR00967/GCF_000005845.2.fasta  FASTA   DNA          1    1,332    1,332    1,332    1,332
    gene_TIGR00967/GCF_000006765.1.fasta  FASTA   DNA          1    1,329    1,329    1,329    1,329
    gene_TIGR00967/GCF_000006945.2.fasta  FASTA   DNA          1    1,332    1,332    1,332    1,332
    gene_TIGR00967/GCF_000013425.1.fasta  FASTA   DNA          1    1,293    1,293    1,293    1,293
    gene_TIGR00967/GCF_000195955.2.fasta  FASTA   DNA          1    1,326    1,326    1,326    1,326
    gene_TIGR00967/GCF_000240185.1.fasta  FASTA   DNA          1    1,332    1,332    1,332    1,332
    gene_TIGR00967/GCF_001457635.1.fasta  FASTA   DNA          1    1,311    1,311    1,311    1,311
    gene_TIGR00967/GCF_008632635.1.fasta  FASTA   DNA          1    1,353    1,353    1,353    1,353
    gene_TIGR00967/GCF_018885085.1.fasta  FASTA   DNA          1    1,269    1,269    1,269    1,269
    gene_TIGR00967/GCF_022869705.1.fasta  FASTA   DNA          1    1,299    1,299    1,299    1,299

Taxonomy data

    csvtk cut -t -f 1 -U genome_info.tsv -o query.ass.txt

## Gold-standard truth

Query assemblies and the species

    csvtk grep -Ht -f 1 -P query.ass.txt ass2taxid2species.tsv \
        | csvtk cut -Ht -f 1,3 -o query.ass2species.tsv

    head -n 2 query.ass2species.tsv
    GCF_000005845.2 Escherichia coli
    GCF_000240185.1 Klebsiella pneumoniae

Indexed genomes (402,534)

    lexicmap utils genomes -d gtdb_complete.lmi/ -o db.ass.txt

    head -n 2 db.ass.txt
    GCA_000006155.2
    GCA_000007325.1

    csvtk grep -Ht -P db.ass.txt ass2taxid2species.tsv -o db.ass2taxid2species.tsv

    csvtk cut -Ht -f 1,3 db.ass2taxid2species.tsv -o db.ass2species.tsv


Query assembly -> species -> other assemblies of the same species.

    mkdir -p gd

    cat query.ass.txt | while read ass; do
        species="$(csvtk grep -Ht -f 1 -p $ass query.ass2species.tsv | csvtk cut -t -f 2)"
        csvtk grep -Ht -f 3 -p "$species" db.ass2taxid2species.tsv -o gd/$ass.txt
    done

Expected number of assemblies for each query.

    ls gd/*.txt | rush -k 'echo -e "{%.}\t$(csvtk nrow -H {})"' -o gd.num_ass.tsv

    $ csvtk -Ht pretty gd.num_ass.tsv
    GCF_000005845.2   33771
    GCF_000006765.1   7027
    GCF_000006945.2   13829
    GCF_000013425.1   14959
    GCF_000195955.2   7132
    GCF_000240185.1   14962
    GCF_001457635.1   8895
    GCF_008632635.1   6874
    GCF_018885085.1   2701
    GCF_022869705.1   2314

## LexicMap

    # search
    ls gene_TIGR00967/*.fasta | grep -v lexicmap | while read f; do
        lexicmap search -d gtdb_complete.lmi/ $f -o $f.lexicmap.tsv.gz
    done

    # filter and add species
    ls gene_TIGR00967/*.lexicmap.tsv.gz \
        | rush 'csvtk filter2 -t -f "\$qcovHSP >= 95" {} \
                | csvtk mutate -t -n species -f sgenome \
                | csvtk replace -t -f species -p "(.+)" -r "{kv}" -k db.ass2species.tsv -o {.}.with_species.tsv.gz'

    # species in results
    ls gene_TIGR00967/*.lexicmap.with_species.tsv.gz \
        | rush 'csvtk uniq -t -f sgenome {} | csvtk freq -t -f species -nr -o {}.freq'

    # stats
    ls gene_TIGR00967/*.fasta \
        | rush -k -v 'ass={%.}' -v 'f={}.lexicmap.with_species.tsv.gz' \
            'hits=$(csvtk head -n 1 {f} | csvtk cut -t -f hits -U); \
            species="$(csvtk grep -Ht -p {ass} db.ass2species.tsv | csvtk cut -t -f 2)";
            tp=$(csvtk uniq -t -f sgenome {f} | csvtk grep -t -f species -p "$species" | csvtk nrow -t);
            expected=$(csvtk grep -Ht -p {ass} gd.num_ass.tsv | csvtk cut -t -f 2);
            echo -e "{ass}\t$species\t$expected\t$hits\t$tp" ' \
        | csvtk add-header -Ht -n query,species,expected,hits,TP \
        | csvtk mutate2 -t -n recall -e '$TP/$expected*100' -w 4 \
            -o gene_TIGR00967/result.lexicmap.tsv

    csvtk pretty -t gene_TIGR00967/result.lexicmap.tsv

## Blastn

Other data

    # mmaping sequence id to assembly
    time cat files.txt \
        | rush --eta 'seqkit seq -ni {} | sed "s/$/\t{%..}/"' \
        | pigz -c > sseqid2ass.tsv.gz

Steps

    # concatenate all queries
    seqkit seq gene_TIGR00967/*.fasta -o gene_TIGR00967/all.fa

    # run blastn with outfmt 7
    blastdb=blastdb/gtdb
    query=gene_TIGR00967/all.fa

    # elapsed time: 1h:25m:32s
    # peak rss: 363.87 GB
    memusg -t -s "blastn -num_threads 48 -max_target_seqs 10000000 -db $blastdb -query $query \
              -outfmt '7 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand qlen qcovs qcovhsp' \
            | csvtk add-header -t -n qseqid,sseqid,pident,length,mismatch,gaps,qstart,qend,sstart,send,evalue,bitscore,sstrand,qlen,qcovs,qcovhsp \
            | gzip -c > $query.blastn.tsv.gz"

    # preprae a subseq of sseqid2ass.tsv.gz
    csvtk grep -t -P <(csvtk cut -tU -f sseqid gene_TIGR00967/all.fa.blastn.tsv.gz) \
        sseqid2ass.tsv.gz -o gene_TIGR00967/all.fa.blastn.tsv.gz.sseqid2ass.tsv.gz

    # filter and add species
    zcat gene_TIGR00967/all.fa.blastn.tsv.gz \
        | csvtk filter2 -t -f '$qcovhsp >= 95' \
        | csvtk mutate -t -n sgenome -f sseqid \
        | csvtk replace -t -f sgenome -p '(.+)' -r '{kv}' -k gene_TIGR00967/all.fa.blastn.tsv.gz.sseqid2ass.tsv.gz \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -p '(.+)' -r '{kv}' -k db.ass2species.tsv \
            -o gene_TIGR00967/all.fa.blastn.with_species.tsv.gz


    ls gene_TIGR00967/*.fasta \
        | rush --eta -k -v 'ass={%.}' -v 'f={/}/all.fa.blastn.with_species.tsv.gz' \
            'hits=$(csvtk grep -t -f qseqid -p {ass} {f} | csvtk uniq -t -f sgenome | csvtk nrow -t); \
            species="$(csvtk grep -Ht -p {ass} db.ass2species.tsv | csvtk cut -t -f 2)";
            tp=$(csvtk grep -t -f qseqid -p {ass} {f} | csvtk uniq -t -f sgenome | csvtk grep -t -f species -p "$species" | csvtk nrow -t);
            expected=$(csvtk grep -Ht -p {ass} gd.num_ass.tsv | csvtk cut -t -f 2);
            echo -e "{ass}\t$species\t$expected\t$hits\t$tp" ' \
        | csvtk add-header -Ht -n query,species,expected,hits,TP \
        | csvtk mutate2 -t -n recall -e '$TP/$expected*100' -w 4 \
            -o gene_TIGR00967/result.blastn.tsv

    csvtk pretty -t gene_TIGR00967/result.blastn.tsv

Compare

    csvtk join -t -f 1,2,3 gene_TIGR00967/result.blastn.tsv gene_TIGR00967/result.lexicmap.tsv \
        | csvtk rename2 -t -f 4-6 -p ^ -r blastn_ \
        | csvtk rename2 -t -f 7-9 -p ^ -r lexicmap_ \
        | tee gene_TIGR00967/result.tsv \
        | csvtk pretty -t
    query             species                      expected   blastn_hits   blastn_TP   blastn_recall   lexicmap_hits   lexicmap_TP   lexicmap_recall
    ---------------   --------------------------   --------   -----------   ---------   -------------   -------------   -----------   ---------------
    GCF_000005845.2   Escherichia coli             33771      83102         33631       99.5854         73803           33633         99.5914
    GCF_000006765.1   Pseudomonas aeruginosa       7027       14598         6972        99.2173         7220            6972          99.2173
    GCF_000006945.2   Salmonella enterica          13829      82567         13769       99.5661         74874           13770         99.5734
    GCF_000013425.1   Staphylococcus aureus        14959      19153         14918       99.7259         15231           14918         99.7259
    GCF_000195955.2   Mycobacterium tuberculosis   7132       10997         7120        99.8317         7146            7120          99.8317
    GCF_000240185.1   Klebsiella pneumoniae        14962      81001         14855       99.2849         73724           14854         99.2782
    GCF_001457635.1   Streptococcus pneumoniae     8895       9815          8893        99.9775         9304            8889          99.9325
    GCF_008632635.1   Acinetobacter baumannii      6874       8962          6852        99.6800         7711            6848          99.6218
    GCF_018885085.1   Clostridioides difficile     2701       2868          2665        98.6672         2725            2665          98.6672
    GCF_022869705.1   Enterococcus faecalis        2314       6100          2311        99.8704         2313            2312          99.9136

    csvtk corr -t -f blastn_recall,lexicmap_recall gene_TIGR00967/result.tsv
    blastn_recall   lexicmap_recall 0.9974

Plot

    cat gene_TIGR00967/result.tsv \
        | csvtk cut -t -f query,blastn_recall,lexicmap_recall \
        | csvtk rename  -t -f blastn_recall,lexicmap_recall -n Blastn,LexicMap \
        > marker_gene.tsv

    Rscript plot.marker_gene.R
