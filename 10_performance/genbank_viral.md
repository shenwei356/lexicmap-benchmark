# Genbank viral genomes

## Data

    time genome_updater.sh -d "refseq,genbank" -g "viral" \
        -f "genomic.fna.gz" -o "genbank_viral" -M "ncbi" -a -t 12 -m -L curl -i

    # for convenience, creat a new directory containing renamed symbol links pointing to the orginal files
    mkdir files; cd files
    find ../files -name "*.fna.gz" \
        | rush --eta -v 'id={%@^(..._.........\.\d+)}' \
                     -v 'dir1={%@^(..._...)}' -v 'dir2={%@^..._...(...)}' \
            'mkdir -p {dir1}/{dir2}; cd {dir1}/{dir2}; ln -s ../../{} {id}.fna.gz'
    cd ..

    fd .fna.gz$ files/ > files.txt

## LexicMap

Indexing

    memusg -t -s "lexicmap index -S -X files.txt -b 50000 -O genbank_viral.lmi" > genbank_viral.lmi.log 2>&1

