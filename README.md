# Practical course in pan-genomics, transcriptomics and functional genomics

**Martin Mascher · IPK Gatersleben · 2020**

## Contents

- [Software installation](#software-installation)
- [Comparing two rice genome assemblies](#comparing-two-rice-genome-assemblies)
  - [Downloading the assemblies](#downloading-the-assemblies)
  - [Chromosome-wise alignment](#chromosome-wise-alignment)
  - [Comparing a single gene in two assemblies](#comparing-a-single-gene-in-two-assemblies)
- [Genebank genomics portal](#genebank-genomics-portal)
- [Analysis of RNA-seq of a wheat mutant and its wildtype](#analysis-of-rna-seq-of-a-wheat-mutant-and-its-wildtype)
  - [Data download](#data-download)
    - [RNAseq reads](#rnaseq-reads)
    - [Wheat annotation](#wheat-annotation)
  - [Creating a kallisto index](#creating-a-kallisto-index)
  - [Quantifying transcript abundance](#quantifying-transcript-abundance)
- [Transcriptome atlases](#transcriptome-atlases)
- [Barley epigenome browser](#barley-epigenome-browser)
- [Cereal pan-genomes at grain genes](#cereal-pan-genomes-at-grain-genes)

---

## Software installation

1. All the commands were tested on an Ubuntu virtual machine (VM) with standard software configuration (Ubuntu 20.04.1 LTS). The VM has 8 GB of main memory and runs on an Intel(R) Xeon(R) Gold 5218 CPU (2.30GHz, 2 physical cores), and has access to 500 GB of hard disk space. Fewer than 8 GB of RAM would prevent one programs we will use (minimap2) to run properly, but most steps should work with 4 GB or even 2 GB RAM. A total of 70 GB of disk space are needed to store the RNA-seq dataset we will analyze.

2. Information on how to setup an Ubuntu VM were given in the [course](https://github.com/eead-csic-compbio/scripting_linux_shell) of Bruno Contreras and Carlos Cantalapiedra. See their [session 0](https://github.com/eead-csic-compbio/scripting_linux_shell/blob/master/session0.md).

3. Install the required tools with the [apt](http://manpages.ubuntu.com/manpages/bionic/man8/apt.8.html) package manager.

    You will be prompted to enter your user password.
    ```bash
    sudo apt install samtools kallisto minimap2 python git r-base ncbi-blast+ libssl-dev curl libcurl4-openssl-dev pandoc emboss mafft snp-sites gmap
    ```

    - [samtools](http://samtools.github.io), [minimap2](https://github.com/lh3/minimap2), [ncbi-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [emboss](http://emboss.sourceforge.net), [mafft](https://mafft.cbrc.jp/alignment/software/), [snp-sites](http://sanger-pathogens.github.io/snp-sites/), [gmap](http://research-pub.gene.com/gmap/) and [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) will be used in the analysis of two rice genomes.

    - [kallisto](https://pachterlab.github.io/kallisto/) will be used to analyse RNAseq data of wheat.

    - [python](https://en.wikipedia.org/wiki/Python_(programming_language)), [git](https://en.wikipedia.org/wiki/Git), [libcurl4-openssl-dev](https://packages.ubuntu.com/xenial/libcurl4-openssl-dev), [curl](https://curl.se), [libssl-dev](https://packages.debian.org/jessie/libssl-dev), and [pandoc](https://pandoc.org) are dependencies of dotPlotly, a program to show alignments of two genome assemblies.

    - r-base is the Ubuntu package of the [R statistical environment](https://www.r-project.org).

4. Now install the required R packages. These are also dependencies of [dotPlotly](https://github.com/tpoorten/dotPlotly).

5. Start R by running `R` on the command line. You will be prompted to enter your password.

    ```r
    sudo R
    install.packages('optparse')
    install.packages('ggplot2')
    install.packages('plotly')
    ```

## Comparing two rice genome assemblies

### Downloading the assemblies

1. Create a folder for downloads and change to it.

    ```bash
    mkdir -p ~/international_masters/downloads
    cd ~/international_masters/downloads
    ```

2. Download the IR64 assembly from the [ROOTomics website](https://rootomics.dna.affrc.go.jp/en/research/IR64). It is described in this [paper](https://www.g3journal.org/content/10/5/1495).

    ```bash
    wget https://rootomics.dna.affrc.go.jp/genomes/1/dl
    ```

3. Take a look at the file, rename, and unzip it.

    ```bash
    less dl
    mv dl IR64.fasta.gz
    gunzip IR64.fasta.gz
    ```

4. Download the Nipponbare assembly from the [Rice Genome Annotation Project website](http://rice.plantbiology.msu.edu). It is described in [this paper](https://thericejournal.springeropen.com/articles/10.1186/1939-8433-6-4).

    ```bash
    wget http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.chrs.con
    ```

5. Take a look at the file and rename it.

    ```bash
    less all.chrs.con
    mv all.chrs.con Nipponbare.fasta
    ```

### Chromosome-wise alignment

1. Create the project directory and change into it.

    ```bash
    mkdir ~/international_masters/rice_assemblies/
    cd ~/international_masters/rice_assemblies/
    ```

2. Create symbolic links to the assembly FASTA files in the downloads folder.

    ```bash
    ln -s ~/international_masters/downloads/IR64.fasta
    ln -s ~/international_masters/downloads/Nipponbare.fasta
    ```

3. Create index files recording the chromosome lengths in each assembly.

    ```bash
    samtools faidx IR64.fasta
    samtools faidx Nipponbare.fasta
    ```

4. Align the two assemblies with [minimap2](https://github.com/lh3/minimap2).

    ```bash
    minimap2 -I500M -x asm5 Nipponbare.fasta IR64.fasta > Nipponbare_IR64.paf 2> Nipponbare_IR64.paf.err
    ```

    > **Notes**
    >
    > - The parameter preset `-x asm5` is suited to align assemblies from the same species.
    > - The value of `-I` is reduced from the default (4 GB) to 500 MB to reduce the memory requirments.
    > - `>` and `2>` are used to direct the standard output (stdout) and standard error (stderr).

5. Send the `minimap2` process to the background with Ctrl-Z. Resume it with `bg`. Check with `top` that it's running.

6. In the meantime, take a look the definition of the PAF format at the end of the minimap2 manpage.

    ```bash
    man minimap2
    ```

7. Once the `minimap2` is done, take a look at output.

    ```bash
    cat Nipponbare_IR64.paf
    head Nipponbare_IR64.paf
    less Nipponbare_IR64.paf
    less -S Nipponbare_IR64.paf
    column -t Nipponbare_IR64.paf | less -S
    ```

8. Install [dotPlotly](https://github.com/tpoorten/dotPlotly) from Github and check that it works.

    ```bash
    git clone https://github.com/tpoorten/dotPlotly.git
    ./dotPlotly/pafCoordsDotPlotly.R -h
    ```

9. Create a dotplot from PAF alignment file.

    ```bash
    ./dotPlotly/pafCoordsDotPlotly.R -i Nipponbare_IR64_sorted.paf -o Nipponbare_IR64 -m 1000 -l
    ```

10. Copy the files to your main system and open them.

### Comparing a single gene in two assemblies

1. Download the genomic sequence of [PLASTOCHRON1](https://dx.doi.org/10.1073%2Fpnas.2636936100) gene from NCBI accession [AB096259.1](https://www.ncbi.nlm.nih.gov/nuccore/AB096259), rename the file `pla1.fasta`  and copy it to the working directory `~/international_masters/rice_assemblies`.

2. Take a look at the sequence.

    ```bash
    less pla1.fasta
    ```

3. Create BLAST database for the two genome assemblies.

    ```bash
    makeblastdb -dbtype nucl -in Nipponbare.fasta
    makeblastdb -dbtype nucl -in IR64.fasta
    ```

4. Run the BLAST alignment and output to [tabular format](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6).

    ```bash
    blastn -db Nipponbare.fasta -query pla1.fasta -outfmt 6 > pla1_Nipponbare.txt
    blastn -db IR64.fasta -query pla1.fasta -outfmt 6 > pla1_IR64.txt
    ```

    > **Notes**
    >
    > - `-outfmt 6` means [tabular](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) output.

5. Compare the results to BLAST web tool provided on the [ROOTomics website](https://rootomics.dna.affrc.go.jp/en/research/IR64).

6. Check the BLAST version.

    ```bash
    blastn -version
    ```

7. Extract the aligned sequence of the first exon with `samtools faidx` [[man page](http://www.htslib.org/doc/samtools-faidx.html)].

    ```bash
    samtools faidx IR64.fasta chr10:10413299-10414334 > pla1_IR64.fasta
    samtools faidx Nipponbare.fasta Chr10:13659508-13660543 > pla1_Nipponbare.fasta
    ```

8. Extract the sequence of all BLAST hits.

    ```bash
    cat pla1_IR64.txt | awk '$9 < $10 {print $2":"$9"-"$10} $10 < $9 {print $2":"$10"-"$9}' \
     | sort | xargs samtools faidx  IR64.fasta  > pla1_IR64_all_hits.fasta
    ```

    > **Notes**
    >
    > - The backslash `\` makes it possible to split long lines into two.

9. Align the first exon of PLA1 sequence of Nipponbare and IR64 using [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/).

10. Align the two sequence using MAFFT and find SNPs between them with SNP-sites:

    ```bash
    cat pla1_IR64.fasta pla1_Nipponbare.fasta | mafft - > pla1_mafft.fasta
    snp-sites -v pla1_mafft.fasta
    snp-sites pla1_mafft.fasta
    ```

11. Now we use [GMAP](https://academic.oup.com/bioinformatics/article/21/9/1859/409207) for spliced alignment

    to extract and compare alignments of of the full transcript, not only the first exon.

12. Build the GMAP indices for both genomes.

    ```bash
    gmap_build Nipponbare.fasta -D . -d Nipponbare_db > Nipponbare_build.out 2> Nipponbare_build.err &
    gmap_build IR64.fasta -D . -d IR64_db > IR64_build.out 2> IR64_build.err &
    ```

    > **Notes**
    >
    > - The `&` at the end of the line sends the process immediately to the background.

13. Align the PLA1 sequence to both genomes.

    ```bash
    gmap -f 2 -D . -d Nipponbare_db pla1.fasta > pla1_Nipponbare_gmap.gff
    gmap -f 2 -D . -d IR64_db pla1.fasta > pla1_IR64_gmap.gff
    ```

    > **Notes**
    >
    > - `-f 2` generate GFF output. A description of the GFF format can be found [here](https://www.ensembl.org/info/website/upload/gff.html).

14. Extract the matched sequence with [gffread](https://github.com/gpertea/gffread) (part of [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)).

    ```bash
    cat pla1_Nipponbare_gmap.gff | gffread -g Nipponbare.fasta -w pla1_Nipponbare_gmap_mRNA.fasta
    cat pla1_IR64_gmap.gff | gffread -g IR64.fasta -w pla1_IR64_gmap_mRNA.fasta
    ```

15. Run the multiple sequence aligment and SNP extraction with the full transcript sequences.

    ```bash
    cat pla1_Nipponbare_gmap_mRNA.fasta pla1_IR64_gmap_mRNA.fasta | mafft - > pla1_gmap_mafft.fasta
    snp-sites -v pla1_gmap_mafft.fasta
    ```

## Genebank genomics portal

1. Visit [BRIDGE](https://bridge.ipk-gatersleben.de/#geomap), the barley genebank genomics portal.

## Analysis of RNA-seq of a wheat mutant and its wildtype

### Data download

#### RNAseq reads

1. We will download RNA-seq for 16 bread wheat samples reported by [Sakuma et al. 2019.](https://www.pnas.org/content/116/11/5182)

2. In a web browser, go to https://www.ebi.ac.uk/ena/browser/view/PRJEB25119.

3. Click on "Show Column Selection".

4. Enable "submitted_md5" and "submitted_ftp"

5. Right-click on the "TSV" link next to "Download report"  and select "Copy link" to copy the download URL to the clipboard.

6. Go to the UNIX terminal and run the following commands.

    ```bash
    cd ~/international_masters/downloads
    wget -O file_list.tsv 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB25119&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_md5,submitted_ftp,sra_ftp&format=tsv&download=true'
    ```

    > **Notes**
    >
    > - The long URL is one you copied to clipboard. You can paste it with Ctrl+Shift+V. The download file will be named `file_list.tsv`.

7. View the columns the file.

    ```bash
    column -t file_list.tsv | less -S
    head file_list.tsv -n 1 | tr '\t' '\n' | nl
    ```

8. Extract the FTP links from the TSV file.

    ```bash
    tail -n +2 file_list.tsv  | cut -f 9 | tr ';' '\n' > ftp_links.txt
    ```

9. Download all the files using `wget`.

    ```bash
    cat ftp_links.txt | while read i; do basename $i | sed 's/.fastq.gz//' | read b; wget -nv $i > $b.out > $b.err; done
    cat ftp_links.txt | while read i; do b=$(basename $i | sed 's/.fastq.gz//'); wget -nv $i > $b.out > $b.err; done
    ```

10. Generate MD5 sums for the download files.

    ```bash
    md5sum GA*gz WA*gz > md5sum.txt
    ```

11. Compare the MD5 check sums reported in the TSV file downloaded from ENA to the ones you have just calcuated.

    ```bash
    cat file_list.tsv | cut -f 8 | tr ';' '\n'  | tail -n +2 | sort | md5sum
    #f0e418bc72ca6433012e0182e54348a4  -

    cat md5sum.txt  | cut -d ' ' -f 1 | sort | md5sum
    #f0e418bc72ca6433012e0182e54348a4  -
    ```

    > **Notes**
    >
    > - The two need to match. Otherwise, there was an error during the transfer and some files have to be downloaded again.

#### Wheat annotation

1. Download the coding sequences of the genes annotated on the wheat reference genones from [URGI](https://wheat-urgi.versailles.inra.fr/Seq-Repository/Annotations).

    ```bash
    cd international_masters/downloads
    wget 'https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa.zip'
    ```

2. Download the functional annotation of the gene models.

    ```bash
    wget 'https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_FunctionalAnnotation_v1.zip'
    ```

### Creating a kallisto index

1. Create the project directory and change to it.

    ```bash
    mkdir ~/international_masters/wheat_rnaseq
    cd ~/international_masters/wheat_rnaseq
    ```

2. Create symbolic links to the annotation files and decompress them.

    ```bash
    ln -s ~/international_masters/downloads/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa.zip .
    ln -s ~/international_masters/downloads/iwgsc_refseqv1.0_FunctionalAnnotation_v1.zip .

    unzip iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa.zip
    unzip iwgsc_refseqv1.0_FunctionalAnnotation_v1.zip
    ```

3. Count the number of sequences in the file.

    ```bash
    grep -c '^>' iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa
    ```

4. Create an index for alignment with [kallisto](https://pachterlab.github.io/kallisto/).

    ```bash
    kallisto index --index wheat_index iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa > kallisto_index.out 2>  kallisto_index.err &
    ```

### Quantifying transcript abundance

1. Create symbolic links to read files.

    ```bash
    ln -s ~/international_masters/downloads/GA*gz .
    ln -s ~/international_masters/downloads/WA*gz .
    ```

2. Run the quantification for a single sample.

    ```bash
    kallisto quant --index wheat_index GA_0908-N_1_R1.fastq.gz GA_0908-N_1_R2.fastq.gz \
     --output GA_0908-N_1_kallisto > GA_0908-N_1_kallisto.out 2> GA_0908-N_1_kallisto.err &
    ```

3. Run the quantification for all samples using a loop.

    ```bash
     find | grep R1 | cut -d _ -f -3 | sort | while read i; do
      kallisto quant --index wheat_index  ${i}_R1.fastq.gz ${i}_R2.fastq.gz --output ${i}_kallisto > ${i}_kallisto.out 2>  ${i}_kallisto.err
     done
    ```

4. Check that there are results for samples.

    ```bash
    head GA_0908-N_1_kallisto/abundance.tsv | column -t
    find -L | grep abundance.tsv | wc  -l
    find -L | grep abundance.tsv | xargs wc -l
    grep -c '^>' iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa
    ```

5. Further analyses will be run [3D RNA-seq](https://github.com/wyguo/ThreeDRNAseq) ([paper](https://www.biorxiv.org/content/10.1101/656686v1)).

6. Create tables with the metadata.

    ```bash
    grep '^>' iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa | tr -d '>' | cut -d ' ' -f 1 > transcripts.txt
    cut -d . -f 1 transcripts.txt > genes.txt
    paste -d , transcripts.txt genes.txt  > transcript_genes.csv

    find -L -type d | grep kallisto  | cut -d / -f 2 | sort > kallisto.txt

    cat kallisto.txt  | tr _- '\t' | awk '{print $1","$3","$4}' \
     | paste -d , - kallisto.txt | awk 'BEGIN{print "stage,allele,rep,folder"} {print}' > sample_info.csv
    ```

    > **Notes**
    >
    > - Assignment of transcript isoforms to genes.
    > - List of Kallisto output directories.
    > - Table with experimental design.

7. Copy the Kallisto output folders for all samples to your local machine (Mac/Windows) and create a ZIP archive containing all output folders.

8. Open the [3D RNA-seq app](https://3drnaseq.hutton.ac.uk/app_direct/3DRNAseq/).

9. Upload the datasets in the *Data generation* tab.

10. Follow the 3D DNA-seq steps. Click on the question mark symbol to get more guidance.

## Transcriptome atlases

1. Visit [BAR](http://bar.utoronto.ca), the Bio-Analytic Resource for Plant Biology ([paper](https://link.springer.com/protocol/10.1007%2F978-1-4939-6658-5_6)).

## Barley epigenome browser

1. Visit the [Barley epigenome browser](https://ics.hutton.ac.uk/barley-epigenome/).

## Cereal pan-genomes at grain genes

1. Cereal pan-genomes are hosted at [GrainGenes](https://wheat.pw.usda.gov/GG3/).
