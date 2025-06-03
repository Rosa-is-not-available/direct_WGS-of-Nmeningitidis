# direct_WGS-of-Nmeningitidis
This repository contains code used to analyse data from the publication entitled "Direct whole-genome sequencing enables strain typing of unculturable Neisseria meningitidis from oropharyngeal carriage specimens". This data analysis was designed for use with short-read libraries prepared using the Agilent SureSelect XT target enrichment system.

## basecalling settings used for processing of dWGS Neisseria meningitidis carriage specimen illumina sequencing reads
```
bcl2fastq --sample-sheet {illumina_samplesheet} --runfolder-dir {sequencing_run_directory} --mask-short-adapter-reads 0 --use-bases-mask Y150,I8,N10,Y150 --no-lane-splitting --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --output-dir {sequencing_run_directory}/BaseCalls
```

## snakemake rules used for processing of dWGS Neisseria meningitidis carriage specimen illumina sequencing reads
### prepare reads for assembly ###
```
rule trim_reads:
     input: 
          r1 = "input/{sample}_R1.fastq.gz",
          r2 = "input/{sample}_R2.fastq.gz"
     output:
          t1 = "input/{sample}_T1.fastq.gz",
          t2 = "input/{sample}_T2.fastq.gz",
          up1 = "input/unpaired/{sample}_UP1.fastq.gz",
          up2 = "input/unpaired/{sample}_UP2.fastq.gz"
     conda: "phesiqcal"
     threads: 4
     params:
          ILLUMINACLIP = "ILLUMINACLIP:/path/to/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10"
     shell: 
          "trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.t1} {output.up1} {output.t2} {output.up2} {params.ILLUMINACLIP} SLIDINGWINDOW:4:20 MINLEN:30"

rule human_removal:
     input:
          t1 = rules.trim_reads.output.t1,
          t2 = rules.trim_reads.output.t2
     output: 
          f1 = "filtered_input/{sample}_F1.fastq.gz",
          f2 = "filtered_input/{sample}_F2.fastq.gz"
     conda: "metaphe"
     threads: 4
     params:
          REF = "/path/to/GRCh38.p14/GRCh38.p14"
     shell:
          "bowtie2 --very-sensitive-local -p {threads} --seed 1000 -x {params.REF} -1 {input.t1} -2 {input.t2} | samtools fastq -1 {output.f1} -2 {output.f2} -f 12 -F 256"
```

### perform denovo assembly and filter contigs based on size ###
```
rule shovill:
     input: 
          f1 = rules.human_removal.output.f1,
          f2 = rules.human_removal.output.f2
     output:
          directory("{sample}/shovill")
     conda: "phesiqcal"
     threads: 4
     shell:
          "shovill --ram 32 --force --depth 80 --gsize 2M --cpus {threads} -R1 {input.f1} -R2 {input.f2} --outdir {output} "

rule remove_small:
     input:
          rules.shovill.output
     output:
          "{sample}/{sample}_raw_contigs.fna"
     conda: "phesiqcal"
     shell:
          "seqtk seq -L 1000 {input}/contigs.fa > {output}"
```
### assess QC of reads ###
```
rule fq_input:
     input: 
          r1 = "input/{sample}_R1.fastq.gz",
          r2 = "input/{sample}_R2.fastq.gz"
     output:
          "{sample}/yield_input.tab"
     conda: "phesiqcal"
     shell:
          "fq {input.r1} {input.r2} > {output}"

rule fq_trim:
     input: 
          t1 = rules.trim_reads.output.t1,
          t2 = rules.trim_reads.output.t2
     output:
          "{sample}/yield_trim.tab"
     conda: "phesiqcal"
     shell:
          "fq {input.t1} {input.t2} > {output}"

rule fq_filtered:
     input:
          f1 = rules.human_removal.output.f1,
          f2 = rules.human_removal.output.f2
     output:
          "{sample}/yield_filtered.tab"
     conda: "phesiqcal"
     shell:
          "fq {input.f1} {input.f2} > {output}"
```
### assess QC species identificantion ###
```
rule kraken_trim:
     input: 
          t1 = rules.trim_reads.output.t1,
          t2 = rules.trim_reads.output.t2
     output:
          "{sample}/kraken2_trim.tab"
     conda: "phesiqcal"
     threads: 4
     params:
          DB = "/path/to/k2_pluspf_20220607/"
     shell:
          "kraken2 --threads {threads} --db {params.DB} --memory-mapping --report {output} --paired {input.t1} {input.t2} "

rule kraken_filtered:
     input:
          f1 = rules.human_removal.output.f1,
          f2 = rules.human_removal.output.f2
     output:
          "{sample}/kraken2_filtered.tab"
     conda: "phesiqcal"
     threads: 4
     params:
          DB = "/path/to/k2_pluspf_20220607/"
     shell:
          "kraken2 --threads {threads} --db {params.DB} --memory-mapping --report {output} --paired {input.f1} {input.f2}"
```
### assess QC of denovo assemblies ### 
```
rule QC_denovo:
   input:
        expand("{sample}/{sample}_raw_contigs.fna", sample=config["samples"])
   output:
        "denovo_raw.tab"
   conda: "phesiqcal"
   shell:
         "fa -e -t {input} >> {output}"
```
### filter contigs that map to N.meningitidis reference ###
```
rule filter_Nmen_contigs:
     input: 
          rules.remove_small.output
     output:
          "filtered_contigs/align_bam/{sample}_lra.bam" 
     conda: "phedir"
     threads: 4
     params: 
          REF = "{/path/to/reference/}ncbi/GCF_008330805.1_ASM833080v1_genomic.fna" 
     shell:
          "lra align -CONTIG {params.REF} {input} -t {threads} -p s | samtools view -F 4 -b | samtools view -F 2048 -b >  {output}"

rule convert_fasta:
     input:
          rules.filter_Nmen_contigs.output
     output: 
          "filtered_contigs/{sample}_lra.fna"
     conda: "phedir"
     shell:
          "samtools fasta {input} > {output}"
```
### assess QC of Nmen filtered denovo assemblies ### 
```
rule QC_Nmen_denovo:
   input:
        expand("filtered_contigs/{sample}_lra.fna", sample=config["samples"])
   output:
        "denovo_filtered.tab"
   conda: "phesiqcal"
   shell:
         "fa -e -t {input} >> {output}"
```
### perform typing on Nmen filtered denovo assemblies ### 
```
rule mlst:
   input:
        expand("filtered_contigs/{sample}_lra.fna", sample=config["samples"])
   output:
        "mlst.tab"
   conda: "phesiqcal"
   shell:
        "mlst {input} >> {output}"

rule meningotype:
     input: 
          expand("filtered_contigs/{sample}_lra.fna", sample=config["samples"])
     output: 
          "meningotype.tab"
     conda: "phetype"
     shell:
          "meningotype --porB --finetype --bast {input} >> {output}"
```

## Citations
Acknowledgments to all the authors of tools used in the pipeline.

Trimmomatic:
https://github.com/usadellab/Trimmomatic

bowtie2:
Langmead B, Salzberg S, https://github.com/BenLangmead/bowtie2

Nullarbor:
Seemann T, Goncalves da Silva A, Bulach DM, Schultz MB, Kwong JC, Howden BP. Nullarbor Github https://github.com/tseemann/nullarbor

kraken2:
Taxonomic sequence classifier that assigns taxonomic labels to DNA sequences Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2 Genome Biol 20, 257 (2019)

Shovill:
Seemann T, Shovill Github https://github.com/tseemann/shovill

LRA:
https://github.com/ChaissonLab/LRA

samtools:
Twelve years of SAMtools and BCFtools
Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

MLST:
Seemann T, mlst Github https://github.com/tseemann/mlst

meningotype:
Kwong J, Stroehlein A, Gonçalves da Silva A, Seemann T, https://github.com/MDU-PHL/meningotype

sankemake:
Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.




