configfile: "config.yaml"
ENVIRONMENT = "environment.yaml"

SAMPLE_IDS = config["SAMPLE_IDS"]  # sample IDs
READS_DIR = config["READS_DIR"]  # folder with reads in fastq
RESULTS_DIR = config["RESULTS_DIR"]  # path where output will be stored
DIAG_DIR = config["DIAG_DIR"]  # path to diagnostic outputs

# input data
REF_RNA = config["REF_RNA"]  # fasta with reference transcripts
ADAPTERS = config["ADAPTERS"]  # illumina and RT-PCR adapters
CONTAMINANTS = config["CONTAMINANTS"]  # radish contaminants
T2G = config["T2G"]  # transcript-to-gene correspondenceq
PHENOTABLE = config["PHENOTABLE"]
ANNOTATION = config["ANNOTATION"]
DESCRIPTION = config["DESCRIPTION"]
AT_TFS = config["AT_TFS"]  # table with annotated At TFs
# MANUAL_TFS = config["MANUAL_TFS"]  # manual filter imposed on TFs (see TF.R for details)
PLANT_GSEA = config["PLANT_GSEA"]  # gene list matrix transposed from PlantGSEA

if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

if not os.path.exists(DIAG_DIR):
    os.makedirs(DIAG_DIR)


# Params
THREADS = config["THREADS"]
MEMORY = config["MEMORY"]


rule all:
    input:
        DIAG_DIR + "/multiqc_report.html",
        DEtable=RESULTS_DIR + "/DEtable.csv",
        up_genes=RESULTS_DIR + "/up_genes.csv",
        down_genes=RESULTS_DIR + "/down_genes.csv",
        normalized_expression=RESULTS_DIR + "/normalized_expression.csv",
        GO_up=RESULTS_DIR + "/GO_up.txt",
        GO_down=RESULTS_DIR + "/GO_down.txt",
        GSEA_up=RESULTS_DIR + "/GSEA_up.txt",
        GSEA_down=RESULTS_DIR + "/GSEA_down.txt",
        GO_dotplot=RESULTS_DIR + "/GO_dotplot.png",
        KEGG_dotplot=RESULTS_DIR + "/KEGG_dotplot.png",
        GSEA_volcano=RESULTS_DIR + "/GSEA_volcano.png",
        figure=RESULTS_DIR + "/TF_figure.png"
    output: touch(".status")


rule bbduk_adapters:
    input:
        in1=READS_DIR + "/{sample}_1.fastq.gz",
        in2=READS_DIR + "/{sample}_2.fastq.gz"
    output:
        out1=temp(READS_DIR + "/{sample}_1.adaptless.fastq.gz"),
        out2=temp(READS_DIR + "/{sample}_2.adaptless.fastq.gz")
    params:
        k=23,   # kmer for finding contaminants
        mink=11,    # look for shorter kmers are read ends
        hdist=1,    # hamming distance from reference (affects memory usage)
        trimq=10,    # regions with average quality below will be trimmed
        qtrim="rl", # trim both ends
        ktrim="r"    # trim to the right to remove reference
    threads: THREADS
    conda: ENVIRONMENT
    shell:
        "bbduk.sh -Xmx{MEMORY} -t={threads}"
        "in1={input.in1} in2={input.in2} "
        "out1={output.out1} out2={output.out2} "
        "ref={ADAPTERS} ktrim={params.ktrim} qtrim={params.qtrim} "
        "trimq={params.trimq} k={params.k} mink={params.mink} "
        "hdist={params.hdist} tpe tbo"    # tpe - trim paired to minimum length, tbo - trim based on where reads overlap


rule trimmomatic:
    input:
        in1=READS_DIR + "/{sample}_1.adaptless.fastq.gz",
        in2=READS_DIR + "/{sample}_2.adaptless.fastq.gz"
    output:
        out1=temp(READS_DIR + "/{sample}_1.paired.fastq.gz"),
        out2=temp(READS_DIR + "/{sample}_2.paired.fastq.gz"),
        out1u=temp(READS_DIR + "/{sample}_1.unpaired.fastq.gz"),
        out2u=temp(READS_DIR + "/{sample}_2.unpaired.fastq.gz")
    threads: THREADS
    conda: ENVIRONMENT
    shell:
        "java -jar trimmomatic-0.39.jar PE -threads {threads} -phred33 {input.in1} {input.in2} "
        "{output.out1} {output.out1u} {output.out2} {output.out2u} "
        "ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"


rule cutadapt:
    input:
        in1=READS_DIR + "/{sample}_1.paired.fastq.gz",
        in2=READS_DIR + "/{sample}_2.paired.fastq.gz"
    output:
        out1=temp(READS_DIR + "/{sample}_1.cut.fastq.gz"),
        out2=temp(READS_DIR + "/{sample}_2.cut.fastq.gz")
    threads: THREADS
    conda: ENVIRONMENT
    shell:
        "cutadapt --cores={threads} -b GTGGTCATTACGGCCGCGGG -B GTGGTCATTACGGCCGCGGG "  # removing additional overrepresented sequences
        "-b GAGTGGCCATTACGGGGGG -B GAGTGGCCATTACGGGGGG "
        "-o {output.out1} -p {output.out2} {input.in1} {input.in2}"


rule bbduk_contaminants:
    input:
        in1=READS_DIR + "/{sample}_1.cut.fastq.gz",
        in2=READS_DIR + "/{sample}_2.cut.fastq.gz"
    output:
        out1=temp(READS_DIR + "/{sample}_1.unmatched.fastq.gz"),
        out2=temp(READS_DIR + "/{sample}_2.unmatched.fastq.gz"),
        outm1=temp(READS_DIR + "/{sample}_1.matched.fastq.gz"),
        outm2=temp(READS_DIR + "/{sample}_2.matched.fastq.gz"),
        stats=READS_DIR + "/{sample}_contamination_stats.txt",
        rpkm=READS_DIR + "/{sample}_rpkm.txt"
    threads: THREADS
    params:
        k=31,
        hdist=1
    conda: ENVIRONMENT
    shell:
        "bbduk.sh -Xmx{MEMORY} t={threads} "
        "in1={input.in1} in2={input.in2} "
        "out1={output.out1} out2={output.out2} "
        "outm1={output.outm1} outm2={output.outm2} "
        "ref={CONTAMINANTS} k={params.k} hdist={params.hdist} "
        "stats={output.stats} rpkm={output.rpkm}"


rule kallisto_idx:
    input: REF_RNA
    output: temp(RESULTS_DIR + "/kallisto.idx")
    conda: ENVIRONMENT
    shell:
        "kallisto index -i {output} {input}"


rule kallisto_quant:
    input:
        index=RESULTS_DIR + "/kallisto.idx",
        reads_1=READS_DIR + "/{sample}_1.unmatched.fastq.gz",
        reads_2=READS_DIR + "/{sample}_2.unmatched.fastq.gz"
    output:
        log=DIAG_DIR + "/{sample}.kallisto.log",
        outfile=RESULTS_DIR + "/{sample}.kallisto/abundance.h5"
    params:
        outdir=directory(RESULTS_DIR + "/{sample}.kallisto")
    threads: THREADS
    conda: ENVIRONMENT
    shell:
        "kallisto quant -i {input.index} -o {params.outdir} -b 50 -t {threads} {input.reads_1} {input.reads_2} 2> {output.log}"

rule DE:
    input:
        annotation=ANNOTATION,
        description=DESCRIPTION,
        abundance=expand(RESULTS_DIR + "/{sample}.kallisto/abundance.h5", sample = SAMPLE_IDS),
        t2g=T2G,
        phenotable=PHENOTABLE
    output:
        DEtable=RESULTS_DIR + "/DEtable.csv",
        up_genes=RESULTS_DIR + "/up_genes.csv",
        down_genes=RESULTS_DIR + "/down_genes.csv",
        normalized_expression=RESULTS_DIR + "/normalized_expression.csv"
    conda: ENVIRONMENT
    script:
        "DE.R"


rule TF_picture:
    input:
        DEtable=RESULTS_DIR + "/DEtable.csv",
        normalized_expression=RESULTS_DIR + "/normalized_expression.csv",
        AtTFs=AT_TFS,
        # manual_TFs=MANUAL_TFS,
        phenotable=PHENOTABLE
    output:
        figure=RESULTS_DIR + "/TF_figure.png"
    conda: ENVIRONMENT
    script:
        "TF.R"


rule enrichment:
    input:
        DEtable=RESULTS_DIR + "/DEtable.csv",
        PlantGSEA=PLANT_GSEA
    output:
        GO_up=RESULTS_DIR + "/GO_up.txt",
        GO_down=RESULTS_DIR + "/GO_down.txt",
        GSEA_up=RESULTS_DIR + "/GSEA_up.txt",
        GSEA_down=RESULTS_DIR + "/GSEA_down.txt",
        GO_dotplot=RESULTS_DIR + "/GO_dotplot.png",
        KEGG_dotplot=RESULTS_DIR + "/KEGG_dotplot.png",
        GSEA_volcano=RESULTS_DIR + "/GSEA_volcano.png"
    conda: ENVIRONMENT
    script:
        "enrichment.R"


rule quality_report:
    input:
        kallisto_log=expand(DIAG_DIR + "/{sample}.kallisto.log", sample = SAMPLE_IDS),
        reads_1=expand(READS_DIR + "/{sample}_1.unmatched.fastq.gz", sample = SAMPLE_IDS),
        reads_2=expand(READS_DIR + "/{sample}_2.unmatched.fastq.gz", sample = SAMPLE_IDS)
    output:
        DIAG_DIR + "/multiqc_report.html"
    conda: ENVIRONMENT
    shell:
        "fastqc {input.reads_1} -o {DIAG_DIR};"
        "fastqc {input.reads_2} -o {DIAG_DIR};"
        "multiqc {DIAG_DIR} -o {DIAG_DIR} -n {output}"
