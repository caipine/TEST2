include: "rules/common_A.smk"


############################



ALL_FASTQ =  expand("trimmed/{sample1}.fastq.gz", sample1 = ALL_SAMPLES)
print("****")
print(ALL_FASTQ)

TARGETS = []
TARGETS.extend(ALL_FASTQ)



localrules: all
rule all:
     input: TARGETS


def get_fastqAAA(wildcards):
    #sample = "-".join(wildcards.sample.split("-")[0:-1])
    #unit = wildcards.sample.split("-")[-1]
    #return units.loc[(sample, unit), ["fq1", "fq2"]].dropna()
    return units.loc[("-".join(wildcards.sample.split("-")[0:-1]), wildcards.sample.split("-")[-1]), ["fq1", "fq2"]].dropna()

#print("-".join(wildcards.sample1.split("-")[0:-1]))

######################
rule trim_reads_se:
    input:
        get_fastqAAA
        #to get SRA/fastqs/SRR1374316_1.fastq.gz, SRA/fastqs/SRR1374316_2.fastq.gz
    output:
        fastq=temp("trimmed/{sample}.fastq.gz"),
        qc="trimmed/{sample}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt"]["se"])
    log:
        "logs/cutadapt/{sample}.log"
    wrapper:
        "0.17.4/bio/cutadapt/se"

######################
rule trim_reads_pe:
    input:
        get_fastq
    output:
        fastq1=temp("trimmed/{sample}-{unit}.1.fastq.gz"),
        fastq2=temp("trimmed/{sample}-{unit}.2.fastq.gz"),
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt"]["pe"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"



