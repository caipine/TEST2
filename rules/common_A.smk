import pandas as pd
from snakemake.utils import validate

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_table(config["ref"]["genome"] + ".fai",
                       header=None, usecols=[0], squeeze=True, dtype=str)
print(samples)
##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),
    contig="|".join(contigs)
##################################################################################
SAMPLES=list(units.index.levels[0])
MARK_SAMPLES=list(samples.index)

UNITS=list(units.index.levels[1])

CONTROL = config["control"]
CONTROLS = [sample for sample in MARK_SAMPLES if CONTROL in sample]
CASES = [sample for sample in MARK_SAMPLES if CONTROL not in sample]


## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))

ALL_SAMPLES = CASES + CONTROLS_UNIQUE     # same as MARK_SAMPLES
#ALL_BAM     = CONTROL_BAM + CASE_BAM


print("SAMPLES:::")
print(SAMPLES)

print("UNITS:::")
print(UNITS)

print("MARK_SAMPLES:::")
print(MARK_SAMPLES)

print("CONTROL:::")
print(CONTROL)

print("CONTROLS:::")
print(CONTROLS)

print("CASES:::")
print(CASES)

print("CONTROLS_UNIQUE:::")
print(CONTROLS_UNIQUE)

print("ALL_SAMPLES:::")
print(ALL_SAMPLES)
print("WWWWWWWWWWWWWWW")
print(units.loc[("BCC2","control"),  ["fq1", "fq2"]].dropna())
print("WWWWWWWWWWWWWWW")
##################################################################################


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)
