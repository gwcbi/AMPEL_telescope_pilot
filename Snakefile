#! /usr/bin/env python
from future import standard_library
standard_library.install_aliases()


import re
import gzip
from os import path
from os.path import abspath, dirname
from subprocess import check_output


import snakemake


def dict_args(d):
    """ Convert a dictionary to command line arguments"""
    ret = ''
    for k,v in d.items():
        if len(k) == 1:
            ret += '-{} {} '.format(k,v)
        else:
            ret += '--{} {} '.format(k,v)
    return ret


""" Constants """


SAMPLES = [l.strip('\n') for l in open('samples.txt')]

configfile: "config.yaml"

wildcard_constraints:
    sampid="\w+\d+",

localrules: all, complete_project, complete_sample, structure, download_telescope_annotation

rule all:
    input:        
        expand("{s}/completed.txt", s=SAMPLES)

rule complete_sample:
    output:
        touch("{sampid}/completed.txt")
    input:
        "{sampid}/inform-telescope_report.tsv",


rule flexbar:
    input:
        "{sampid}/read_1.fq",
        "{sampid}/read_2.fq"
    output:
        "{sampid}/clean_1.fastq.gz",
        "{sampid}/clean_2.fastq.gz"
    params:
        outdir = "{sampid}"
    threads: snakemake.utils.available_cpu_count()
    run:
        flexbar_args = dict_args(config['flexbar'])
        adapters = config['adapters']['PE']
        shell("module load flexbar/3.0.3 && "
            "flexbar "
            "--threads {threads} "
            "{flexbar_args} "
            "--adapters {adapters} "
            "--reads {input[0]} "
            "--reads2 {input[1]} "            
            "--target {params.outdir}/clean "
        )

rule sortbam:
    input:
        "{f}.unsorted.bam"
    output:
        "{f}.sorted.bam",
        "{f}.sorted.bam.bai"
    threads: snakemake.utils.available_cpu_count()
    shell:
        "samtools sort -@ {threads} {input} > {output[0]} && "
        "samtools index {output[0]}"


rule bowtie2:
    input:
        "{sampid}/clean_1.fastq.gz",
        "{sampid}/clean_2.fastq.gz",
        expand(config['bt2idx'] + ".{i}.bt2", i = range(1, 5)),
        expand(config['bt2idx'] + ".rev.{i}.bt2", i = range(1, 3))
    output:
        "{sampid}/bt2_{preset, \w+}.unsorted.bam",
        "{sampid}/bt2_{preset, \w+}.summary.txt"
    threads: snakemake.utils.available_cpu_count()        
    run:
        if wildcards.preset in config['bowtie2']:
            bowtie2_args = dict_args(config['bowtie2'][wildcards.preset])
        else:
            bowtie2_args = ""
        shell(
            "(bowtie2 "
            "-p {threads} "
            "{bowtie2_args} "
            "--rg-id $(basename $(dirname {input[0]})) "
            "-x {config[bt2idx]} "
            "-1 {input[0]} -2 {input[1]} | "
            "samtools view -b > {output[0]}"
            ") 3>&1 1>&2 2>&3 | tee {output[1]} "
        )


rule telescope:
    input:
        "{sampid}/bt2_multi.unsorted.bam",
        config['herv_annotation']
    output:
        "{sampid}/{preset,\w+}-telescope_report.tsv",
        "{sampid}/{preset,\w+}.log"
    run:
        if wildcards.preset in config['telescope']:
            telescope_args = dict_args(config['telescope'][wildcards.preset])
        else:
            telescope_args = ""
        shell(
            "telescope "
            "assign "
            "{telescope_args} "            
            "--outdir $(dirname {output[0]}) "
            "{input[0]} "
            "{input[1]} "
            "2>&1 | tee {output[1]}"
        )

# rule kallisto:
#     input:
#         "{sraproj}/{sraid}/clean_1.fastq.gz",
#         "{sraproj}/{sraid}/clean_2.fastq.gz",
#     output:
#         "{sraproj}/{sraid}/abundance.tsv",
#         "{sraproj}/{sraid}/abundance.h5"        
#     threads: snakemake.utils.available_cpu_count()
#     run:
#         if ispaired(input[0], input[1]):
#             rargs = '{} {}'.format(input[0], input[1])        
#         else:
#             rargs = '-l 250 -s 50 --single {}'.format(input[0])
#         
#         shell(
#             "kallisto "
#             "quant "
#             "-t {threads} "
#             "-b 100 "
#             "-i {config[kallistoidx]} "            
#             "-o $(dirname {output[0]}) "
#             "{rargs} "
#         )
# 
