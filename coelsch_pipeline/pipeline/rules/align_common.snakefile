import os
import math
import gzip
from glob import glob
from collections import defaultdict


SAMPLE_NAME_DATASET_MAPPING = {}
DATASET_SAMPLE_NAME_MAPPING = defaultdict(list)
for dataset_name in config['datasets']:
    # config only stores dataset_name -> sample_name_glob relationship
    # use globbing to determine the sample_name -> dataset_name relationship
    for basename_glob in config['datasets'][dataset_name]['input_file_basenames']:
        # create glob.glob style globbing string for identifying files matching the input_file_basenames pattern
        read1_glob = raw_data(f'{basename_glob}{config["file_suffixes"]["read1"]}')
        # create a snakemake glob_wildcards style glob for extracting the fully expanded sample_names
        sn_wildcards_read1_glob = raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}')
        sample_names = glob_wildcards(sn_wildcards_read1_glob, files=glob(read1_glob)).sample_name
        # add the sample names to the reference mappings
        for sn in sample_names:
            SAMPLE_NAME_DATASET_MAPPING[sn] = dataset_name
            DATASET_SAMPLE_NAME_MAPPING[dataset_name].append(sn)


def get_star_index_input(wc):
    '''input for STAR genome index generation'''
    input_ = {
        'fasta': get_fasta(wc.ref),
        'gtf': get_gtf(wc.ref),
    }

    if wc.ref != wc.qry:
        input_['vcf'] = annotation('vcf/star_consensus/msyd/{geno_group}.{qry}.vcf')
    return input_


def get_read_length(fastq_fn):
    '''get assumed length of reads in a fastq file from the first read'''
    if fastq_fn.endswith('.gz'):
        open_method = gzip.open
    else:
        open_method = open
    with open_method(fastq_fn) as f:
        next(f) # skip first read id
        first_read_seq = next(f).decode().strip()
    return len(first_read_seq)


def get_sj_overhang_size(wc):
    '''get the splice junction overhang size for a genome index from the data'''
    if config['alignment']['star'].get('sjdb_overhang_length') is not None:
        return config['alignment']['star']['sjdb_overhang_length']
    fastq_fns = glob(
        raw_data(f'*{config["file_suffixes"]["read2"]}')
    )[:5] # just check a few random files or this gets slow
    max_overhang = 0
    for fq_fn in fastq_fns:
        max_overhang = max(max_overhang, get_read_length(fq_fn))
    return max_overhang - 1


def calculate_genome_sa_index_nbases(wc, input):
    '''estimate the suffix array nbases for star from the size of the reference genome'''
    fai_path = str(input.fasta) + '.fai'
    if not os.path.exists(fai_path):
        if config['alignment']['star'].get('genome_sa_index_nbases') is not None:
            return config['alignment']['star']['genome_sa_index_nbases']
        try:
            import pysam
            pysam.index(input.fasta)
        except:
            raise ValueError(
                f'Cannot find/create fai file for "{input.fasta}". Please provide either .fai or '
                'specify alignment:star:genome_sa_index_nbases in the config file'
            )
    total_length = 0
    with open(fai_path) as f:
        for line in f:
            _, ln, *_ = line.strip().split('\t')
            total_length += int(ln)
    return min(14, math.floor(math.log2(total_length) / 2 - 1))


rule build_STAR_index:
    """
    Build a STAR genome index for alignment, optionally using STAR genome transform.

    Creates a STAR index either from the reference genome alone, or (when ref genotype != qry genotype)
    from the reference genome plus a per-query VCF to generate a haploid “STAR consensus”
    transformed genome index.

    Parameters:
      Defined in config section alignment: star (and uses VCFs produced under variants: star_consensus).

    Inputs:
      - reference genome fasta (annotations)
      - reference genome annotation GTF (annotations)
      - optional query VCF for genome transformation (annotations/vcf/star_consensus)

    Outputs:
      - STAR genome index directory for the reference × query combination (annotations/star_indexes)
    """
    input:
        unpack(get_star_index_input)
    output:
        directory(annotation('star_indexes/{geno_group}/{ref}.{qry}.star_idx'))
    threads: 12
    resources:
        mem_mb=12 * 1024,
    params:
        overhang=get_sj_overhang_size,
        genome_sa_index_nbases=calculate_genome_sa_index_nbases,
        vcf_flag=lambda wc, input: f'--genomeTransformVCF {input.vcf}' if hasattr(input, 'vcf') else '',
        transform_flag=lambda wc, input: '--genomeTransformType Haploid' if hasattr(input, 'vcf') else '',
    conda:
        get_conda_env('star')
    shell:
        format_command('''
        mkdir {output};
        STAR
          --outTmpDir _STARtmp_{wildcards.geno_group}.{wildcards.ref}.{wildcards.qry}
          --runThreadN {threads}
          --runMode genomeGenerate
          --genomeDir {output}
          --genomeFastaFiles {input.fasta}
          --sjdbGTFfile {input.gtf}
          {params.vcf_flag}
          {params.transform_flag}
          --genomeSAindexNbases {params.genome_sa_index_nbases}
          --sjdbOverhang {params.overhang}
        ''')


def get_geno_group(dataset_name):
    '''creates a geno_group wildcard from the config by joining together all the genotypes'''
    dataset = config['datasets'][dataset_name]
    ref = dataset['reference_genotype']
    qry_names = set()
    for geno in dataset['genotypes'].values():
        for qry in geno['founder_haplotypes'].values():
            if qry != ref:
                qry_names.add(qry)
    qry_names = '_'.join(sorted(qry_names))
    return f'{ref}_{qry_names}'


def get_star_fastq_input(dataset_name, fastq_type):
    '''use globbing to identify all the input fastqs for a dataset from input_file_basenames'''
    sample_names = DATASET_SAMPLE_NAME_MAPPING[dataset_name]
    return expand(
        raw_data('{sample_name}{file_suffix}'),
        sample_name=sample_names,
        file_suffix=config['file_suffixes'][fastq_type]
    )