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
        read1_glob = raw_data(f'{basename_glob}{cfg.file_suffixes.read1}')
        # create a snakemake glob_wildcards style glob for extracting the fully expanded sample_names
        sn_wildcards_read1_glob = raw_data(f'{{sample_name}}{cfg.file_suffixes.read1}')
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
    vcf_fns = config['annotations']['vcf_fns']
    ref, *qry_names = wc.geno_group.split('_')
    if wc.ref != wc.qry:
        if ref in vcf_fns and vcf_fns[ref] is not None:
            input_['vcf'] = annotation('vcf/star_consensus/user_defined/{geno_group}.{qry}.vcf')
        else:
            input_['vcf'] = annotation('vcf/star_consensus/syri/{geno_group}.{qry}.vcf')
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
    if cfg.alignment.star.sjdb_overhang_length is not None:
        return cfg.alignment.star.sjdb_overhang_length
    fastq_fns = glob(
        raw_data(f'*{cfg.file_suffixes.read2}')
    )[:5] # just check a few random files or this gets slow
    max_overhang = 0
    for fq_fn in fastq_fns:
        max_overhang = max(max_overhang, get_read_length(fq_fn))
    return max_overhang - 1


def calculate_genome_sa_index_nbases(wc, input):
    '''estimate the suffix array nbases for star from the size of the reference genome'''
    fai_path = str(input.fasta) + '.fai'
    if not os.path.exists(fai_path):
        if cfg.alignment.star.genome_sa_index_nbases is not None:
            return cfg.alignment.star.genome_sa_index_nbases
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


star_index_resources = rule_resources('star/index')


rule build_STAR_index:
    """
    Builds a STAR genome index for alignment. Uses either reference genome alone, or in conjunction with a 
    VCF file to produce a "STAR consensus" variant transformed genome.
    """
    input:
        unpack(get_star_index_input)
    output:
        directory(annotation('star_indexes/{geno_group}/{ref}.{qry}.star_idx'))
    threads: star_index_resources.threads
    resources:
        mem_mb=star_index_resources.mem_mb,
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

        STARTMP="_STARtmp_{wildcards.geno_group}.{wildcards.ref}.{wildcards.qry}";
        [ -n "$STARTMP" ] && rm -rf -- "$STARTMP";
        trap '[ -n "$STARTMP" ] && rm -rf -- "$STARTMP"' EXIT;

        STAR
          --outTmpDir "$STARTMP"
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


def get_dataset_haplotypes(dataset_name):
    '''return founder haplotypes represented by a dataset crossing design.'''
    ds = dataset(dataset_name)
    return founder_haplotypes_from_crossing_combinations(ds.crossing_combinations)


def get_geno_group(dataset_name):
    '''creates a geno_group wildcard from the config by joining together all the genotypes'''
    ds = dataset(dataset_name)
    ref = ds.reference_genotype
    qry_names = sorted(hap for hap in get_dataset_haplotypes(dataset_name) if hap != ref)
    return '_'.join([ref, *qry_names])


def get_star_fastq_input(dataset_name, fastq_type):
    '''use globbing to identify all the input fastqs for a dataset from input_file_basenames'''
    sample_names = DATASET_SAMPLE_NAME_MAPPING[dataset_name]
    return expand(
        raw_data('{sample_name}{file_suffix}'),
        sample_name=sample_names,
        file_suffix=getattr(cfg.file_suffixes, fastq_type)
    )