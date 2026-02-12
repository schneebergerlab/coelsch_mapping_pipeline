import math
import gzip
from glob import glob
import itertools as it


include: './align_common.snakefile'


rule cell_barcode_rc:
    """
    Reverse-complement barcode FASTQs for sequencing layouts that require it.

    Some snATAC-seq datasets (e.g. certain BGI/Nanoball runs) encode the barcode read
    in the opposite orientation. This rule reverse-complements the barcode FASTQ and
    bgzips the result for downstream alignment.

    Parameters:
      Defined in config section preprocessing: atac (rev_comp_barcode) and file_suffixes.

    Inputs:
      - barcode fastq for each sample (raw_data)

    Outputs:
      - reverse-complemented barcode fastq (raw_data)
    """
    input:
        barcode=raw_data(f'{{sample_name}}{config["file_suffixes"]["barcode"]}'),
    output:
        barcode=raw_data(f'{{sample_name}}{config["file_suffixes"]["barcode_rc"]}'),
    conda:
        '../env_yamls/seqtk.yaml'
    shell:
        format_command('''
        seqtk seq -r {input.barcode} |
          bgzip > {output.barcode};
        ''')


def STAR_consensus_input(wc):
    '''
    full input to STAR consensus (droplet mode with STARsolo)
    input is:
      - index: the STAR index generated from the reference genome plus VCF transform
      - barcode_whitelist: the whitelist file(s) for the sequencing technology type
      - read, mate, barcode: the fastq files for the dataset - contains ALL cells/barcodes
    '''
    dataset = config['datasets'][wc.dataset_name]
    tech_type = dataset['technology']
    ref_name = dataset['reference_genotype']
    geno_group = get_geno_group(wc.dataset_name)
    barcode_whitelist = get_barcode_whitelist(tech_type)
    if isinstance(barcode_whitelist, str):
        barcode_whitelist = [barcode_whitelist,]
    input_ = {
        'index': annotation(f'star_indexes/{geno_group}/{ref_name}.{{qry}}.star_idx'),
        'barcode_whitelist': barcode_whitelist,
    }
    if tech_type == "10x_atac":
        input_['read'] = get_star_fastq_input(wc.dataset_name, 'read1')
        input_['mate'] = get_star_fastq_input(wc.dataset_name, 'read2')
        if config['preprocessing']['atac']['rev_comp_barcode']:
            input_['barcode'] = get_star_fastq_input(wc.dataset_name, 'barcode_rc')
        else:
            input_['barcode'] = get_star_fastq_input(wc.dataset_name, 'barcode')
    else:
        input_['read'] = get_star_fastq_input(wc.dataset_name, 'read2')
        input_['barcode'] = get_star_fastq_input(wc.dataset_name, 'read1')
    return input_


def get_adapter_parameters(wc, input):
    '''
    Adapter parameters for the different sequencing modalities
    '''
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type == '10x_atac':
        params = '''
          --soloType "CB_samTagOut"
          --soloCBwhitelist {whitelist}
          --soloBarcodeReadLength 0
          --soloCBmatchWLtype "1MM"
          --outSAMattributes "NH" "HI" "AS" "nM" "RG" "CB"
        '''
    elif tech_type.startswith('10x_rna'):
        params = '''
          --soloType "CB_UMI_Simple"
          --soloCBwhitelist {whitelist}
          --soloBarcodeReadLength {barcode_read_length}
          --soloUMIdedup {umi_dedup_method}
          --soloCBmatchWLtype "1MM"
          --soloCBlen 16
          --soloCBstart 1
          --soloUMIlen 12
          --soloUMIstart 17
          --outSAMattributes "NH" "HI" "AS" "nM" "RG" "CB" "UB"
        '''
    elif tech_type == 'bd_rna':
        params = '''
          --soloType "CB_UMI_Complex"
          --soloCBwhitelist {whitelist}
          --soloBarcodeReadLength {barcode_read_length}
          --soloUMIdedup {umi_dedup_method}
          --soloAdapterSequence "NNNNNNNNNGTGANNNNNNNNNGACA"
          --soloCBmatchWLtype "1MM"
          --soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9
          --soloUMIposition 3_10_3_17
          --outSAMattributes "NH" "HI" "AS" "nM" "RG" "CB" "UB"
        '''
    else:
        raise NotImplementedError()

    whitelist = ' '.join(f'${{RELPATH}}/{fn}' for fn in input.barcode_whitelist)
    params = params.format(
        whitelist=whitelist,
        barcode_read_length=get_read_length(input.barcode[0]),
        umi_dedup_method=config["alignment"]["star"]["rna"]["umi_dedup_method"],
    )

    return format_command(params.lstrip())


def get_transform_flag(wc):
    '''STAR genome transform flag - necessary if index has a VCF transformation'''
    dataset = config['datasets'][wc.dataset_name]
    ref = dataset['reference_genotype']
    if wc.qry != ref:
        return '--genomeTransformOutput SAM'
    else:
        return ''


def get_temp_dir(wc, output):
    '''
    create a private directory for each dataset_name/query combination, since STAR does not work well 
    when multiple processes are run in the same directory
    '''
    output_dir = os.path.split(output.bam)[0]
    tmp_dir = os.path.join(output_dir, f'{wc.dataset_name}.{wc.qry}.tmpdir')
    return tmp_dir


def get_input_flags(wc, input):
    '''
    create flags for multiple input files, check if they are gzipped and add readFilesCommand
    '''
    tech_type = config['datasets'][wc.dataset_name]['technology']
    all_fastqs = it.chain(*(
        (input.read, input.mate, input.barcode)
        if tech_type == '10x_atac'
        else (input.read, input.barcode)
    ))
    if all([fn.endswith('.gz') for fn in all_fastqs]):
        flag = '''
          --readFilesCommand "zcat" 
          --readFilesIn'''
    else:
        flag = '''
          --readFilesIn'''
    flag += ' '
    flag += ','.join(f'${{RELPATH}}/{fn}' for fn in input.read)
    if tech_type == '10x_atac':
        flag += ' '
        flag += ','.join(f'${{RELPATH}}/{fn}' for fn in input.mate)
    flag += ' '
    flag += ','.join(f'${{RELPATH}}/{fn}' for fn in input.barcode)
    return format_command(flag.lstrip())
    

def get_spliced_alignment_params(wc):
    '''splicing/fragment size parameters for STAR'''
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type == "10x_atac":
        params = '''
          --alignIntronMax 1
          --alignMatesGapMax {mates_gap_max}
        '''
    else:
        params = '''
          --outFilterIntronMotifs RemoveNoncanonical
          --alignSJoverhangMin {align_sj_overhang_min}
          --alignSJDBoverhangMin {align_sjdb_overhang_min}
          --alignIntronMin {align_intron_min}
          --alignIntronMax {align_intron_max}
        '''
    params = params.format(
        mates_gap_max=config["alignment"]["star"]["atac"]["mates_gap_max"],
        align_sj_overhang_min=config["alignment"]["star"]["rna"]["align_sj_overhang_min"],
        align_sjdb_overhang_min=config["alignment"]["star"]["rna"]["align_sjdb_overhang_min"],
        align_intron_min=config["alignment"]["star"]["rna"]["align_intron_min"],
        align_intron_max=config["alignment"]["star"]["rna"]["align_intron_max"]
    )
    return format_command(params.lstrip())


rule STAR_consensus:
    """
    Align droplet-based single-cell reads with STARsolo against haplotype transformed indices.

    Runs STARsolo in droplet mode for one dataset_name × haplotype combination. Uses the
    appropriate barcode/UMI parsing configuration for the dataset technology, and
    emits a coordinate-sorted BAM plus STAR logs. If the index was built with a VCF
    transform (qry != ref), emits SAM records with genome transform output enabled.

    Parameters:
      Defined in config sections:
        - alignment: star (shared alignment settings)
        - alignment: star: rna / atac (technology-specific settings)
        - datasets: <dataset_name> (technology type, reference genotype, genotypes)
        - preprocessing: atac (barcode reverse-complement behaviour)
        - file_suffixes (fastq naming)

    Inputs:
      - STAR index directory for reference × query combination (annotations/star_indexes)
      - barcode whitelist file(s) for the sequencing technology
      - dataset fastqs (reads, and mates/barcodes depending on technology) (raw_data)

    Outputs:
      - temporary coordinate-sorted BAM and index for this dataset × haplotype (results/aligned_data/haploid)
      - STAR run logs (results/logs)
    """
    input:
        unpack(STAR_consensus_input)
    output:
        bam=temp(results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam')),
        bai=temp(results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam.bai')),
    params:
        star_tmp_dir=get_temp_dir,
        input_flag=get_input_flags,
        transform_flag=get_transform_flag,
        adapter_parameters=get_adapter_parameters,
        splicing_parameters=get_spliced_alignment_params,
        filter_multimap_nmax=config['alignment']['star']['filter_multimap_nmax'],
        filter_mismatch_nmax=config['alignment']['star']['filter_mismatch_nmax'],
        sort_mem=lambda wc, resources: (resources.mem_mb - 4096) * 1_000_000,
        n_files=lambda wc, threads: threads * 150 + 200,
    log:
        progress=results('logs/{dataset_name}.{qry}.STAR_progress.log'),
        final=results('logs/{dataset_name}.{qry}.STAR_final.log'),
        main=results('logs/{dataset_name}.{qry}.STAR.log')
    threads: 24
    resources:
        mem_mb=lambda wildcards, threads: (threads + 4) * 2048,
    conda:
        get_conda_env('star')
    shell:
        format_command('''
        mkdir -p {params.star_tmp_dir};
        RELPATH=$(realpath --relative-to="{params.star_tmp_dir}" ".");
        cd {params.star_tmp_dir};
        ulimit -n {params.n_files};
        STAR
          --runThreadN {threads}
          --genomeDir "${{RELPATH}}/{input.index}"
          {params.input_flag}
          {params.adapter_parameters}
          {params.splicing_parameters}
          --outFilterMultimapNmax {params.filter_multimap_nmax}
          --outFilterMismatchNmax {params.filter_mismatch_nmax}
          --outSAMtype "BAM" "SortedByCoordinate"
          --outBAMsortingBinsN 150
          --limitBAMsortRAM {params.sort_mem}
          {params.transform_flag}
          --outSAMattrRGline "ID:{wildcards.qry}";

        cd $RELPATH;
        mv {params.star_tmp_dir}/Aligned.sortedByCoord.out.bam
          {output.bam};
        samtools index {output.bam};
        mv {params.star_tmp_dir}/Log.progress.out
          {log.progress};
        mv {params.star_tmp_dir}/Log.final.out
          {log.final};
        mv {params.star_tmp_dir}/Log.out
          {log.main};

        rm -rf {params.star_tmp_dir};
        ''')


rule sort_bam_by_name:
    """
    Name-sort STAR output to guarantee consistent read order across haplotypes.

    Converts coordinate-sorted STAR BAMs into name-sorted BAMs so that
    read-id ordering is consistent across haplotype-specific alignments. This is required
    for deterministic merging and subsequent per-read haplotype selection.

    Fixmate is also run to update mate information for paired-ended snATAC-seq data.

    Inputs:
      - coordinate-sorted haplotype BAM and index (results/aligned_data/haploid)

    Outputs:
      - temporary name-sorted haplotype BAM (results/aligned_data/haploid)
    """
    input:
        bam=results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam'),
        bai=results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam.bai')
    output:
        bam=temp(results('aligned_data/haploid/{dataset_name}.{qry}.namesorted.bam'))
    threads: 12
    resources:
        mem_mb=20_000,
        threads_per_cmd=lambda wc, threads: threads // 2
    conda:
        get_conda_env('htslib')
    shell:
        format_command('''
        samtools sort
          -T ${{TMPDIR}}/{wildcards.dataset_name}.{wildcards.qry}
          -n -@ {resources.threads_per_cmd}
          {input.bam} |
        samtools fixmate -@ {resources.threads_per_cmd} -m - {output.bam};
        ''')


def get_merge_input(wc):
    '''Expand all haplotypes that have been aligned to for a dataset to create merge input'''
    dataset = config['datasets'][wc.dataset_name]
    qrys = set()
    for geno in dataset['genotypes'].values():
        for qry in geno['founder_haplotypes'].values():
            qrys.add(qry)
    return {
        'bams': expand(results('aligned_data/haploid/{dataset_name}.{qry}.namesorted.bam'),
                       dataset_name=wc.dataset_name, qry=qrys),
    }


rule merge_name_sorted_bams:
    """
    Merge per-haplotype name-sorted BAMs into a single name-sorted BAM.

    Combines all haplotype-specific alignments for a dataset into one BAM while
    preserving name-sorted order, producing the input required for haplotype
    collapsing / best-alignment selection.

    Inputs:
      - name-sorted BAMs for each haplotype aligned for this dataset (results/aligned_data/haploid)

    Outputs:
      - temporary merged name-sorted BAM for the dataset (results/aligned_data)
    """
    input:
        unpack(get_merge_input)
    output:
        bam=temp(results('aligned_data/{dataset_name}.namesorted.bam')),
    threads: 24
    resources:
        mem_mb=lambda wildcards, threads: threads * 1024,
    conda:
        get_conda_env('htslib')
    shell:
        format_command('''
        samtools merge -@ {threads} -n
          {output.bam}
          {input.bams};
        ''')


rule collapse_alignments:
    """
    Collapse haplotype-specific alignments into a single best alignment per read.

    Uses the coelsch collapse_ha_specific_alns.py script to choose the best haplotype
    alignment(s) for each read from the merged name-sorted BAM, then sorts and indexes
    the resulting BAM. Adds/updates a haplotype tag [ha] indicating which haplotype(s) were
    selected as best for each read.

    Inputs:
      - merged name-sorted BAM containing alignments to all haplotypes (results/aligned_data)

    Outputs:
      - final coordinate-sorted haplotype-labelled BAM and index for the dataset (results/aligned_data)
    """
    input:
        bam=results('aligned_data/{dataset_name}.namesorted.bam')
    output:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai')
    resources:
        mem_mb=20_000,
    threads: 12
    conda:
        get_conda_env('coelsch')
    shell:
        format_command('''
        collapse_ha_specific_alns.py
          -o {output.bam}.unsorted.bam
          {input.bam};

        samtools sort -@ {threads}
          -T ${{TMPDIR}}/{wildcards.dataset_name}
          -o {output.bam}
          {output.bam}.unsorted.bam;

        samtools index {output.bam};

        rm {output.bam}.unsorted.bam
        ''')