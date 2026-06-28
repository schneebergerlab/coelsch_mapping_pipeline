import math
import gzip
from glob import glob
import itertools as it


include: './align_common.snakefile'


seqtk_resources = rule_resources('seqtk')


rule cell_barcode_rc:
    '''reverse complement a fastq file, necessary for some snATAC-seq datasets sequenced on BGI nanoball'''
    input:
        barcode=raw_data(f'{{sample_name}}{cfg.file_suffixes.barcode}'),
    output:
        barcode=raw_data(f'{{sample_name}}{cfg.file_suffixes.barcode_rc}'),
    threads: seqtk_resources.threads
    resources:
        mem_mb=seqtk_resources.mem_mb,
    conda:
        get_conda_env('seqtk')
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
    ds = dataset(wc.dataset_name)
    tech_type = ds.technology
    ref_name = ds.reference_genotype
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
        if cfg.preprocessing.atac.rev_comp_barcode:
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
    tech_type = dataset(wc.dataset_name).technology
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
        umi_dedup_method=cfg.alignment.star.rna.umi_dedup_method,
    )

    return format_command(params.lstrip())


def get_transform_flag(wc):
    '''STAR genome transform flag - necessary if index has a VCF transformation'''
    ref = dataset(wc.dataset_name).reference_genotype
    if wc.qry != ref:
        return '--genomeTransformOutput SAM'
    else:
        return ''


def get_input_flags(wc, input):
    '''
    create flags for multiple input files, check if they are gzipped and add readFilesCommand
    '''
    tech_type = dataset(wc.dataset_name).technology
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
    tech_type = dataset(wc.dataset_name).technology
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
        mates_gap_max=cfg.alignment.star.atac.mates_gap_max,
        align_sj_overhang_min=cfg.alignment.star.rna.align_sj_overhang_min,
        align_sjdb_overhang_min=cfg.alignment.star.rna.align_sjdb_overhang_min,
        align_intron_min=cfg.alignment.star.rna.align_intron_min,
        align_intron_max=cfg.alignment.star.rna.align_intron_max
    )
    return format_command(params.lstrip())


star_align_resources = rule_resources('star/align/droplet')


rule STAR_consensus:
    '''
    map reads for a single-cell dataset using STARsolo and a VCF-transformed genome index.
    '''
    input:
        unpack(STAR_consensus_input)
    output:
        bam=temp(results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam')),
        bai=temp(results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam.bai')),
    params:
        input_flag=get_input_flags,
        transform_flag=get_transform_flag,
        adapter_parameters=get_adapter_parameters,
        splicing_parameters=get_spliced_alignment_params,
        filter_multimap_nmax=cfg.alignment.star.filter_multimap_nmax,
        filter_mismatch_nmax=cfg.alignment.star.filter_mismatch_nmax,
        sort_mem=lambda wc, resources: (resources.mem_mb - 4096) * 1_000_000,
        n_files=lambda wc, threads: threads * 150 + 200,
    log:
        progress=results('logs/{dataset_name}.{qry}.STAR_progress.log'),
        final=results('logs/{dataset_name}.{qry}.STAR_final.log'),
        main=results('logs/{dataset_name}.{qry}.STAR.log')
    threads: star_align_resources.threads
    resources:
        mem_mb=star_align_resources.mem_mb,
    conda:
        get_conda_env('star')
    shell:
        format_command(setup_tmpdir() + '''
        RELPATH=$(realpath --relative-to="$tmpdir" ".");
        cd $tmpdir;
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
        mv $tmpdir/Aligned.sortedByCoord.out.bam
          {output.bam};
        samtools index {output.bam};
        mv $tmpdir/Log.progress.out
          {log.progress};
        mv $tmpdir/Log.final.out
          {log.final};
        mv $tmpdir/Log.out
          {log.main};

        ''')


samtools_sort_resources = rule_resources('samtools/sort')


rule sort_bam_by_name:
    '''
    name-sort the output of star consensus to ensure consistent order of read-ids across haplotypes
    '''
    input:
        bam=results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam'),
        bai=results('aligned_data/haploid/{dataset_name}.{qry}.sorted.bam.bai')
    output:
        bam=temp(results('aligned_data/haploid/{dataset_name}.{qry}.namesorted.bam'))
    threads: samtools_sort_resources.threads
    resources:
        mem_mb=samtools_sort_resources.mem_mb,
        threads_per_cmd=max(1, samtools_sort_resources.threads // 2)
    conda:
        get_conda_env('htslib')
    shell:
        format_command(setup_tmpdir() + '''
        samtools sort
          -T ${{tmpdir}}/{wildcards.dataset_name}.{wildcards.qry}.sort_tmp
          -n -@ {resources.threads_per_cmd}
          {input.bam} |
        samtools fixmate -@ {resources.threads_per_cmd} -m - {output.bam};
        ''')


def get_merge_input(wc):
    '''Expand all haplotypes that have been aligned to for a dataset to create merge input'''
    qrys = get_dataset_haplotypes(wc.dataset_name)
    return {
        'bams': expand(results('aligned_data/haploid/{dataset_name}.{qry}.namesorted.bam'),
                       dataset_name=wc.dataset_name, qry=qrys),
    }


samtools_merge_resources = rule_resources('samtools/merge/dataset')


rule merge_name_sorted_bams:
    '''
    merge all haplotype-specific alignments into a single name-sorted bam file for haplotype collapsing
    '''
    input:
        unpack(get_merge_input)
    output:
        bam=temp(results('aligned_data/{dataset_name}.namesorted.bam')),
    threads: samtools_merge_resources.threads
    resources:
        mem_mb=samtools_merge_resources.mem_mb,
    conda:
        get_conda_env('htslib')
    shell:
        format_command('''
        samtools merge -@ {threads} -n
          {output.bam}
          {input.bams};
        ''')


collapse_aln_resources = rule_resources('samtools/sort')


rule collapse_alignments:
    '''
    Uses coelsch script collapse_ha_specific_alns.py to select the best haplotype alignment(s) for each read
    and outputs a single alignment with a new ha tag that indicates which haplotype(s) is/are best
    '''
    input:
        bam=results('aligned_data/{dataset_name}.namesorted.bam')
    output:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai')
    threads: collapse_aln_resources.threads
    resources:
        mem_mb=collapse_aln_resources.mem_mb,
    conda:
        get_conda_env('coelsch')
    shell:
        format_command(setup_tmpdir() + '''
        collapse_ha_specific_alns.py
          -o ${{tmpdir}}/{wildcards.dataset_name}.unsorted.bam
          {input.bam};

        samtools sort -@ {threads}
          -T ${{tmpdir}}/{wildcards.dataset_name}.sort_tmp
          -o {output.bam}
          ${{tmpdir}}/{wildcards.dataset_name}.unsorted.bam;

        samtools index {output.bam};
        ''')