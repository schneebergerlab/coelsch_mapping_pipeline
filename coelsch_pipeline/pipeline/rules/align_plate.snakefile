from glob import glob

include: './align_common.snakefile'


def STAR_consensus_input(wc):
    '''
    full input to STAR consensus (plate mode WITHOUT STARsolo).
    input is:
      - index: the STAR index generated from the reference genome plus VCF transform
      - read, mate: the fastq files for the sample - represents a single individual/barcode
    '''
    dataset_name = SAMPLE_NAME_DATASET_MAPPING[wc.sample_name]
    ds = dataset(dataset_name)
    ref_name = ds.reference_genotype
    geno_group = get_geno_group(dataset_name)
    return {
        'index': annotation(f'star_indexes/{geno_group}/{ref_name}.{{qry}}.star_idx'),
        'read': raw_data(f'{{sample_name}}{cfg.file_suffixes.read1}'),
        'mate': raw_data(f'{{sample_name}}{cfg.file_suffixes.read2}'),
    }


def get_readfiles_cmd(wc, input):
    '''
    create flags files, checking if they are gzipped and add readFilesCommand
    '''
    if all([fn.endswith('.gz') for fn in (input.read, input.mate)]):
        return '--readFilesCommand "zcat"'
    else:
        return ''


def get_transform_flag(wc):
    '''STAR genome transform flag - necessary if index has a VCF transformation'''
    dataset_name = SAMPLE_NAME_DATASET_MAPPING[wc.sample_name]
    ref = dataset(dataset_name).reference_genotype
    if wc.qry != ref:
        return '--genomeTransformOutput "SAM"'
    else:
        return ''


star_align_resources = rule_resources('star/align/plate')


rule STAR_consensus:
    '''
    map reads for a single individual/barcode using STARsolo and a VCF-transformed genome index.
    '''
    input:
        unpack(STAR_consensus_input)
    output:
        bam=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam')),
        bai=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam.bai')),
    params:
        readfiles_cmd=get_readfiles_cmd,
        transform_flag=get_transform_flag,
        filter_multimap_nmax=cfg.alignment.star.filter_multimap_nmax,
        filter_mismatch_nmax=cfg.alignment.star.filter_mismatch_nmax,
        align_mates_gap_max=cfg.alignment.star.atac.mates_gap_max,
    log:
        progress=results('logs/{sample_name}.{qry}.STAR_progress.log'),
        final=results('logs/{sample_name}.{qry}.STAR_final.log'),
        main=results('logs/{sample_name}.{qry}.STAR.log')
    threads: star_align_resources.threads
    resources:
        mem_mb=star_align_resources.mem_mb,
    conda:
        get_conda_env('star')
    shell:
        format_command(setup_tmpdir() + '''
        RELPATH=$(realpath --relative-to="$tmpdir" ".");
        cd $tmpdir;

        STAR
          --runThreadN {threads}
          --genomeDir "${{RELPATH}}/{input.index}"
          {params.readfiles_cmd}
          --readFilesIn "${{RELPATH}}/{input.read}" "${{RELPATH}}/{input.mate}"
          --alignIntronMax 1
          --alignMatesGapMax {params.align_mates_gap_max}
          --outFilterMultimapNmax {params.filter_multimap_nmax}
          --outFilterMismatchNmax {params.filter_mismatch_nmax}
          --outSAMtype "BAM" "Unsorted"
          --outSAMattrRGline "ID:{wildcards.qry}"
          {params.transform_flag}
          --outSAMattributes "NH" "HI" "AS" "nM" "RG";

        cd $RELPATH;
        mv $tmpdir/Log.progress.out {log.progress};
        mv $tmpdir/Log.final.out {log.final};
        mv $tmpdir/Log.out {log.main};

        samtools sort
          -m 1G -@ {threads}
          -o {output.bam}
          $tmpdir/Aligned.out.bam;

        samtools index {output.bam};
        ''')


samtools_sort_resources = rule_resources('samtools/sort')


rule sort_bam_by_name:
    '''
    name-sort the output of star consensus to ensure consistent order of read-ids across haplotypes
    '''
    input:
        bam=results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam'),
        bai=results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam.bai')
    output:
        bam=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'))
    threads: samtools_sort_resources.threads
    resources:
        mem_mb=samtools_sort_resources.mem_mb,
    conda:
        get_conda_env('htslib')
    shell:
        format_command(setup_tmpdir() + '''
        samtools sort -n -@ {threads}
          -T ${{tmpdir}}/{wildcards.sample_name}.{wildcards.qry}.sort_tmp
          {input.bam}
          > {output.bam};
        ''')


def get_merge_haps_input(wc):
    '''Expand all haplotypes that have been aligned to for a sample/barcode to create merge input'''
    dataset_name = SAMPLE_NAME_DATASET_MAPPING[wc.sample_name]
    qrys = get_dataset_haplotypes(dataset_name)
    return expand(
        results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'),
        sample_name=wc.sample_name, qry=qrys
    )


samtools_merge_resources = rule_resources('samtools/merge/sample')


rule merge_name_sorted_hap_bams:
    '''
    merge all haplotype-specific alignments for a sample/barcode into a single name-sorted bam file for haplotype collapsing
    '''
    input:
        bams=get_merge_haps_input
    output:
        bam=temp(results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')),
    threads: samtools_merge_resources.threads
    resources:
        mem_mb=samtools_merge_resources.mem_mb,
    conda:
        get_conda_env('htslib')
    shell:
        format_command('''
        samtools merge -@ {threads}
          -n {output.bam}
          {input.bams};
        ''')


collapse_aln_resources = rule_resources('collapse_alignments')


rule collapse_alignments:
    '''
    Uses coelsch script collapse_ha_specific_alns.py to select the best haplotype alignment(s) for each read
    and outputs a single alignment with a new ha tag that indicates which haplotype(s) is/are best.
    Also adds new RG tag to bam file to indicate the sample name, so that this can be identified later
    after sample-wise merging
    '''
    input:
        bam=results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')
    output:
        bam=results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        bai=results('aligned_data/single_barcodes/{sample_name}.sorted.bam.bai')
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

        samtools addreplacerg
          -r "ID:{wildcards.sample_name}"
          ${{tmpdir}}/{wildcards.dataset_name}.unsorted.bam |
        samtools sort -@ {threads}
          -T ${{tmpdir}}/{wildcards.sample_name}.sort_tmp
          -o {output.bam} - ;

        samtools index {output.bam};
        ''')


def merge_samples_input(wc):
    '''input for sample/barcode-wise merging'''
    sample_names = DATASET_SAMPLE_NAME_MAPPING[wc.dataset_name]
    return expand(
        results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        sample_name=sample_names
    )


rule list_bam:
    '''create a list of bam files representing all the samples for a dataset, for merging'''
    input:
        bams=merge_samples_input
    output:
        txt=temp(results("{dataset_name}.list"))
    run:
        with open(output.txt, 'w') as f:
            f.write('\n'.join(input.bams))


samtools_merge_ds_resources = rule_resources('samtools/merge/dataset')


rule merge_bams:
    '''sample-wise barcode merging'''
    input:
        results("{dataset_name}.list")
    output:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai'),
    conda:
        get_conda_env('htslib')
    threads: samtools_merge_ds_resources.threads
    resources:
        mem_mb=samtools_merge_ds_resources.mem_mb,
    shell:
        format_command('''
        ulimit -n 5000;

        samtools merge -@ {threads} -b {input} -o {output.bam};

        samtools index {output.bam};
        ''')
