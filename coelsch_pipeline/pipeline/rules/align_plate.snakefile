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
    dataset = config['datasets'][dataset_name]
    tech_type = dataset['technology']
    ref_name = dataset['reference_genotype']
    geno_group = get_geno_group(dataset_name)
    return {
        'index': annotation(f'star_indexes/{geno_group}/{ref_name}.{{qry}}.star_idx'),
        'read': raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}'),
        'mate': raw_data(f'{{sample_name}}{config["file_suffixes"]["read2"]}'),
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
    dataset = config['datasets'][dataset_name]
    ref = dataset['reference_genotype']
    if wc.qry != ref:
        return '--genomeTransformOutput "SAM"'
    else:
        return ''


def get_temp_dir(wc, output):
    '''
    create a private directory for each dataset_name/query combination, since STAR does not work well 
    when multiple processes are run in the same directory
    '''
    output_dir = os.path.split(output.bam)[0]
    tmp_dir = os.path.join(output_dir, f'{wc.sample_name}.{wc.qry}.tmpdir')
    return tmp_dir


rule STAR_consensus:
    """
    Align plate-based (single-individual/-barcode) reads with STAR against per-haplotype indices.

    Runs STAR in non-solo mode for one individual/barcode × haplotype combination, using the
    appropriate haplotype-variant-transformed STAR index. Emits a sorted BAM plus STAR logs.
 
    This is the plate-based analogue of the align_droplet STAR_consensus rule, but each sample_name
    corresponds to one individual/barcode, rather than the whole dataset.

    Parameters:
      Defined in config sections:
        - alignment: star (shared alignment settings)
        - alignment: star: atac (technology-specific settings)
        - datasets: <dataset_name> (technology type, reference genotype, genotypes)
        - file_suffixes (fastq naming)

    Inputs:
      - STAR index directory for reference × query combination (annotations/star_indexes)
      - paired-end FASTQs for the individual/barcode (read1/read2) (raw_data)

    Outputs:
      - temporary coordinate-sorted BAM and index for this individual × haplotype (results/aligned_data/single_barcodes/haploid)
      - STAR run logs (results/logs)
    """
    input:
        unpack(STAR_consensus_input)
    output:
        bam=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam')),
        bai=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam.bai')),
    params:
        star_tmp_dir=get_temp_dir,
        readfiles_cmd=get_readfiles_cmd,
        transform_flag=get_transform_flag,
        filter_multimap_nmax=config['alignment']['star']['filter_multimap_nmax'],
        filter_mismatch_nmax=config['alignment']['star']['filter_mismatch_nmax'],
        align_mates_gap_max=config['alignment']['star']['atac']['mates_gap_max'],
    log:
        progress=results('logs/{sample_name}.{qry}.STAR_progress.log'),
        final=results('logs/{sample_name}.{qry}.STAR_final.log'),
        main=results('logs/{sample_name}.{qry}.STAR.log')
    threads: 6 # usually small files, fewer threads needed than for droplets
    resources:
        mem_mb=lambda wildcards, threads: (threads + 4) * 2048,
    conda:
        get_conda_env('star')
    shell:
        format_command('''
        mkdir -p {params.star_tmp_dir};
        RELPATH=$(realpath --relative-to="{params.star_tmp_dir}" ".");
        cd {params.star_tmp_dir};

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
        mv {params.star_tmp_dir}/Log.progress.out {log.progress};
        mv {params.star_tmp_dir}/Log.final.out {log.final};
        mv {params.star_tmp_dir}/Log.out {log.main};

        samtools sort
          -m 1G -@ {threads}
          -o {output.bam}
          {params.star_tmp_dir}/Aligned.out.bam;

        samtools index {output.bam};

        rm -rf {params.star_tmp_dir}
        ''')


rule sort_bam_by_name:
    """
    Name-sort plate-based STAR BAMs to standardise read order across haplotypes.

    Converts coordinate-sorted STAR BAMs into name-sorted BAMs so that
    read-id ordering is consistent across haplotype-specific alignments. This is required
    for deterministic merging and subsequent per-read haplotype selection.

    Fixmate is also run to update mate information for paired-ended snATAC-seq data.

    This is the plate-based analogue of the align_droplet sort_bam_by_name rule, but each sample_name
    corresponds to one individual/barcode, rather than the whole dataset.

    Parameters:
      Tool behaviour controlled by samtools options in the rule.

    Inputs:
      - coordinate-sorted haplotype BAM and index (results/aligned_data/single_barcodes/haploid)

    Outputs:
      - temporary name-sorted haplotype BAM (results/aligned_data/single_barcodes/haploid)
    """
    input:
        bam=results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam'),
        bai=results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam.bai')
    output:
        bam=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'))
    threads: 4
    resources:
        mem_mb=20_000,
        threads_per_cmd=lambda wc, threads: threads // 2
    conda:
        get_conda_env('htslib')
    shell:
        format_command('''
        samtools sort
          -T ${{TMPDIR}}/{wildcards.sample_name}.{wildcards.qry}
          -n -@ {resources.threads_per_cmd}
          {input.bam} |
        samtools fixmate -@ {resources.threads_per_cmd} -m - {output.bam};
        ''')


def get_merge_haps_input(wc):
    '''Expand all haplotypes that have been aligned to for a sample/barcode to create merge input'''
    dataset_name = SAMPLE_NAME_DATASET_MAPPING[wc.sample_name]
    dataset = config['datasets'][dataset_name]
    qrys = set()
    for geno in dataset['genotypes'].values():
        for qry in geno['founder_haplotypes'].values():
            qrys.add(qry)
    return expand(
        results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'),
        sample_name=wc.sample_name, qry=qrys
    )


rule merge_name_sorted_bams:
    """
    Merge per-haplotype name-sorted for each individual/barcode BAMs into a single name-sorted BAM.

    Combines all haplotype-specific alignments for a dataset into one BAM while
    preserving name-sorted order, producing the input required for haplotype
    collapsing / best-alignment selection.

    This is the plate-based analogue of the align_droplet merge_name_sorted_bams rule, but each sample_name
    corresponds to one individual/barcode, rather than the whole dataset.

    Inputs:
      - name-sorted BAMs for each haplotype aligned for this individual (results/aligned_data/single_barcodes/haploid)

    Outputs:
      - temporary merged name-sorted BAM for the dataset (results/aligned_data/single_barcodes)
    """
    input:
        bams=get_merge_haps_input
    output:
        bam=temp(results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')),
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1024,
    conda:
        get_conda_env('htslib')
    shell:
        format_command('''
        samtools merge -@ {threads}
          -n {output.bam}
          {input.bams};
        ''')


rule collapse_alignments:
    """
    Collapse haplotype-specific alignments into a single best alignment per read.

    Uses the coelsch collapse_ha_specific_alns.py script to choose the best haplotype
    alignment(s) for each read from the merged name-sorted BAM, then sorts and indexes
    the resulting BAM. Adds/updates a haplotype tag [ha] indicating which haplotype(s) were
    selected as best for each read.

    This is the plate-based analogue of the align_droplet collapse_alignments rule, but each sample_name
    corresponds to one individual/barcode, rather than the whole dataset.

    Inputs:
      - merged name-sorted BAM containing alignments to all haplotypes (results/aligned_data/single_barcodes)

    Outputs:
      - final coordinate-sorted haplotype-labelled BAM and index for the individual/barcode (results/aligned_data/single_barcodes)
    """
    input:
        bam=results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')
    output:
        bam=results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        bai=results('aligned_data/single_barcodes/{sample_name}.sorted.bam.bai')
    resources:
        mem_mb=5_000,
    threads: 4
    conda:
        get_conda_env('coelsch')
    shell:
        format_command('''
        collapse_ha_specific_alns.py
          -o {output.bam}.unsorted.bam
          {input.bam};

        samtools addreplacerg
          -r "ID:{wildcards.sample_name}"
          {output.bam}.unsorted.bam |
        samtools sort -@ {threads}
          -T ${{TMPDIR}}/{wildcards.sample_name}
          -o {output.bam} - ;

        samtools index {output.bam};

        rm {output.bam}.unsorted.bam
        ''')


def merge_samples_input(wc):
    '''input for sample/barcode-wise merging'''
    sample_names = DATASET_SAMPLE_NAME_MAPPING[wc.dataset_name]
    return expand(
        results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        sample_name=sample_names
    )


rule list_bam:
    """
    Write a per-dataset BAM list file for sample-wise merging.

    Collects the haplotype-labelled per-individual/barcode BAMs belonging to a dataset and writes them
    to a newline-delimited list file, used as input to `samtools merge -b`.

    Inputs:
      - haplotype-labelled per-individual/barcode BAMs for the dataset (results/aligned_data/single_barcodes)

    Outputs:
      - temporary text file listing BAM paths for the dataset (results)
    """
    input:
        bams=merge_samples_input
    output:
        txt=temp(results("{dataset_name}.list"))
    run:
        with open(output.txt, 'w') as f:
            f.write('\n'.join(input.bams))


rule merge_bams:
    """
    Merge haplotype-labelled per-individual/barcode BAMs into a single dataset-level BAM.

    Merges all per-individual/barcode collapsed BAMs for a dataset into one coordinate-sorted
    BAM and creates its index. This produces the final dataset-level alignment output.

    Inputs:
      - newline-delimited list of per-sample BAMs for the dataset (results)

    Outputs:
      - final merged dataset BAM and index (results/aligned_data)
    """
    input:
        results("{dataset_name}.list")
    output:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai'),
    conda:
        get_conda_env('htslib')
    threads: 24
    shell:
        format_command('''
        ulimit -n 5000;

        samtools merge -@ {threads} -b {input} -o {output.bam};

        samtools index {output.bam};
        ''')
