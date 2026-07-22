def get_cb_tag_id(wc):
    if dataset(wc).technology in ('takara_dna', 'wgs'):
        return 'RG'
    return 'CB'


samtools_markdup_resources = rule_resources('samtools/dedup')


rule remove_pcr_duplicates:
    '''Remove PCR duplicates from non-UMI BAMs when configured.'''
    input:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai'),
    output:
        bam=temp(results('aligned_data/{dataset_name}.deduped.bam')),
        bai=temp(results('aligned_data/{dataset_name}.deduped.bam.bai')),
    params:
        cb_tag=get_cb_tag_id,
    threads: samtools_markdup_resources.threads
    resources:
        mem_mb=samtools_markdup_resources.mem_mb,
    conda:
        get_conda_env('htslib')
    shell:
        format_command('''
        samtools markdup -r -@ {threads}
          --barcode-tag {params.cb_tag}
          {input.bam} {output.bam};

        samtools index {output.bam}
        ''')


def filter_input(wc):
    tech = dataset(wc).technology
    if tech not in ('10x_rna', 'bd_rna') and \
            cfg.haplotyping.clean.bulk.filter_pcr_duplicates:
        input_file_type = 'deduped'
    else:
        # skips PCR deduplication by taking collapsed output as direct input
        input_file_type = 'sorted'
    return dict(
        bam=results(f'aligned_data/{wc.dataset_name}.{input_file_type}.bam'),
        bai=results(f'aligned_data/{wc.dataset_name}.{input_file_type}.bam.bai'),
    )


def get_filter_tag(wc):
    haps = founder_haplotypes_from_crossing_combinations(
        dataset(wc).crossing_combinations
    )
    return ','.join(sorted(haps))


filter_reads_resources = rule_resources('samtools/view')


rule filter_informative_reads:
    '''Keep informative haplotype-tagged reads and create an initial barcode whitelist.'''
    input:
        unpack(filter_input)
    output:
        bam=results('aligned_data/{dataset_name}.filtered.bam'),
        bai=results('aligned_data/{dataset_name}.filtered.bam.bai'),
        cb_whitelist=results('aligned_data/{dataset_name}.initial_whitelist.txt')
    conda:
        get_conda_env('htslib')
    params:
        filter_tag=get_filter_tag,
        cb_tag=get_cb_tag_id,
        min_reads=cfg.haplotyping.clean.min_markers_per_barcode,
    threads: filter_reads_resources.threads
    resources:
        mem_mb=filter_reads_resources.mem_mb,
    shell:
        format_command('''
        samtools view -b -@ {threads}
          -e '[ha] != "{params.filter_tag}" && [{params.cb_tag}] != "-"'
          {input.bam}
        > {output.bam};

        samtools index {output.bam};

        samtools view -@ {threads} --keep-tag "{params.cb_tag}" {output.bam} |
        awk '{{cb_count[substr($12, 6)]++}}
             END {{for (cb in cb_count)
                 {{if (cb_count[cb] > {params.min_reads})
                 {{print cb}}}}}}'
        > {output.cb_whitelist}
        ''')


def get_load_input(wc):
    ds = dataset(wc)
    input_ = {
        'bam': results(f'aligned_data/{wc.dataset_name}.filtered.bam'),
        'bai': results(f'aligned_data/{wc.dataset_name}.filtered.bam.bai'),
        'cb_whitelist': results(f'aligned_data/{wc.dataset_name}.initial_whitelist.txt'),
    }

    for i, parent_json in enumerate(ds.recombinant_parent_jsons or []):
        input_[f'recombinant_parent_json_{i}'] = annotation(parent_json)

    return input_


def get_clean_input(wc):
    ds = dataset(wc)
    input_ = {
        'markers': results(f'haplotypes/{wc.dataset_name}.markers_init.json'),
    }

    mask_bed_fn = config['annotations']['mask_bed_fns'][ds.reference_genotype]
    if mask_bed_fn is not None:
        input_['bed'] = annotation(mask_bed_fn)

    return input_


def run_genotype_flag(ds):
    recombinant_jsons = ds.recombinant_parent_jsons or []
    pooled_crosses = len(ds.crossing_combinations or []) > 1

    if recombinant_jsons or pooled_crosses:
        return '--genotype'
    return '--no-genotype'


def recombinant_parent_json_flags(ds):
    return ' '.join(
        f'--recombinant-parent-json {annotation(fn)}'
        for fn in (ds.recombinant_parent_jsons or [])
    )


def load_tag_params(ds):
    tech = ds.technology

    if tech in ('wgs', 'takara_dna'):
        return '''
          --cb-tag RG
          --hap-tag ha
          --hap-tag-type multi_haplotype
          --cb-correction-method none
          --umi-collapse-method none
        '''

    if tech == '10x_atac':
        return '''
          --cb-tag CB
          --hap-tag ha
          --hap-tag-type multi_haplotype
          --cb-correction-method exact
          --umi-collapse-method none
        '''

    return '''
      --cb-tag CB
      --umi-tag UB
      --hap-tag ha
      --hap-tag-type multi_haplotype
      --cb-correction-method exact
      --umi-collapse-method exact
    '''


def get_load_params(wc):
    ds = dataset(wc)
    load = cfg.haplotyping.load

    params = [
        f'-x {ds.technology}',
        f'--sample-unit {ds.sample_unit}',
        f'--lifecycle-stage {ds.lifecycle_stage}',
        f'--crossing-strategy {ds.crossing_strategy}',
        f'--genotyping-strategy {ds.genotyping_strategy}',
        '-X "' + ','.join(ds.crossing_combinations) + '"',
        run_genotype_flag(ds),
        recombinant_parent_json_flags(ds),
        f'--min-alignment-score {load.min_alignment_score}',
        optional_flag('--min-mapq', load.min_mapq),
        optional_flag('--exclude-contigs', load.exclude_contigs),
        auto_bool_flag(
            load.validate_barcodes,
            '--validate-barcodes',
            '--no-validate',
        ),
        load_tag_params(ds),
    ]

    return format_params(params, indent=' ' * 10)


def get_clean_params(wc):
    clean = cfg.haplotyping.clean
    single_cell = clean.single_cell

    params = [
        f'--min-markers-per-cb {clean.min_markers_per_barcode}',
        f'--min-markers-per-chrom {clean.min_markers_per_chrom}',
        f'--min-genotyping-prob {clean.min_genotyping_probability}',
        f'--max-genotyping-error {clean.max_genotyping_error_rate}',
        bool_flag(clean.normalise_bins, '--normalise-bins', '--no-normalise-bins'),
        f'--bin-shrinkage-quantile {clean.bin_shrinkage_quantile}',
        auto_bool_flag(
            clean.normalise_depth,
            '--normalise-depth',
            '--no-normalise-depth',
        ),
        optional_flag('--max-bin-count', clean.max_bin_count),
        bool_flag(clean.mask_imbalanced, '--mask-imbalanced', '--no-mask-imbalanced'),
        f'--max-marker-imbalance {clean.max_marker_imbalance}',
        bool_flag(clean.clean_by_genotype, '--clean-by-genotype', '--clean-all'),
        bool_flag(single_cell.clean_background, '--clean-bg', '--no-clean-bg'),
        f'--max-frac-bg {single_cell.max_fraction_background}',
    ]

    return format_params(params, indent=' ' * 10)


def get_predict_params(wc):
    pred = cfg.haplotyping.predict

    params = [
        f'-R {pred.segment_size}',
        f'-t {pred.terminal_segment_size}',
        f'-C {pred.cm_per_mb}',
        f'--interference-half-life {pred.interference_half_life}',
        f'--distribution-type {pred.distribution_type}',
        bool_flag(pred.sample_co_locs, '--sample-co-locs', '--no-sample-cos'),
        f'--n-samples {pred.n_samples}',
        auto_bool_flag(
            pred.predict_doublets,
            '--predict-doublets',
            '--no-predict-doublets',
        ),
        f'--n-doublets {pred.n_doublets}',
        f'--nco-min-prob-change {pred.nco_min_prob_change}',
        f'--batch-size {pred.batch_size}',
    ]

    return format_params(params, indent=' ' * 10)


coelsch_loadbam_resources = rule_resources('coelsch/loadbam')


rule run_loadbam:
    '''Load a haplotype-tagged BAM into an initial MarkerRecords JSON.'''
    input:
        unpack(get_load_input)
    output:
        markers=results('haplotypes/{dataset_name}.markers_init.json')
    conda:
        get_conda_env('coelsch')
    params:
        bin_size=cfg.haplotyping.predict.genomic_bin_size,
        load_params=get_load_params,
    threads: coelsch_loadbam_resources.threads
    resources:
        mem_mb=coelsch_loadbam_resources.mem_mb
    shell:
        format_command('''
        export OPENBLAS_NUM_THREADS=1;
        export OMP_NUM_THREADS=1;
        export MKL_NUM_THREADS=1;

        coelsch loadbam -v debug -p {threads}
          --cb-whitelist-fn {input.cb_whitelist}
          -N {params.bin_size}
          {params.load_params}
          -o {output.markers}
          {input.bam}
        ''')


coelsch_clean_resources = rule_resources('coelsch/clean')


rule run_clean:
    '''Filter, mask and normalise MarkerRecords before prediction.'''
    input:
        unpack(get_clean_input)
    output:
        markers=results('haplotypes/{dataset_name}.markers_filt.json')
    conda:
        get_conda_env('coelsch')
    params:
        bin_size=cfg.haplotyping.predict.genomic_bin_size,
        clean_params=get_clean_params,
        mask_bed_flag=lambda wc, input: f'-m {input.bed}' if hasattr(input, 'bed') else '',
    threads: coelsch_clean_resources.threads
    resources:
        mem_mb=coelsch_clean_resources.mem_mb
    shell:
        format_command('''
        coelsch clean -v debug
          -N {params.bin_size}
          {params.mask_bed_flag}
          {params.clean_params}
          -o {output.markers}
          {input.markers}
        ''')


coelsch_predict_resources = rule_resources('coelsch/predict')


rule run_predict:
    '''Infer haplotype probabilities and crossovers from cleaned MarkerRecords.'''
    input:
        markers=results('haplotypes/{dataset_name}.markers_filt.json')
    output:
        preds=results('haplotypes/{dataset_name}.pred.json'),
        stats=results('haplotypes/{dataset_name}.pred.stats.tsv')
    conda:
        get_conda_env('coelsch')
    params:
        bin_size=cfg.haplotyping.predict.genomic_bin_size,
        predict_params=get_predict_params,
    threads: coelsch_predict_resources.threads
    resources:
        mem_mb=coelsch_predict_resources.mem_mb
    shell:
        format_command('''
        export OPENBLAS_NUM_THREADS=1;
        export OMP_NUM_THREADS=1;
        export MKL_NUM_THREADS=1;

        coelsch predict -v debug -p {threads}
          -N {params.bin_size}
          {params.predict_params}
          -o {output.preds}
          {input.markers}
        ''')


report_resources = rule_resources('haplotyping_report')


rule haplotyping_report:
    '''Render the dataset haplotyping report notebook.'''
    input:
        markers_init=results('haplotypes/{dataset_name}.markers_init.json'),
        markers_filt=results('haplotypes/{dataset_name}.markers_filt.json'),
        preds=results('haplotypes/{dataset_name}.pred.json'),
        stats=results('haplotypes/{dataset_name}.pred.stats.tsv'),
    log:
        notebook=results('analysis/{dataset_name}.haplotyping_report.py.ipynb')
    threads: report_resources.threads
    resources:
        mem_mb=report_resources.mem_mb
    conda:
        get_conda_env('coelsch')
    notebook:
        '../notebook_templates/haplotyping_report.py.ipynb'