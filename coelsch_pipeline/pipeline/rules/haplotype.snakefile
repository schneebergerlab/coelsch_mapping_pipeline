import re


def get_cb_tag_id(wc):
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type in ('takara_dna', 'plate_wgs'):
        return 'RG'
    else:
        return 'CB'


rule remove_pcr_duplicates:
    '''filter out PCR duplicates in single-cell ATAC libraries'''
    input:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai'),
    output:
        bam=temp(results('aligned_data/{dataset_name}.deduped.bam')),
        bai=temp(results('aligned_data/{dataset_name}.deduped.bam.bai')),
    params:
        cb_tag=get_cb_tag_id,
    resources:
        mem_mb=10_000,
    conda:
        get_conda_env('htslib')
    threads: 12
    shell:
        format_command('''
        samtools markdup -r -@ {threads}
          --barcode-tag {params.cb_tag}
          {input.bam} {output.bam};

        samtools index {output.bam}
        ''')


def filter_input(wc):
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type not in ('10x_rna_v4', '10x_rna_v3', 'bd_rna') and \
            config['haplotyping']['preprocessing']['filter_pcr_duplicates']:
        input_file_type = 'deduped'
    else:
        # skips PCR deduplication by taking collapsed output as direct input
        input_file_type = 'sorted'
    return dict(
        bam=results(f'aligned_data/{wc.dataset_name}.{input_file_type}.bam'),
        bai=results(f'aligned_data/{wc.dataset_name}.{input_file_type}.bam.bai'),
    )


def get_filter_tag(wc):
    dataset = config['datasets'][wc.dataset_name]
    all_haplos = set()
    for geno in dataset['genotypes'].values():
        all_haplos.update(list(geno['founder_haplotypes'].values()))
    return ','.join(sorted(all_haplos))


rule filter_informative_reads:
    input:
        unpack(filter_input)
    output:
        bam=results('aligned_data/{dataset_name}.filtered.bam'),
        bai=results('aligned_data/{dataset_name}.filtered.bam.bai'),
        cb_whitelist=results('aligned_data/{dataset_name}.initial_whitelist.txt')
    conda:
        '../env_yamls/htslib.yaml'
    params:
        filter_tag=get_filter_tag,
        cb_tag=get_cb_tag_id,
        min_reads=config['haplotyping']['preprocessing']['min_informative_reads_per_barcode'],
    resources:
        mem_mb=20_000
    threads: 4
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


def get_hap_filter_exprs(wc):
    genos = config['datasets'][wc.dataset_name]['genotypes']
    if len(genos) > 1:
        raise ValueError(
            'Cannot split haplotypes for ploidy "diploid_f1*f1" when there is more than one genotype in dataset'
        )
    haps = genos[list(genos)[0]]['founder_haplotypes']
    if wc.haplo == 'hap1':
        return f'\'[ha] == "{haps["parent1"]}" || [ha] == "{haps["parent2"]}"\''
    elif wc.haplo == 'hap2':
        return f'\'[ha] == "{haps["parent3"]}" || [ha] == "{haps["parent4"]}"\''


rule split_f1_x_f1_haplotypes:
    input:
        bam=results('aligned_data/{dataset_name}.filtered.bam'),
        bai=results('aligned_data/{dataset_name}.filtered.bam.bai'),
    output:
        bam=results('aligned_data/{dataset_name}_{haplo}.filtered.bam'),
        bai=results('aligned_data/{dataset_name}_{haplo}.filtered.bam.bai'),
    conda:
        '../env_yamls/htslib.yaml'
    params:
        filter_exprs=get_hap_filter_exprs
    shell:
        format_command('''
        samtools view -b
          -e {params.filter_exprs}
          {input.bam} > {output.bam};
        samtools index {output.bam};
        ''')


dataset_name_plus_haplo_regex = re.compile(
    r'^(?P<ds>' + '|'.join(re.escape(d) for d in config['datasets']) + r')(?:_(?P<hap>hap[12]))?$'
)

def parse_dataset_name(tok):
    m = dataset_name_plus_haplo_regex.fullmatch(tok)
    if not m:
        raise ValueError(f"Bad token: {tok}")
    ds, hap = m.group('ds'), m.group('hap')
    if hap is None and config['datasets'][ds]['ploidy'] == 'diploid_f1*f1':
        raise ValueError('diploid_f1*f1 requires hap suffix (hap1|hap2)')
    return ds, hap


def get_haplotyping_input(wc):
    dataset_name, haplo = parse_dataset_name(wc.dataset_name_plus_optional_haplo)
    dataset = config['datasets'][dataset_name]
    if dataset['ploidy'] != 'diploid_f1*f1':
        input_ = {
            'bam': results(f'aligned_data/{dataset_name}.filtered.bam'),
            'bai': results(f'aligned_data/{dataset_name}.filtered.bam.bai'),
            'cb_whitelist': results(f'aligned_data/{dataset_name}.initial_whitelist.txt')
        }
    else:
        input_ = {
            'bam': results(f'aligned_data/{dataset_name}_{haplo}.filtered.bam'),
            'bai': results(f'aligned_data/{dataset_name}_{haplo}.filtered.bam.bai'),
            'cb_whitelist': results(f'aligned_data/{dataset_name}.initial_whitelist.txt')
        }
    if dataset['ploidy'] == 'haploid_f1*f1':
        genos = dataset['genotypes']
        hap_fns = genos[list(genos)[0]]['recombinant_haplotypes']
        input_['hap1_json_fn'] = annotation(hap_fns['parent12'])
        input_['hap2_json_fn'] = annotation(hap_fns['parent34'])
    ref = config['datasets'][dataset_name]['reference_genotype']
    mask_bed_fn = config['annotations']['mask_bed_fns'][ref]
    if mask_bed_fn is not None:
        input_['bed'] = annotation(mask_bed_fn)
    return input_


def get_genotyping_params(wc, input):
    dataset_name, haplo = parse_dataset_name(wc.dataset_name_plus_optional_haplo)
    dataset = config['datasets'][dataset_name]
    ploidy_mode = dataset['ploidy']
    genos = dataset['genotypes']
    if ploidy_mode != 'haploid_f1*f1':
        crosses = []
        if ploidy_mode != 'diploid_f1*f1':
            for geno in dataset['genotypes'].values():
                crosses.append(':'.join(sorted(geno['founder_haplotypes'].values())))
        else:
            if len(genos) > 1:
                raise ValueError(
                    'Cannot split haplotypes for ploidy "diploid_f1*f1" when there is more than one genotype in dataset'
                )
            haps = genos[list(genos)[0]]['founder_haplotypes']
            if haplo == 'hap1':
                crosses.append(":".join(sorted([haps["parent1"], haps["parent2"]])))
            elif haplo == 'hap2':
                crosses.append(":".join(sorted([haps["parent3"], haps["parent4"]])))
        if len(crosses) > 1:
            genotype_params = '--genotype --clean-by-genotype -X '
        else:
            genotype_params = '-X '
        return genotype_params + ','.join(crosses)
    else:
        # special haploid_f1*f1 mode, provide haplotype jsons for genotyping
        parent12_hap_fn = input.hap1_json_fn
        parent34_hap_fn = input.hap2_json_fn
        return f'--genotype --clean-by-genotype --recombinant-parent-jsons {parent12_hap_fn} {parent34_hap_fn}'


def get_tech_specific_params(wc):
    dataset_name, haplo = parse_dataset_name(wc.dataset_name_plus_optional_haplo)
    tech_type = config['datasets'][dataset_name]['technology']
    if tech_type == "10x_atac":
        params = '''
          -x {x_flag}
          -y {ploidy}
          --cb-tag CB
          --hap-tag ha
          --hap-tag-type "multi_haplotype"
          --cb-correction-method "exact"
          --umi-collapse-method "none"
          --no-clean-bg
       '''
    elif tech_type in ("takara_dna", "plate_wgs"):
        params = '''
          -x {x_flag}
          -y {ploidy}
          --cb-tag RG
          --hap-tag ha
          --hap-tag-type "multi_haplotype"
          --cb-correction-method "none"
          --no-validate
          --umi-collapse-method "none"
          --no-clean-bg
          '''
        if tech_type == 'plate_wgs':
            params += '--no-predict-doublets'
    else:
        params = '''
          -x {x_flag}
          -y {ploidy}
          --cb-tag CB
          --umi-tag UB
          --hap-tag ha
          --hap-tag-type "multi_haplotype"
          --cb-correction-method "exact"
          --umi-collapse-method "exact"
       '''

    if tech_type.startswith('10x_rna'):
        x_flag = '10x_rna'
    elif tech_type == 'plate_wgs':
        x_flag = 'wgs'
    else:
        x_flag = tech_type
    ploidy = config['datasets'][dataset_name]['ploidy']
    if ploidy == 'diploid_f1*f1' or ploidy.startswith('haploid'):
        # switch to haploid since haplotypes have already been separated
        ploidy = 'haploid'
    params = params.format(
        x_flag=x_flag,
        ploidy=ploidy
    )
    return format_command(params.lstrip())


rule run_haplotyping:
    input:
        unpack(get_haplotyping_input)
    output:
        markers_init=results('haplotypes/{dataset_name_plus_optional_haplo}.markers_init.json'),
        markers_filt=results('haplotypes/{dataset_name_plus_optional_haplo}.markers_filt.json'),
        preds=results('haplotypes/{dataset_name_plus_optional_haplo}.pred.json'),
        stats=results('haplotypes/{dataset_name_plus_optional_haplo}.pred.stats.tsv'),
    conda:
        get_conda_env('coelsch')
    params:
        bin_size=config['haplotyping']['coelsch']['genomic_bin_size'],
        max_frac_bg=config['haplotyping']['preprocessing'].get('max_fraction_background', 0.25),
        min_geno_prob=config['haplotyping']['preprocessing'].get('min_genotyping_probability', 0.9),
        max_geno_error=config['haplotyping']['preprocessing'].get('max_genotyping_error_rate', 0.25),
        max_marker_imbalance=config['haplotyping']['preprocessing'].get('max_marker_imbalance', 0.75),
        rfactor=config['haplotyping']['coelsch']['segment_size'],
        term_rfactor=config['haplotyping']['coelsch']['terminal_segment_size'],
        cm_per_mb=config['haplotyping']['coelsch']['cm_per_mb'],
        mask_bed_flag=lambda wc, input: f'-m {input.bed}' if hasattr(input, 'bed') else '',
        genotyping_params=get_genotyping_params,
        tech_specific_params=get_tech_specific_params,
        min_reads_per_cb=config['haplotyping']['preprocessing']['min_informative_reads_per_barcode'],
        min_reads_per_chrom=config['haplotyping']['preprocessing']['min_informative_reads_per_chrom'],
        output_prefix=lambda wc: results(f'haplotypes/{wc.dataset_name_plus_optional_haplo}')
    threads: 50
    resources:
        mem_mb=150_000
    shell:
        format_command('''
        export OPENBLAS_NUM_THREADS=1;
        export OMP_NUM_THREADS=1;
        export MKL_NUM_THREADS=1;

        coelsch bam2pred -v debug -p {threads}
          --cb-whitelist-fn {input.cb_whitelist}
          {params.mask_bed_flag}
          -N {params.bin_size}
          -R {params.rfactor}
          -t {params.term_rfactor}
          -C {params.cm_per_mb}
          {params.tech_specific_params}
          --min-markers-per-cb {params.min_reads_per_cb}
          --min-markers-per-chrom {params.min_reads_per_chrom}
          --min-genotyping-prob {params.min_geno_prob}
          --max-genotyping-error {params.max_geno_error}
          --max-frac-bg {params.max_frac_bg}
          --max-marker-imbalance {params.max_marker_imbalance}
          --batch-size 128
          {params.genotyping_params}
          -o {params.output_prefix}
          {input.bam}
        ''')


rule haplotyping_report:
    input:
        markers_init=results('haplotypes/{dataset_name_plus_optional_haplo}.markers_init.json'),
        markers_filt=results('haplotypes/{dataset_name_plus_optional_haplo}.markers_filt.json'),
        preds=results('haplotypes/{dataset_name_plus_optional_haplo}.pred.json'),
        stats=results('haplotypes/{dataset_name_plus_optional_haplo}.pred.stats.tsv'),
    log:
        notebook=results('analysis/{dataset_name_plus_optional_haplo}.haplotyping_report.py.ipynb')
    resources:
        mem_mb=50_000
    conda:
        get_conda_env('coelsch')
    notebook:
        '../notebook_templates/haplotyping_report.py.ipynb'