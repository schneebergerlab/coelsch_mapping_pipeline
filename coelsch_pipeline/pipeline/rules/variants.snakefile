rule mappability:
    '''
    Compute genome mappability for each input genome with GenMap and soft-mask low-mappability regions.

    Uses GenMap to build a k-mer mappability bedGraph for a genotype, then:
      - selects positions with mappability < min_mappability,
      - expands intervals by kmer_size,
      - merges nearby intervals within mask_gap_size,
      - soft-masks those intervals in the FASTA and indexes it.

    Parameters:
      - Defined in config section variants: mappability

    Inputs:
      - genotype fasta file plus index (annotations)

    Outputs:
      - genotype bedgraph of mappability (annotations)
      - genotype softmasked fasta (annotations)
    '''
    input:
        fasta=lambda wc: get_fasta(wc.geno),
        fai=lambda wc: get_fasta(wc.geno) + '.fai',
    output:
        bedgraph=annotation('{geno}.mappability.bedgraph'),
        fasta=annotation('{geno}.softmasked.fa'),
        fai=annotation('{geno}.softmasked.fa.fai'),
    params:
        kmer_size=config['variants']['mappability']['kmer_size'],
        edit_dist=config['variants']['mappability']['edit_dist'],
        min_mappability=config['variants']['mappability']['min_mappability'],
        mask_gap_size=config['variants']['mappability']['mask_gap_size'],
        bedgraph_prefix=lambda wc: annotation(f'{wc.geno}.mappability')
    threads: 16
    conda:
        get_conda_env('genmap')
    shell:
        format_command('''
        genmap index -F {input.fasta} -I {input.fasta}.genmap_idx;

        genmap map -bg -T {threads}
          -K {params.kmer_size}
          -E {params.edit_dist}
          -I {input.fasta}.genmap_idx
          -O {params.bedgraph_prefix};

        awk '$4 < {params.min_mappability}' {output.bedgraph} |
        bedtools slop -l 0 -r {params.kmer_size}
          -g {input.fai}
          -i stdin |
        bedtools merge -d {params.mask_gap_size}
          -i stdin |
        bedtools maskfasta -soft
          -fi {input.fasta}
          -bed stdin
          -fo {output.fasta};

        samtools faidx {output.fasta};
        rm -rf {input.fasta}.genmap_idx;
        ''')


rule minimap2_wga:
    """
    Whole-genome alignment of two chromosome-scale assemblies using minimap2.

    Produces a coordinate-sorted BAM (with index) suitable as input for SyRI.

    Parameters:
      Defined in config section variants: minimap2.

    Inputs:
      - reference genome softmasked fasta (annotations)
      - query softmasked fasta (annotations)

    Outputs:
      - whole genome alignment of query against reference (annotations/wga)
    """
    input:
        ref=annotation('{ref}.softmasked.fa'),
        qry=annotation('{qry}.softmasked.fa'),
    output:
        annotation('wga/{ref}.{qry}.bam')
    threads:
        12
    resources:
        mem_mb=lambda wc, threads: 1024 * threads,
    params:
        preset=config['variants']['minimap2']['preset'],
        zdrop=config['variants']['minimap2']['zdrop'],
    conda:
        get_conda_env('minimap2')
    shell:
        format_command('''
        minimap2 --eqx -t {threads} -N 100
          -ax {params.preset}
          -z{params.zdrop}
          {input.ref} {input.qry} |
        samtools view -bS |
        samtools sort -o - - > {output};

        samtools index {output}
        ''')


rule run_syri:
    """
    Call synteny and small variants between two assemblies with SyRI.

    Consumes a minimap2 WGA BAM and the two softmasked assemblies, producing:
      - SyRI structural/synteny calls (.syri.out)
      - a VCF summarising synteny blocks and SNP/indel calls

    Parameters:
      Defined in config section variants: syri.

    Inputs:
      - reference genome softmasked fasta (annotations)
      - query softmasked fasta (annotations)
      - whole genome alignment of query against reference (annotations/wga)

    Outputs:
      - syri synteny calls (annotations/vcf/syri)
      - syri vcf output (annotations/vcf/syri)
    """
    input:
        ref=annotation('{ref}.softmasked.fa'),
        qry=annotation('{qry}.softmasked.fa'),
        bam=annotation('wga/{ref}.{qry}.bam')
    output:
        vcf=annotation('vcf/syri/{ref}.{qry}.syri.vcf'),
        syri=annotation('vcf/syri/{ref}.{qry}.syri.out'),
    resources:
        mem_mb=5_000,
    conda:
        get_conda_env('msyd')
    params:
        filter_alns='' if config['variants']['syri']['use_low_qual_filters'] else '-f',
        out_dir=annotation('vcf/syri')
    shell:
        format_command('''
        syri -F B --hdrseq
          {params.filter_alns}
          --dir {params.out_dir}
          --prefix {wildcards.ref}.{wildcards.qry}.
          -c {input.bam}
          -q {input.qry}
          -r {input.ref}
          --samplename {wildcards.qry};
        ''')


def get_msyd_input(wc):
    '''full input for msyd, requires genomes, wga bams, and syri output as vcf/syri.out files'''
    ref, *qry_names = wc.geno_group.split('_')
    return {
        'bams': expand(annotation('wga/{ref}.{qry}.bam'), ref=ref, qry=qry_names),
        'syri': expand(annotation('vcf/syri/{ref}.{qry}.syri.out'), ref=ref, qry=qry_names),
        'vcf': expand(annotation('vcf/syri/{ref}.{qry}.syri.vcf'), ref=ref, qry=qry_names),
        'qry_fasta': expand(annotation('{qry}.softmasked.fa'), qry=qry_names),
        'ref_fasta': annotation(f'{ref}.softmasked.fa')
    }


rule msyd_input:
    """
    Create an msyd config TSV for a multi-assembly comparison.

    Writes vcf/msyd/{geno_group}.msyd_config.tsv with one row per query assembly,
    pointing to the WGA BAM, SyRI output, SyRI VCF, and the query genome FASTA.

    Inputs:
      - WGA BAMs for each reference–query comparison (annotations/wga)
      - SyRI structural outputs and VCFs (annotations/vcf/syri)
      - query and reference genome FASTAs (annotations)

    Outputs:
      - msyd configuration TSV (annotations/vcf/msyd)
    """
    input:    
        unpack(get_msyd_input)
    output:
        cfg=temp(annotation('vcf/msyd/{geno_group}.msyd_config.tsv'))
    params:
        annot=config['annotation_dir']
    run:
        ref, *qry_names = wildcards.geno_group.split('_')
        with open(output.cfg, 'w') as f:
            f.write('#name\taln\tsyri\tvcf\tgenome\n')
            for qry in qry_names:
                f.write(
                    f'{qry}\t'
                    f'{params.annot}/wga/{ref}.{qry}.bam\t'
                    f'{params.annot}/vcf/syri/{ref}.{qry}.syri.out\t'
                    f'{params.annot}/vcf/syri/{ref}.{qry}.syri.vcf\t'
                    f'{get_fasta(qry)}\n'
                )


rule run_msyd:
    """
    Identify core synteny across multiple assemblies with msyd and emit SNP/indel VCF.

    Runs `msyd call --core --impute` and post-processes the output VCF to:
      - drop CORESYN ALT records,
      - remove msyd-specific FORMAT/INFO fields,
      - sort, normalise multiallelics, and set missing GTs to reference.

    Inputs:
      - msyd configuration TSV (annotations/vcf/msyd)
      - WGA BAMs for each reference–query comparison (annotations/wga)
      - SyRI structural outputs and VCFs (annotations/vcf/syri)
      - query and reference genome FASTAs (annotations)

    Outputs:
      - msyd pan-synteny feature file (pff) (annotations/vcf/msyd)
      - msyd SNP/indel vcf for core regions (annotations/vcf/msyd)
    """
    input:
        unpack(get_msyd_input),
        cfg=annotation('vcf/msyd/{geno_group}.msyd_config.tsv'),
    output:
        pff=annotation(r'vcf/msyd/{geno_group,\w+}.pansyn.pff'),
        vcf=annotation(r'vcf/msyd/{geno_group,\w+}.vcf'),
    conda:
        get_conda_env('msyd')
    threads: lambda wc, input: min((len(input.bams) + 1) // 2, 12)
    shell:
        format_command('''
        msyd call -c {threads} --core
          --impute
          -i {input.cfg}
          -r {input.ref_fasta}
          -o {output.pff}
          -m {output.vcf}.tmp.vcf;

        bcftools annotate
          --exclude 'ALT ~ "CORESYN"'
          --remove "FORMAT/CHR,FORMAT/START,FORMAT/END,INFO/PID"
          {output.vcf}.tmp.vcf |
        grep -v "^##ALT" |
        grep -v "^##INFO" |
        bcftools sort |
        bcftools norm --multiallelics "+any" |
        bcftools filter -S0 -e 'GT=="."' > {output.vcf};

        rm {output.vcf}.tmp.vcf;
        ''')


def blacklist_input(wc):
    '''input for the blacklist command - the original bam files for minimap2 wga'''
    ref, *qry_names = wc.geno_group.split('_')
    return {
        'bams': expand(annotation('wga/{ref}.{qry}.bam'), ref=ref, qry=qry_names),
        'fai': annotation(f'{ref}.softmasked.fa.fai')
    }


rule blacklist_nonsyntenic_overlapping:
    '''
    Create a blacklist of ambiguous or overlapping non-syntenic regions.

    Regions showing multiple overlapping alignments across assemblies are
    detected using alignment coverage and merged into consolidated intervals
    that can be excluded from downstream variant analyses.

    This overcomes the issue that syri/msyd syntenic alignmets can overlap 
    with non-syntenic alignments.

    Parameters:
      mask_gap_size defined in config section variants: mappability.

    Inputs:
      - WGA BAMs for each reference–query comparison (annotations/wga)
      - reference genome fasta index (annotations)

    Outputs:
      - merged blacklist BED of problematic regions (annotations/bed)
    '''
    input:
        unpack(blacklist_input)
    output:
        blacklist=annotation('bed/{geno_group}.blacklist.bed')
    conda:
        get_conda_env('genmap')
    params:
        mask_gap_size=config['variants']['mappability']['mask_gap_size'],
    shell:
        format_command('''
        for bam in {input.bams}; do
          bedtools bamtobed -i $bam |
          bedtools genomecov -g {input.fai} -bg -i stdin |
          awk '$4 > 1';
        done >> {output.blacklist}.tmp.bed;

        sort -k1,1 -k2,2n {output.blacklist}.tmp.bed |
        bedtools merge -d {params.mask_gap_size} -i stdin > {output.blacklist};

        rm {output.blacklist}.tmp.bed;
        ''')



def filter_snps_vcf_input(wc):
    '''
    input for vcf filtering, can be either msyd output created from wga or 
    a user defined vcf file
    '''
    vcf_fns = config['annotations']['vcf_fns']
    ref, *qry_names = wc.geno_group.split('_')
    if ref in vcf_fns and vcf_fns[ref] is not None:
        return {
            'vcf': ancient(annotation(vcf_fns[ref])),
            'fasta': annotation(f'{ref}.softmasked.fa')
        }
    else:
        return {
            'vcf': annotation('vcf/msyd/{geno_group}.vcf'),
            'blacklist': annotation('bed/{geno_group}.blacklist.bed'),
            'fasta': annotation(f'{ref}.softmasked.fa')
        }
    


rule filter_msyd_snps_for_star_consensus:
    '''
    Filter SNP/indel variants for use in STAR genome consensus generation.

    Extracts variants for a single query sample, removes softmasked or
    blacklisted regions, filters large indels, retains standard nucleotide
    alleles only, and normalises the resulting VCF.

    Parameters:
      max_indel_size defined in config section variants: star_consensus.

    Inputs:
      - msyd or user-specified VCF of structural and sequence variants (annotations)
      - optional blacklist BED of problematic regions (annotations/bed)
      - reference genome softmasked fasta (annotations)

    Outputs:
      - filtered per-genotype SNP/indel VCF suitable for STAR consensus (annotations/vcf/star_consensus)
   
    '''
    input:
        unpack(filter_snps_vcf_input)
    output:
        vcf=annotation('vcf/star_consensus/msyd/{geno_group}.{qry}.vcf'),
    conda:
        get_conda_env('htslib')
    params:
        max_indel_size=config['variants']['star_consensus']['max_indel_size'],
        blacklist_flag=lambda wc, input: f'-T "^{input.blacklist}"' if hasattr(input, 'blacklist') else ''
    shell:
        format_command(r'''
        bcftools view 
          -s {wildcards.qry}
          {input.vcf} |
        bcftools view
          -i 'GT=="alt"'
          {params.blacklist_flag} |
        bcftools view -G
          -e "STRLEN(REF)>{params.max_indel_size} || STRLEN(ALT)>{params.max_indel_size}" |
        awk '$0 ~ /^#/ || ($4 ~ /^[ACGT]+$/ && $5 ~ /^[ACGT]+$/)' |
        bcftools norm --remove-duplicates 
          -f {input.fasta}
        > {output.vcf};
        ''')
