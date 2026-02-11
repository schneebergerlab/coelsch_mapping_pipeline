rule mappability:
    ''''''
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
    '''align two chromosome-scale genome assemblies with minimap2 for analysis with syri'''
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
    '''use syri to convert a wga of two genome assemblies into a  vcf file of synteny, as well as snps/indels'''
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
    '''builds a config file for msyd'''
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
    '''
    runs msyd to identify "core" synteny between two or more assemblies, 
    outputs a vcf of SNPs/indels in core regions
    '''
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
    Use bedtools to create a blacklist of regions which are non-syntenic.
    Since msyd outputs syntenic regions that can overlap non-syntenic ones, this can help remove noisy markers
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
    '''input for vcf filtering, can be either msyd output created from wga or a user defined vcf file'''
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
    filter and reformat a vcf file (msyd output or user defined) to produce a VCF
    for a single query vs reference comparison suitable for STAR consensus.
    Removes softmasked records, indels > max_indel_size, and blacklisted regions
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
