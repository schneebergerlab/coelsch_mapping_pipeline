def geno_group_members(wc):
    return wc.geno_group.split('_')


def geno_group_ref(wc):
    return geno_group_members(wc)[0]


def geno_group_qrys(wc):
    return geno_group_members(wc)[1:]


mappability_resources = rule_resources('mappability')


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
        kmer_size=cfg.variants.mappability.kmer_size,
        edit_dist=cfg.variants.mappability.edit_dist,
        min_mappability=cfg.variants.mappability.min_mappability,
        mask_gap_size=cfg.variants.mappability.mask_gap_size,
        bedgraph_prefix=lambda wc: annotation(f'{wc.geno}.mappability')
    threads: mappability_resources.threads
    resources:
        mem_mb=mappability_resources.mem_mb,
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


minimap2_resources = rule_resources('minimap2')


rule minimap2_wga:
    '''align two chromosome-scale genome assemblies with minimap2 for analysis with syri'''
    input:
        ref=annotation('{ref}.softmasked.fa'),
        qry=annotation('{qry}.softmasked.fa'),
    output:
        annotation('wga/{ref}.{qry}.bam')
    threads: minimap2_resources.threads
    resources:
        mem_mb=minimap2_resources.mem_mb,
    params:
        preset=cfg.variants.minimap2.preset,
        zdrop=cfg.variants.minimap2.zdrop,
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


syri_resources = rule_resources('syri')


rule run_syri:
    '''use syri to convert a wga of two genome assemblies into a  vcf file of synteny, as well as snps/indels'''
    input:
        ref=annotation('{ref}.softmasked.fa'),
        qry=annotation('{qry}.softmasked.fa'),
        bam=annotation('wga/{ref}.{qry}.bam')
    output:
        vcf=annotation('vcf/syri/{ref}.{qry}.syri.vcf'),
        syri=annotation('vcf/syri/{ref}.{qry}.syri.out'),
    threads: syri_resources.threads
    resources:
        mem_mb=syri_resources.mem_mb,
    conda:
        get_conda_env('syri')
    params:
        filter_alns='' if cfg.variants.syri.use_low_qual_filters else '-f',
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


bedtools_resources = rule_resources('bedtools')


rule get_whitelist_regions:
    '''find suitable syntentic and inverted regions that can be used for colinear marker analysis'''
    input:
        vcf=annotation('vcf/syri/{ref}.{qry}.syri.vcf'),
        fai=annotation('{ref}.softmasked.fa.fai'),
    output:
        bed=temp(annotation('bed/{ref}.{qry}.syri.bed')),
        blacklist=temp(annotation('bed/{ref}.{qry}.blacklist.bed')),
        whitelist=annotation('bed/{ref}.{qry}.whitelist.bed'),
    threads: bedtools_resources.threads
    resources:
        mem_mb=bedtools_resources.mem_mb,
    conda:
        get_conda_env('genmap')
    shell:
        format_command('''
        awk -v OFS='\t'
          '/^#/ {{next}} 
           ($5 ~ "^<" && $8 !~ "VarType=SR" && $8 ~ /(^|;)END=[0-9]+(;|$)/) {{
             info=$8; sub(/^.*(^|;)END=/,"",info); sub(/;.*/,"",info); print $1, $2-1, info, $3
           }}' {input.vcf} > {output.bed};

        bedtools genomecov 
          -i {output.bed} -bg -g {input.fai} |
        awk -v OFS='\t' '$4 > 1 {{print $1, $2, $3}}' > {output.blacklist};

        awk -v OFS='\t'
           '$4 ~ "^SYNAL" || $4 ~ "^INVAL"' {output.bed} |
        bedtools subtract 
          -a stdin -b {output.blacklist} > {output.whitelist}
        ''')


def full_whitelist_input(wc):
    ref = geno_group_ref(wc)
    qry_names = geno_group_qrys(wc)
    return {
        'beds': expand(annotation('bed/{ref}.{qry}.whitelist.bed'), ref=ref, qry=qry_names),
    }


rule get_full_dataset_whitelist_regions:
    input:
        unpack(full_whitelist_input)
    output:
        bed=annotation('bed/merged/{geno_group}.whitelist.bed'),
    threads: bedtools_resources.threads
    resources:
        mem_mb=bedtools_resources.mem_mb,
    params:
        n_qrys=lambda wc: len(geno_group_qrys(wc))
    conda:
        get_conda_env('genmap')
    shell:
        format_command('''
        if [ {params.n_qrys} -eq 1 ]; then
          awk -v OFS='\t' '(($3 - $2) > 50) {{print $1, $2, $3}}' {input.beds};
        else
          bedtools multiinter -i {input.beds} |
          awk -v OFS='\t' '
            $4 == {params.n_qrys} && (($3 - $2) > 50) {{print $1, $2, $3}}
            ';
        fi > {output.bed};
        ''')


def filter_syri_vcf_input(wc):
    ref = geno_group_ref(wc)
    return {
        'vcf': annotation(f'vcf/syri/{ref}.{wc.qry}.syri.vcf'),
        'bed': annotation(f'bed/merged/{wc.geno_group}.whitelist.bed'),
        'fasta': annotation(f'{ref}.softmasked.fa')
    }


bcftools_resources = rule_resources('bcftools')


rule filter_syri_vcf:
    input:
        unpack(filter_syri_vcf_input)
    output:
        vcf=annotation('vcf/star_consensus/syri/{geno_group}.{qry}.vcf')
    threads: bcftools_resources.threads
    resources:
        mem_mb=bcftools_resources.mem_mb,
    params:
        max_indel_size=cfg.variants.star_consensus.max_indel_size,
    conda:
        get_conda_env('htslib')
    shell:
        format_command(r'''
        bcftools filter -i 'INFO/VarType=="ShV" && ALT !~ "<"' 
          -T {input.bed} {input.vcf} |
        bcftools view -G
          -e "STRLEN(REF)>{params.max_indel_size} || STRLEN(ALT)>{params.max_indel_size}" |
        awk '$0 ~ /^#/ || ($4 ~ /^[ACGT]+$/ && $5 ~ /^[ACGT]+$/)' |
        bcftools norm --remove-duplicates 
          -f {input.fasta} |
        bcftools sort
        > {output.vcf};
        ''')


def extract_sample_vcf_input(wc):
    vcf_fns = config['annotations']['vcf_fns']
    ref = geno_group_ref(wc)
    return {
        'vcf': ancient(annotation(vcf_fns[ref])),
        'fasta': annotation(f'{ref}.softmasked.fa')
    }


extract_sample_vcf_resources = rule_resources('bcftools')


rule extract_sample_vcf:
    input:
        unpack(extract_sample_vcf_input)
    output:
        vcf=annotation('vcf/star_consensus/user_defined/{geno_group}.{qry}.vcf'),
    conda:
        get_conda_env('htslib')
    threads: bcftools_resources.threads
    resources:
        mem_mb=bcftools_resources.mem_mb,
    params:
        max_indel_size=cfg.variants.star_consensus.max_indel_size,
    conda:
        get_conda_env('htslib')
    shell:
        format_command(r'''
        bcftools view 
          -s {wildcards.qry}
          {input.vcf} |
        bcftools view
          -i 'GT=="alt"' |
        bcftools view -G
          -e "STRLEN(REF)>{params.max_indel_size} || STRLEN(ALT)>{params.max_indel_size}" |
        awk '$0 ~ /^#/ || ($4 ~ /^[ACGT]+$/ && $5 ~ /^[ACGT]+$/)' |
        bcftools norm --remove-duplicates 
          -f {input.fasta} |
        bcftools sort
        > {output.vcf};
        ''')
