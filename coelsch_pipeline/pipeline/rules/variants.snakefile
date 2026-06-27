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


rule get_whitelist_regions:
    '''find suitable syntentic and inverted regions that can be used for colinear marker analysis'''
    input:
        vcf=annotation('vcf/syri/{ref}.{qry}.syri.vcf'),
        fai=annotation('{ref}.softmasked.fa.fai'),
    output:
        bed=temp(annotation('bed/{ref}.{qry}.syri.bed')),
        blacklist=temp(annotation('bed/{ref}.{qry}.blacklist.bed')),
        whitelist=annotation('bed/{ref}.{qry}.whitelist.bed'),
    resources:
        mem_mb=5_000,
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
    ref, *qry_names = wc.geno_group.split('_')
    return {
        'beds': expand(annotation('bed/{ref}.{qry}.whitelist.bed'), ref=ref, qry=qry_names),
    }


rule get_full_dataset_whitelist_regions:
    input:
        unpack(full_whitelist_input)
    output:
        bed=annotation('bed/merged/{geno_group}.whitelist.bed'),
    resources:
        mem_mb=5_000,
    params:
        nsamples=lambda wc: len(wc.geno_group.split('_')) - 1
    conda:
        get_conda_env('genmap')
    shell:
        format_command('''
        bedtools multiinter -i {input.beds} |
        awk -v OFS='\t' '
          $4 == {params.nsamples} && (($3 - $2) > 50) {{print $1, $2, $3}}
        ' > {output.bed}
        ''')


def filter_syri_vcf_input(wc):
    ref, *_ = wc.geno_group.split('_')
    return {
        'vcf': expand(annotation('vcf/syri/{ref}.{qry}.syri.vcf'), ref=ref, qry=wc.qry),
        'bed': expand(annotation('bed/merged/{geno_group}.whitelist.bed'), geno_group=wc.geno_group),
        'fasta': annotation(f'{ref}.softmasked.fa')
    }


rule filter_syri_vcf:
    input:
        unpack(filter_syri_vcf_input)
    output:
        vcf=annotation('vcf/star_consensus/syri/{geno_group}.{qry}.vcf')
    params:
        max_indel_size=config['variants']['star_consensus']['max_indel_size'],
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
    ref, *qry_names = wc.geno_group.split('_')
    return {
        'vcf': ancient(annotation(vcf_fns[ref])),
        'fasta': annotation(f'{ref}.softmasked.fa')
    }


rule extract_sample_vcf:
    input:
        unpack(extract_sample_vcf_input)
    output:
        vcf=annotation('vcf/star_consensus/user_defined/{geno_group}.{qry}.vcf'),
    conda:
        get_conda_env('htslib')
    params:
        max_indel_size=config['variants']['star_consensus']['max_indel_size'],
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
