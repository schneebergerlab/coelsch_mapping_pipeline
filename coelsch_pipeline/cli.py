import os
import shutil
import textwrap
import sys
from pathlib import Path
import yaml

import click
from jinja2 import Environment, FileSystemLoader
from coelsch_pipeline.pipeline_utils import (
    parse_crossing_combinations,
    founder_haplotypes_from_crossing_combinations,
)


BASE_DIR = Path(__file__).parent
PIPELINE_DIR = BASE_DIR / "pipeline"
CONFIG_DIR = PIPELINE_DIR / "config"
ANNOTATIONS_DIR = PIPELINE_DIR / "annotations"
SNAKEFILE = PIPELINE_DIR / "Snakefile"

CONFIG_VARIABLES = [
    'annotation_dir', 'raw_data_dir', 'results_dir',
    'dataset_name', 'technology', 'sample_unit', 'lifecycle_stage',
    'crossing_combinations', 'reference_name', 'input_fastq_basename',
]

BARCODE_WHITELIST_FNS = {
    '10x_rna': ['3M-3pgex-may-2023.txt'],
    '10x_atac': ['737K-cratac-v1.txt'],
    'bd_rna': ['BD_CLS1.txt', 'BD_CLS2.txt', 'BD_CLS3.txt'],
    'bd_atac': None,
    'takara_dna': None,
    'wgs': None,
    'other': None,
}


PROMPTS = {
    'annotation_dir': {
        'default': 'annotations',
        'help': 'Directory containing reference annotations such as FASTA, GTF, VCF and mask BED files.',
    },
    'raw_data_dir': {
        'default': 'raw_data',
        'help': 'Directory containing input FASTQ files.',
    },
    'results_dir': {
        'default': 'results',
        'help': 'Directory where pipeline outputs will be written.',
    },
    'dataset_name': {
        'default': 'example_dataset',
        'help': 'Name of this dataset in the generated YAML and output paths.',
    },
    'technology': {
        'default': '10x_rna',
        'help': ('Sequencing technology preset passed to coelsch --seq-type. '
                 'Choices: 10x_rna, 10x_atac, bd_rna, bd_atac, takara_dna, wgs, other.'),
    },
    'sample_unit': {
        'default': 'auto',
        'help': ('Unit represented by each barcode/sample. Choices: auto, single_cell, bulk. '
                 'Use bulk for plate or WGS individuals.'),
    },
    'lifecycle_stage': {
        'default': 'gametes',
        'help': 'Biological stage of sampled material. Choices: gametes, progeny.',
    },
    'crossing_combinations': {
        'default': '(col0*ler)',
        'help': ('Genotype expression(s) passed to coelsch --crossing-combinations. '
                 'Use comma-separated expressions for pooled crosses, e.g. (A*B),((C*D)*(E*F)).'),
    },
    'reference_name': {
        'help': 'Reference genotype used for mapping. Must be one of the founder haplotypes in crossing_combinations.',
    },
    'input_fastq_basename': {
        'help': 'Base name of input FASTQ files, before suffixes such as .1.fastq.gz and .2.fastq.gz.',
    },
}

if not SNAKEFILE.exists():
    click.echo(f"Error: Snakefile not found at {SNAKEFILE}", err=True)
    sys.exit(1)


def _raise_template_error(message):
    raise click.ClickException(message)


def _format_founder_fasta_fns(founder_haplotypes):
    return '\n'.join(
        f'    {hap}: "{hap}.fa"' for hap in founder_haplotypes
    )


def _format_crossing_combinations(crossing_combinations):
    return '\n'.join(
        f'      - "{expr}"' for expr in crossing_combinations
    )


def _wrapped_prompt_lines(text, width):
    lines = []
    for line in text.splitlines():
        lines.extend(textwrap.wrap(line, width=width) or [''])
    return lines


def _prompt_config_value(name, default, help_text=None, clear=True):
    lines_printed = 0
    width = max(20, shutil.get_terminal_size((80, 20)).columns)

    if help_text:
        for line in _wrapped_prompt_lines(help_text, width):
            click.secho(line, fg='bright_black')
            lines_printed += 1

    value = click.prompt(
        click.style(name, fg='bright_yellow'),
        default=default,
    )
    prompt_width = len(f'{name} [{default}]: ')
    lines_printed += max(1, (prompt_width + width - 1) // width)

    if clear and sys.stdout.isatty():
        for _ in range(lines_printed):
            click.echo('\033[F\033[K', nl=False)

    return value


def _echo_config_summary(context):
    click.secho('Config values:', fg='green')
    for key in (
        'annotation_dir', 'raw_data_dir', 'results_dir', 'dataset_name',
        'technology', 'sample_unit', 'lifecycle_stage', 'reference_name',
        'input_fastq_basename',
    ):
        click.echo(f'  {key}: {context[key]}')

    click.echo('  crossing_combinations:')
    for expr in context['crossing_combinations']:
        click.echo(f'    - {expr}')
    click.echo('  founder_haplotypes: ' + ', '.join(context['founder_haplotypes']))


@click.group()
def cli():
    """CLI for the pipeline."""
    pass


@cli.command("init")
@click.argument("destination", type=click.Path(), default='coelsch_mapping_config.yaml')
@click.option("-f", "--force", is_flag=True, default=False)
@click.option(
    '--clear-prompts/--keep-prompts',
    default=True,
    help='Clear each explanatory prompt after it is answered in interactive terminals.',
)
def init_config(destination, force, clear_prompts):
    """Prompt for missing values and render config.yaml from template."""
    dest = Path(destination)
    if dest.exists() and not force:
        click.secho(f"{dest} already exists. Refusing to overwrite.", err=True, fg="red")
        sys.exit(1)

    env = Environment(loader=FileSystemLoader(CONFIG_DIR))
    env.globals.update(
        parse_founder_haplotypes=founder_haplotypes_from_crossing_combinations,
        format_founder_fasta_fns=_format_founder_fasta_fns,
        format_crossing_combinations=_format_crossing_combinations,
        raise_error=_raise_template_error,
    )
    template_name = "default_config.yaml.j2"
    context = {}
    for var in CONFIG_VARIABLES:
        prompt = PROMPTS[var]
        if var == 'reference_name':
            founders = founder_haplotypes_from_crossing_combinations(context['crossing_combinations'])
            default = founders[0]
            help_text = prompt['help'] + f' Available founders: {", ".join(founders)}.'
        elif var == 'input_fastq_basename':
            default = context.get('dataset_name', PROMPTS['dataset_name']['default'])
            help_text = prompt['help']
        else:
            default = prompt['default']
            help_text = prompt['help']

        value = _prompt_config_value(
            var,
            default=default,
            help_text=help_text,
            clear=clear_prompts,
        )

        if var == 'crossing_combinations':
            value = parse_crossing_combinations(value)
            founders = founder_haplotypes_from_crossing_combinations(value)
            if not founders:
                raise click.ClickException(
                    'At least one valid crossing combination is required'
                )
        elif var == 'reference_name':
            founders = founder_haplotypes_from_crossing_combinations(context['crossing_combinations'])
            if value not in founders:
                raise click.ClickException(
                    f'reference_name "{value}" is not one of the founder '
                    f'haplotypes in crossing_combinations: {", ".join(founders)}'
                )
        context[var] = value

    context['founder_haplotypes'] = founder_haplotypes_from_crossing_combinations(
        context['crossing_combinations']
    )

    tech_type = context.get('technology', '10x_rna')
    context['barcode_whitelist_fns'] = str(BARCODE_WHITELIST_FNS.get(tech_type, None))

    _echo_config_summary(context)

    # Render and write output
    template = env.get_template(template_name)
    rendered = template.render(context)
    dest.write_text(rendered)
    click.secho(f"Config written to {dest}", fg="green")
    config = yaml.safe_load(rendered)
    for directory in ['annotation_dir', 'raw_data_dir', 'results_dir']:
        dir_name = config[directory]
        if not os.path.exists(dir_name):
            click.secho(f'Creating {directory}: {dir_name}', fg="bright_yellow")
            os.mkdir(dir_name)


@cli.command("run", context_settings={"ignore_unknown_options": True, "allow_extra_args": True})
@click.option("--configfile", "-c", type=click.Path(exists=True), default="coelsch_mapping_config.yaml",
              help="Path to config file")
@click.pass_context
def run_pipeline(ctx, configfile):
    """Run the pipeline."""

    args = [
        "--snakefile", str(SNAKEFILE.resolve()),
        "--configfile", str(configfile),
    ]

    # Add any unknown extra CLI args
    args += ctx.args

    click.secho(f"Running: snakemake {' '.join(args)}", fg="bright_yellow")

    try:
        # Snakemake ≥7.x
        from snakemake.cli import main as snakemake
    except ImportError:
        # Snakemake ≤6.x
        from snakemake import main as snakemake

    snakemake(args)

if __name__ == "__main__":
    cli()
