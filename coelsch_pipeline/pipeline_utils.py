import re
import os
from snakemake.io import ancient, expand


def format_command(cmd):
    cmd = [line.rstrip() for line in cmd.splitlines()]
    # replace trailing ";" with ";\n" for readability with --printshellcmds
    formatted = []
    for line in cmd:
        if not line.strip():
            continue
        elif re.search(';$', line):
            formatted.append(f'{line}\n\n')
        else:
            formatted.append(f'{line} \\\n')
    return ''.join(formatted).strip('\n').strip('\\')


def setup_tmpdir(prefix="{rule}"):
    return r'''
        tmpdir=$(mktemp -d "${{TMPDIR:-/tmp}}/''' + prefix + r'''.{jobid}.XXXXXX");
        trap 'rm -rf -- "$tmpdir"' EXIT;
    '''


def resolve_path(basedir, path):
    if isinstance(path, list):
        return [os.path.join(basedir,  p) for p in path]
    return os.path.join(basedir,  path)


def make_path_getter(config_key):
    def factory(config):
        base = config
        for part in config_key.split('/'):
            base = base[part]
        def _wrapped(path):
            return resolve_path(base, path)
        return _wrapped
    return factory


annotations_getter = make_path_getter("annotation_dir")
raw_data_getter = make_path_getter("raw_data_dir")
results_getter = make_path_getter("results_dir")


def make_annotation_getter(config_key):
    def factory(config):
        base = config
        for part in config_key.split('/'):
            base = base[part]
        annotations = annotations_getter(config)
        def _wrapped(annot_key):
            return ancient(annotations(base[annot_key]))
        return _wrapped
    return factory


fasta_getter = make_annotation_getter("annotations/fasta_fns")
gtf_getter = make_annotation_getter("annotations/gtf_fns")
barcode_whitelist_getter = make_annotation_getter("annotations/barcode_whitelist_fns")


def conda_env_getter(config):
    def _wrapped(env_name):
        user_supplied_env = config['conda_envs'].get(env_name)
        if user_supplied_env is None:
            # use default prespecified environment
            return f'../env_yamls/{env_name}.yaml'
        if os.path.splitext(user_supplied_env)[1] in ('.yaml', '.yml'):
            return os.path.abspath(user_supplied_env)
        return user_supplied_env
    return _wrapped