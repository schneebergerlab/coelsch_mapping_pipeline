import re
import os
from dataclasses import dataclass
from types import SimpleNamespace


_GENOTYPE_TOKEN_RE = re.compile(r'[^*()[\],:\s]+')
_SEX_LABEL_RE = re.compile(r'\[[^\]]*\]')


def parse_crossing_combinations(value):
    if isinstance(value, str):
        return [token.strip() for token in value.split(',') if token.strip()]
    return list(value)


def founder_haplotypes_from_crossing_combinations(crossing_combinations):
    crossing_combinations = parse_crossing_combinations(crossing_combinations)

    founders = []
    seen = set()
    for expr in crossing_combinations:
        expr = expr.strip()
        if expr.count('(') != expr.count(')'):
            raise ValueError(
                f'Unbalanced parentheses in genotype expression: {expr!r}'
            )

        if ':' in expr and '(' in expr and expr.index(':') < expr.index('('):
            expr = expr.split(':', 1)[1]

        expr = _SEX_LABEL_RE.sub('', expr)
        for hap in _GENOTYPE_TOKEN_RE.findall(expr):
            if hap not in seen:
                founders.append(hap)
                seen.add(hap)
    return founders


def config_namespace(obj, keep_dict_keys=frozenset({'datasets'})):
    def convert(value, key=None):
        if isinstance(value, dict):
            converted = {k: convert(v, k) for k, v in value.items()}
            if key in keep_dict_keys:
                return converted
            return SimpleNamespace(**converted)
        if isinstance(value, list):
            return [convert(v) for v in value]
        return value
    return convert(obj)


def bool_flag(value, true_flag, false_flag):
    if value is True:
        return true_flag
    if value is False:
        return false_flag
    if isinstance(value, str):
        if value.lower() == 'true':
            return true_flag
        if value.lower() == 'false':
            return false_flag
    raise ValueError(f'Expected boolean for {true_flag}/{false_flag}, got {value!r}')


def auto_bool_flag(value, true_flag, false_flag):
    if value is None:
        return ''
    if isinstance(value, str) and value.lower() == 'auto':
        return ''
    return bool_flag(value, true_flag, false_flag)


def optional_flag(flag, value):
    if value is None:
        return ''
    return f'{flag} {value}'


def format_params(params, indent):
    lines = []
    for param in params:
        if not param:
            continue
        lines.extend(
            line.strip()
            for line in str(param).splitlines()
            if line.strip()
        )
    return format_command((f'\n{indent}').join(lines))


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


def setup_tmpdir_getter(cfg):
    runtime = _get_attr_or_key(cfg, 'runtime')
    tmpdir = _get_attr_or_key(runtime, 'tmpdir')

    def _wrapped(prefix="{rule}"):
        base = '${TMPDIR:-/tmp}' if tmpdir is None else tmpdir
        return f"""
        tmpdir=$(mktemp -d "{base}/{prefix}.{{jobid}}.XXXXXX");
        trap 'rm -rf -- "$tmpdir"' EXIT;
    """
    return _wrapped


def resolve_path(basedir, path):
    if isinstance(path, list):
        return [os.path.join(basedir,  p) for p in path]
    return os.path.join(basedir,  path)



@dataclass
class RuleResources:
    threads: int
    mem_mb: int = None
    mem_per_thread: int = None

    def __post_init__(self):
        self.threads = int(self.threads)
        if self.mem_per_thread is not None:
            self.mem_per_thread = int(self.mem_per_thread)
            self.mem_mb = self.threads * self.mem_per_thread
        elif self.mem_mb is not None:
            self.mem_mb = int(self.mem_mb)
        else:
            raise ValueError('Either mem_mb or mem_per_thread must be set')


def _get_attr_or_key(obj, key, default=None):
    if obj is None:
        return default
    if isinstance(obj, dict):
        return obj.get(key, default)
    return getattr(obj, key, default)


def _get_path_or_key(obj, path, default=None):
    current = obj
    for part in path.split('/'):
        current = _get_attr_or_key(current, part)
        if current is None:
            return default
    return current


def _make_rule_resources(conf, default):
    threads = _get_attr_or_key(conf, 'threads', default.threads)
    mem_per_thread = _get_attr_or_key(conf, 'mem_per_thread')
    if mem_per_thread is not None:
        mem_mb = None
    else:
        mem_mb = _get_attr_or_key(conf, 'mem_mb', default.mem_mb)
    return RuleResources(
        threads=threads,
        mem_mb=mem_mb,
        mem_per_thread=mem_per_thread,
    )


def rule_resources_getter(cfg):
    resources = _get_attr_or_key(
        _get_attr_or_key(cfg, 'runtime'),
        'resources',
        {},
    )
    default = _make_rule_resources(
        _get_path_or_key(resources, 'default', {}),
        RuleResources(threads=4, mem_mb=8000),
    )

    def _wrapped(rule_name):
        return _make_rule_resources(
            _get_path_or_key(resources, rule_name, {}),
            default,
        )
    return _wrapped


def dataset_getter(cfg):
    def _wrapped(wc):
        if isinstance(wc, str):
            return cfg.datasets[wc]
        return cfg.datasets[wc.dataset_name]
    return _wrapped


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
        from snakemake.io import ancient

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
    runtime = config.get('runtime', {})
    conda_envs = runtime.get('conda_envs', config.get('conda_envs', {}))

    def _wrapped(env_name):
        user_supplied_env = conda_envs.get(env_name)
        if user_supplied_env is None:
            # use default prespecified environment
            return f'../env_yamls/{env_name}.yaml'
        if os.path.splitext(user_supplied_env)[1] in ('.yaml', '.yml'):
            return os.path.abspath(user_supplied_env)
        return user_supplied_env
    return _wrapped
