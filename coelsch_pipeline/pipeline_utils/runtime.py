import os
from dataclasses import dataclass

from .config import _get_attr_or_key, _get_path_or_key


def setup_tmpdir_getter(cfg):
    """Create a shell snippet factory for rule-local temporary directories.

    Parameters
    ----------
    cfg : object
        Attribute-access pipeline configuration. If ``cfg.runtime.tmpdir`` is
        set, temporary directories are created under that path. Otherwise the
        generated shell code uses ``${TMPDIR:-/tmp}``.

    Returns
    -------
    callable
        Function accepting an optional filename prefix and returning shell code
        that creates ``$tmpdir`` and registers a cleanup trap.
    """
    runtime = _get_attr_or_key(cfg, 'runtime')
    tmpdir = _get_attr_or_key(runtime, 'tmpdir')

    def _wrapped(prefix="{rule}"):
        base = '${TMPDIR:-/tmp}' if tmpdir is None else tmpdir
        return f"""
        tmpdir=$(mktemp -d "{base}/{prefix}.{{jobid}}.XXXXXX");
        trap 'rm -rf -- "$tmpdir"' EXIT;
    """
    return _wrapped


@dataclass
class RuleResources:
    """Resolved Snakemake resource settings for a rule.

    Parameters
    ----------
    threads : int
        Number of threads requested by the rule.
    mem_mb : int, optional
        Total memory requested by the rule in megabytes.
    mem_per_thread : int, optional
        Memory requested per thread in megabytes. When provided, this overrides
        ``mem_mb`` and ``mem_mb`` is resolved to ``threads * mem_per_thread``.
    """

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
    """Create a getter for configured rule resources.

    Parameters
    ----------
    cfg : object
        Attribute-access pipeline configuration containing
        ``runtime.resources``. Resource names can be nested and accessed using
        slash-delimited paths such as ``"samtools/sort"``.

    Returns
    -------
    callable
        Function accepting a resource path and returning a
        :class:`RuleResources` instance. Missing resources fall back to
        ``runtime.resources.default`` or the built-in default.
    """
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


def conda_env_getter(config):
    """Create a getter for configured conda environments.

    Parameters
    ----------
    config : dict
        Raw pipeline configuration. The preferred location is
        ``runtime.conda_envs``; the legacy top-level ``conda_envs`` key is also
        supported.

    Returns
    -------
    callable
        Function accepting an environment name. It returns the configured
        environment path/name, or the bundled ``../env_yamls/<name>.yaml`` path
        when the config value is ``None``. YAML paths are converted to absolute
        paths.
    """
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
