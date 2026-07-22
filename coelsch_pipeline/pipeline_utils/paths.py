import os


def resolve_path(basedir, path):
    """Resolve one or more paths relative to a base directory.

    Parameters
    ----------
    basedir : str
        Base directory.
    path : str or list of str
        Relative path or paths to resolve.

    Returns
    -------
    str or list of str
        Path or paths joined to ``basedir``. The return type matches the input
        ``path`` type.
    """
    if isinstance(path, list):
        return [os.path.join(basedir, p) for p in path]
    return os.path.join(basedir, path)


def make_path_getter(config_key):
    """Create a config-backed path resolver factory.

    Parameters
    ----------
    config_key : str
        Slash-delimited path to the config value containing the base directory,
        for example ``"annotation_dir"`` or ``"some/nested/key"``.

    Returns
    -------
    callable
        Factory accepting a raw config dictionary and returning a resolver that
        joins paths to the configured base directory.
    """
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
    """Create a getter for annotation files registered in the config.

    Parameters
    ----------
    config_key : str
        Slash-delimited path to the config mapping containing annotation file
        names, for example ``"annotations/fasta_fns"``.

    Returns
    -------
    callable
        Factory accepting a raw config dictionary and returning a resolver that
        maps annotation keys to files under ``annotation_dir``. Returned paths
        are wrapped with Snakemake ``ancient``.
    """
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
