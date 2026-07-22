from types import SimpleNamespace


def config_namespace(obj, keep_dict_keys=frozenset({'datasets'})):
    """Convert nested config dictionaries to attribute-access namespaces.

    Parameters
    ----------
    obj : object
        Configuration object to convert. Dictionaries are converted
        recursively, lists are converted element-wise, and scalar values are
        returned unchanged.
    keep_dict_keys : frozenset of str, optional
        Dictionary keys that should remain dictionaries instead of being
        converted to :class:`types.SimpleNamespace`. This is useful for config
        sections such as ``datasets`` where arbitrary user-defined keys are
        expected.

    Returns
    -------
    object
        Converted configuration object. Nested dictionaries are accessible by
        attribute unless their key is listed in ``keep_dict_keys``.
    """
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


def dataset_getter(cfg):
    """Create a getter for dataset-level configuration.

    Parameters
    ----------
    cfg : object
        Pipeline configuration, usually produced by :func:`config_namespace`.
        It must contain a ``datasets`` mapping.

    Returns
    -------
    callable
        Function accepting either a dataset name string or a Snakemake
        wildcards object with ``dataset_name`` and returning the corresponding
        dataset configuration object.
    """
    def _wrapped(wc):
        if isinstance(wc, str):
            return cfg.datasets[wc]
        return cfg.datasets[wc.dataset_name]
    return _wrapped
