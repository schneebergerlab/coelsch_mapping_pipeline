import re
from collections import defaultdict
from dataclasses import dataclass, field
from glob import glob


_GENOTYPE_TOKEN_RE = re.compile(r'[^*()[\],:\s]+')
_SEX_LABEL_RE = re.compile(r'\[[^\]]*\]')


def parse_crossing_combinations(value):
    """Parse one or more crossing-combination expressions.

    Parameters
    ----------
    value : str or iterable of str
        Crossing-combination expression or collection of expressions. Strings
        are split on commas and stripped of surrounding whitespace.

    Returns
    -------
    list of str
        Parsed crossing-combination expressions.
    """
    if isinstance(value, str):
        return [token.strip() for token in value.split(',') if token.strip()]
    return list(value)


def founder_haplotypes_from_crossing_combinations(crossing_combinations):
    """Extract founder haplotype names from crossing expressions.

    Parameters
    ----------
    crossing_combinations : str or iterable of str
        Crossing-combination expression(s), for example ``"(col0*ler)"`` or
        ``"((c24*tsu0)*(cvi0*no0))"``. Optional genotype-name prefixes before
        a colon and sex labels such as ``[f]`` or ``[m]`` are ignored.

    Returns
    -------
    list of str
        Founder haplotype names in first-seen order, with duplicates removed.

    Raises
    ------
    ValueError
        If an expression has unbalanced parentheses.
    """
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


@dataclass
class SampleDatasetMap:
    """Bidirectional mapping between sample names and dataset names.

    Parameters
    ----------
    sample_to_dataset : dict, optional
        Mapping from sample name to dataset name.
    dataset_to_samples : dict, optional
        Mapping from dataset name to the list of sample names assigned to that
        dataset.
    """

    sample_to_dataset: dict = field(default_factory=dict)
    dataset_to_samples: dict = field(default_factory=dict)

    def dataset_for_sample(self, sample_name):
        """Return the dataset containing a sample.

        Parameters
        ----------
        sample_name : str
            Sample name to resolve.

        Returns
        -------
        str
            Dataset name containing ``sample_name``.
        """
        return self.sample_to_dataset[sample_name]

    def samples_for_dataset(self, dataset_name):
        """Return sample names assigned to a dataset.

        Parameters
        ----------
        dataset_name : str
            Dataset name to resolve.

        Returns
        -------
        list of str
            Sample names assigned to ``dataset_name``.
        """
        return self.dataset_to_samples[dataset_name]


def sample_dataset_map_getter(config, cfg, raw_data, glob_wildcards):
    """Build a sample/dataset map from configured FASTQ basename globs.

    Parameters
    ----------
    config : dict
        Raw pipeline configuration dictionary containing the ``datasets``
        section.
    cfg : object
        Attribute-access pipeline configuration. It must contain
        ``file_suffixes.read1``.
    raw_data : callable
        Path resolver for files under ``raw_data_dir``.
    glob_wildcards : callable
        Snakemake ``glob_wildcards`` function used to extract sample names from
        matched read-one FASTQ paths.

    Returns
    -------
    SampleDatasetMap
        Bidirectional sample-to-dataset and dataset-to-samples mapping.

    Raises
    ------
    ValueError
        If a dataset basename pattern matches no read-one FASTQs, or if one
        sample name matches multiple datasets.
    """
    sample_to_dataset = {}
    dataset_to_samples = defaultdict(list)

    read1_suffix = cfg.file_suffixes.read1
    wildcard_pattern = raw_data(f'{{sample_name}}{read1_suffix}')

    for dataset_name, ds in config['datasets'].items():
        seen_for_dataset = set()

        for basename_glob in ds['input_file_basenames']:
            read1_glob = raw_data(f'{basename_glob}{read1_suffix}')
            sample_names = glob_wildcards(
                wildcard_pattern,
                files=glob(read1_glob),
            ).sample_name

            if not sample_names:
                raise ValueError(
                    f'No read1 FASTQs matched dataset {dataset_name!r} '
                    f'pattern {basename_glob!r}'
                )

            for sample_name in sample_names:
                previous = sample_to_dataset.get(sample_name)
                if previous is not None and previous != dataset_name:
                    raise ValueError(
                        f'Sample {sample_name!r} matched multiple datasets: '
                        f'{previous!r} and {dataset_name!r}'
                    )

                sample_to_dataset[sample_name] = dataset_name
                if sample_name not in seen_for_dataset:
                    dataset_to_samples[dataset_name].append(sample_name)
                    seen_for_dataset.add(sample_name)

    return SampleDatasetMap(
        sample_to_dataset=sample_to_dataset,
        dataset_to_samples=dict(dataset_to_samples),
    )
