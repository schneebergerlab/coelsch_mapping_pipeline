import re


def bool_flag(value, true_flag, false_flag):
    """Format a boolean value as one of two command-line flags.

    Parameters
    ----------
    value : bool or str
        Boolean value to format. String values ``"true"`` and ``"false"`` are
        accepted case-insensitively.
    true_flag : str
        Flag returned when ``value`` is true.
    false_flag : str
        Flag returned when ``value`` is false.

    Returns
    -------
    str
        ``true_flag`` or ``false_flag``.

    Raises
    ------
    ValueError
        If ``value`` is not a recognised boolean value.
    """
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
    """Format a true/false/auto value as a command-line flag.

    Parameters
    ----------
    value : bool, str or None
        Value to format. ``None`` and ``"auto"`` return an empty string.
        Boolean values are passed to :func:`bool_flag`.
    true_flag : str
        Flag returned when ``value`` is true.
    false_flag : str
        Flag returned when ``value`` is false.

    Returns
    -------
    str
        ``true_flag``, ``false_flag`` or an empty string for auto/default
        behavior.
    """
    if value is None:
        return ''
    if isinstance(value, str) and value.lower() == 'auto':
        return ''
    return bool_flag(value, true_flag, false_flag)


def optional_flag(flag, value):
    """Format an optional command-line flag/value pair.

    Parameters
    ----------
    flag : str
        Command-line flag name.
    value : object or None
        Value to attach to ``flag``. If ``None``, no flag is emitted.

    Returns
    -------
    str
        Formatted ``"<flag> <value>"`` string, or an empty string when
        ``value`` is ``None``.
    """
    if value is None:
        return ''
    return f'{flag} {value}'


def format_params(params, indent):
    """Format a list of command parameter fragments for shell output.

    Parameters
    ----------
    params : iterable of object
        Parameter fragments. Empty fragments are skipped. Multi-line fragments
        are split, stripped and re-indented.
    indent : str
        Indentation inserted before each continuation line.

    Returns
    -------
    str
        Formatted shell command fragment with line continuations.
    """
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
    """Format a multi-line shell command for readable Snakemake output.

    Parameters
    ----------
    cmd : str
        Shell command template.

    Returns
    -------
    str
        Command with blank lines removed and line continuations added. Lines
        ending in ``;`` are followed by a blank line to improve readability in
        ``--printshellcmds`` output.
    """
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
