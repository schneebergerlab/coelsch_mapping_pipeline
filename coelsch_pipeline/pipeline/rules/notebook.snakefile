rule render_notebook_template:
    """
    Render Jinja2 templating in a notebook template (markdown cells only).

    Loads a notebook from ../notebook_templates and renders Jinja2 only within
    markdown cells. The template context exposes workflow values via
    top-level names (config/input/output/wildcards/params/threads/resources/log)

    Inputs:
      - notebook template (../notebook_templates)

    Outputs:
      - rendered notebook (results/analysis)
    """
    input:
        template='../notebook_templates/{notebook_template}.{nb_lang}.ipynb'
    output:
        rendered=temp(results('analysis/{dataset_name}.{notebook_template}.template.{nb_lang,py|r}.ipynb'))
    run:
        import json
        from jinja2 import Environment, StrictUndefined

        def _to_native(x):
            # Convert Snakemake IO/params objects into plain python types for Jinja.
            if x is None:
                return None
            if isinstance(x, (str, int, float, bool)):
                return x
            if isinstance(x, dict):
                return {k: _to_native(v) for k, v in x.items()}
            if isinstance(x, (list, tuple, set)):
                return [_to_native(v) for v in x]
            # Namedlist / IOFile / other path-like objects
            try:
                return str(x)
            except Exception:
                return x

        # Build a flat context the templates can reference directly.
        ctx = {
            "config": config,
            "wildcards": {k: getattr(wildcards, k) for k in getattr(wildcards, "keys", lambda: [])()},
            "input": _to_native({k: v for k, v in getattr(input, "items", lambda: [])()}),
            "output": _to_native({k: v for k, v in getattr(output, "items", lambda: [])()}),
            "params": _to_native({k: v for k, v in getattr(params, "items", lambda: [])()}),
            "threads": threads,
            "resources": _to_native(resources),
            "log": _to_native({k: v for k, v in getattr(log, "items", lambda: [])()}),
        }

        env = Environment(undefined=StrictUndefined)

        with open(input.template, "r", encoding="utf-8") as f:
            nb = json.load(f)

        for cell in nb.get("cells", []):
            if cell.get("cell_type") != "markdown":
                continue
            src = cell.get("source", [])
            if isinstance(src, list):
                txt = "".join(src)
                rendered = env.from_string(txt).render(ctx)
                cell["source"] = rendered.splitlines(keepends=True)
            elif isinstance(src, str):
                cell["source"] = env.from_string(src).render(ctx)

        with open(output.rendered, "w", encoding="utf-8") as f:
            json.dump(nb, f, ensure_ascii=False, indent=1)
