from setuptools import setup, find_packages

setup(
    name="coelsch_pipeline",
    version="0.1.1",
    description="Pipeline tool for mapping reads for crossover analysis with coelsch",
    author="Matthew Parker",
    author_email="mparker@mpipz.mpg.de",
    packages=find_packages(),
    install_requires=[
        "click>=8.0",
        "jinja2",
        "snakemake",
    ],
    entry_points={
        "console_scripts": [
            "coelsch_pipeline = coelsch_pipeline.cli:cli"
        ]
    },
    include_package_data=True,
    zip_safe=False,
)
