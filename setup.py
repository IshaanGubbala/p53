from setuptools import setup, find_packages

setup(
    name="p53cad",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "torch>=2.0.0",
        "transformers>=4.30.0",
        "click",
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn"
    ],
    entry_points={
        'console_scripts': [
            'p53cad=p53cad.cli.main:cli',
        ],
    },
)
