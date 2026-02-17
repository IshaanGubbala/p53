from setuptools import setup, find_packages

setup(
    name="p53cad",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "torch>=2.0.0",
        "transformers>=4.30.0,<5.0.0",
        "click>=8.0",
        "pandas>=2.0",
        "numpy>=1.26.0",
        "matplotlib>=3.7",
        "seaborn>=0.11",
        "scipy>=1.10",
        "pyarrow>=12.0",
        "requests>=2.31",
        "streamlit>=1.30",
        "plotly>=5.0",
    ],
    extras_require={
        "drug": ["rdkit", "meeko"],
        "md": ["openmm", "openff-toolkit", "openff-interchange", "openmmforcefields"],
        "explainability": ["scipy>=1.10"],
    },
    entry_points={
        'console_scripts': [
            'p53cad=p53cad.cli.main:cli',
        ],
    },
)
