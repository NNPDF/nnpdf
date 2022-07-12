from setuptools import setup, find_packages

setup(
    name="evolven3fit",
    version="1.0",
    packages=find_packages(include=["evolven3fit", "evolven3fit.*"]),
    include_package_data=True,
    entry_points={"console_scripts": ["evolven3fit=evolven3fit.evolven3fit_cli:cli"]},
)
