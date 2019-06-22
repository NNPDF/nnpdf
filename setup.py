from setuptools import setup, find_packages

setup(
        name="n3fit",
        version="0.9",
        entry_points = {'console_scripts':
            ['n3fit = n3fit.n3fit:main']
            },
        package_dir = {'':'src/'},
        packages=find_packages('src'),
)
