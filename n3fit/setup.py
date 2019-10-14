from setuptools import setup, find_packages

setup(
        name="n3fit",
        version="0.9",
        package_dir = {'':'src'},
        packages=find_packages('src'),

        entry_points = {'console_scripts':
            ['n3fit = n3fit.n3fit:main',
             'n3Hyperplot = n3fit.hyper_optimization.plotting:main',
             ]
            },
)
