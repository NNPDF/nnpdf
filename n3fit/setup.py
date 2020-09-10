from setuptools import setup, find_packages

setup(
        name="n3fit",
        version="0.9",
        package_dir = {'':'src'},
        packages=find_packages('src'),
        zip_safe=False,
        package_data = {
            '':['*.fitinfo', '*.yml'],
            'tests/regressions': ['*'],
        },
        extras_require={
          'ga' : [
            'evolutionary-keras',
            ],
          },

        entry_points = {'console_scripts':
            ['n3fit = n3fit.scripts.n3fit_exec:main',
             'n3Hyperplot = n3fit.hyper_optimization.plotting:main',
             ]
            },
)
