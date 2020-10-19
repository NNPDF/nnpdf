from setuptools import setup, find_packages

setup(
        name="n3fit",
        version="0.9",
        package_dir = {'':'src'},
        packages=find_packages('src'),
        zip_safe=False,
        package_data = {
            '':['*.fitinfo', '*.yml', '*.h5'],
            'tests/regressions': ['*'],
        },

        entry_points = {'console_scripts':
            ['n3fit = n3fit.scripts.n3fit_exec:main',
             ]
            },
)
