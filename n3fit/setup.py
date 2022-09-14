from setuptools import setup, find_packages

setup(
        name="n3fit",
        version="4.0",
        package_dir = {'':'src'},
        packages=find_packages('src'),
        zip_safe=False,
        package_data = {
            '':['*.json', '*.yml', '*.h5'],
            'tests/regressions': ['*'],
        },

        entry_points = {'console_scripts':
            ['n3fit = n3fit.scripts.n3fit_exec:main',
             'evolven3fit = evolven3fit.cli:main',
             'vp-setupfit = n3fit.scripts.vp_setupfit:main',
             'varflavors = n3fit.scripts.varflavors:main',
             ]
            },
)
