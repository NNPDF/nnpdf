"""
rebuild_data.py

A script for rebuilding the filtered closure data into the format validphys/n3fit
expect. Run this script on a closure fit either before running n3fit to avoid
crashes due to multiple replicas attempting to rebuild filtered data at same time

Also this script can be used after running a closure fit with nnfit, and will
eradicate the need for the data to be rebuilt by validphys when used in analysis.

Example
-------

If running a closure test with n3fit, simply run script on a filtered closure
fit output directory:

```
$ vp-setupfit fit_name.yaml
$ rebuild-data fit_name
```

and the data is good to go.

If running a closure test with nnfit DO NOT run scipt until after all replicas
have finished.

"""

import logging
import argparse

from reportengine import api, colors

from validphys.app import providers
from validphys.config import Environment
from n3fit.n3fit import N3FitConfig

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log.addHandler(colors.ColorHandler())

loader_log = logging.getLogger('validphys.loader')
loader_log.setLevel(logging.INFO)
loader_log.addHandler(colors.ColorHandler())

# We want to have the Config from n3fit to accept a fit as a directory

REBUILD_CONFIG = dict(
    experiments={"from_": "fit"},
    theory={"from_": "fit"},
    theoryid={"from_": "theory"},
    use_cuts="fromfit",
    # TODO: add namespace specifications to API
    closuretest={"from_": "fit"},
    fakedata={"from_": "closuretest"}
)

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'fit',
        type=str,
        help=(
            "output directory of a closure fit"
        ),
    )
    args = parser.parse_args()
    API = api.API(providers, N3FitConfig, Environment, output=args.fit)
    exps = API.experiments(**REBUILD_CONFIG)

if __name__ == "__main__":
    main()
