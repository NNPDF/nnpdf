import argparse
import pathlib
import logging

from validphys.loader import Loader
from validphys.utils import yaml_safe
from validphys.photon.compute import Photon

def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value." % value)
    return ivalue

log = logging.getLogger(__name__)
loader = Loader()



def main():
    parser = argparse.ArgumentParser(description="FitLux-exec - compute the photon PDF using FiatLux")
    parser.add_argument("config_yml", help="Path to the configuration file", type=pathlib.Path)
    parser.add_argument("replica", help="MC replica number", type=check_positive)
    parser.add_argument(
        "-r",
        "--replica_range",
        help="End of the range of replicas to compute",
        type=check_positive,
    )
    args = parser.parse_args()
    config_path = args.config_yml.absolute()
    replica = args.replica
    if args.replica_range:
        replicas = list(range(replica, args.replica_range + 1))
    else:
        replicas = [replica]

    runcard = yaml_safe.load(config_path.read_text(encoding="UTF-8"))
    theoryID = loader.check_theoryID(runcard["theory"]["theoryid"])
    fiatlux_params = runcard.get("fiatlux", None)
    fiatlux_params['luxset'] = loader.check_pdf(fiatlux_params['luxset'])

    Photon(theoryID, fiatlux_params, replicas=replicas, save_to_disk=True)

if __name__ == "__main__":
    main()