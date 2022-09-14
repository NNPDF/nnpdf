"""
This module contains the CLI for evolven3fit
"""
import logging
import pathlib
import numpy as np
from argparse import ArgumentParser

from . import evolve
from validphys.core import CommonDataMetadata

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def main():
    parser = ArgumentParser(
        description="evolven3fit - a script with tools to evolve PDF fits"
    )
    subparsers = parser.add_subparsers(title="actions", dest="actions")
    evolvefit_parser = subparsers.add_parser(
        "evolve",
        help="Evolves the fitted PDFs. The q_grid starts at the Q0 given by the theory but the last point is q_fin and its number of points can be specified by q_points. If a path is given for the dump option, the eko will be dumped in that path after the computation. If a path is given for the load option, the eko to be used for the evolution will be loaded from that path. The two options are mutually exclusive.",
    )
    eko_parser = subparsers.add_parser(
        "produce_eko",
        help="Produce the eko for the specified theory. The q_grid starts at the Q0 given by the theory but the last point is q_fin and its number of points can be specified by q_points. The x_grid starts at x_grid_ini and ends at 1.0 and contains the provided number of points. The eko will be dumped in the provided path.",
    )
    evolvefit_parser.add_argument(
        "configuration_folder",
        help="Path to the folder containing the (pre-DGLAP) fit result",
    )
    parser.add_argument(
        "-q",
        "--q_fin",
        type=float,
        default=None,
        help="Final q-value of the evolution",
    )
    parser.add_argument(
        "-p",
        "--q_points",
        type=int,
        default=None,
        help="Number of q points for the evolution",
    )
    parser.add_argument(
        "-n",
        "--n_cores",
        type=int,
        default=1,
        help="Number of cores to be used",
    )
    evolvefit_parser.add_argument(
        "-l",
        "--load",
        type=pathlib.Path,
        default=None,
        help="Path of the EKO to be loaded",
    )
    evolvefit_parser.add_argument(
        "-d",
        "--dump",
        type=pathlib.Path,
        default=None,
        help="Path where the EKO is dumped (optional)",
    )
    eko_parser.add_argument(
        "theoryID", type=int, help="ID of the theory used to produce the eko"
    )
    eko_parser.add_argument(
        "dump",
        type=pathlib.Path,
        help="Path where the EKO is dumped",
    )
    eko_parser.add_argument(
        "-i" "--x_grid_ini",
        default=None,
        type=float,
        help="Starting point of the x-grid",
    )
    eko_parser.add_argument(
        "-p" "--x_grid_points",
        default=None,
        type=int,
        help="Number of points of the x-grid",
    )
    args = parser.parse_args()
    op_card_info = {
        "ev_op_max_order": 10,
        "ev_op_iterations": 1,
        "n_integration_cores": args.n_cores,
        "backward_inversion": "expanded",
    }
    t_card_info = {}
    if args.actions == "evolve":
        cli_evolven3fit(
            args.configuration_folder,
            args.q_fin,
            args.q_points,
            op_card_info,
            t_card_info,
            args.dump,
            args.load,
        )
    elif args.actions == "produce_eko":
        if args.x_grid_ini is None:
            if args.x_grid_points is None:
                x_grid = np.array(
                    [
                        1.00000000000000e-09,
                        1.29708482343957e-09,
                        1.68242903474257e-09,
                        2.18225315420583e-09,
                        2.83056741739819e-09,
                        3.67148597892941e-09,
                        4.76222862935315e-09,
                        6.17701427376180e-09,
                        8.01211109898438e-09,
                        1.03923870607245e-08,
                        1.34798064073805e-08,
                        1.74844503691778e-08,
                        2.26788118881103e-08,
                        2.94163370300835e-08,
                        3.81554746595878e-08,
                        4.94908707232129e-08,
                        6.41938295708371e-08,
                        8.32647951986859e-08,
                        1.08001422993829e-07,
                        1.40086873081130e-07,
                        1.81704331793772e-07,
                        2.35685551545377e-07,
                        3.05703512595323e-07,
                        3.96522309841747e-07,
                        5.14321257236570e-07,
                        6.67115245136676e-07,
                        8.65299922973143e-07,
                        1.12235875241487e-06,
                        1.45577995547683e-06,
                        1.88824560514613e-06,
                        2.44917352454946e-06,
                        3.17671650028717e-06,
                        4.12035415232797e-06,
                        5.34425265752090e-06,
                        6.93161897806315e-06,
                        8.99034258238145e-06,
                        1.16603030112258e-05,
                        1.51228312288769e-05,
                        1.96129529349212e-05,
                        2.54352207134502e-05,
                        3.29841683435992e-05,
                        4.27707053972016e-05,
                        5.54561248105849e-05,
                        7.18958313632514e-05,
                        9.31954227979614e-05,
                        1.20782367731330e-04,
                        1.56497209466554e-04,
                        2.02708936328495e-04,
                        2.62459799331951e-04,
                        3.39645244168985e-04,
                        4.39234443000422e-04,
                        5.67535660104533e-04,
                        7.32507615725537e-04,
                        9.44112105452451e-04,
                        1.21469317686978e-03,
                        1.55935306118224e-03,
                        1.99627451141338e-03,
                        2.54691493736552e-03,
                        3.23597510213126e-03,
                        4.09103436509565e-03,
                        5.14175977083962e-03,
                        6.41865096062317e-03,
                        7.95137940306351e-03,
                        9.76689999624100e-03,
                        1.18876139251364e-02,
                        1.43298947643919e-02,
                        1.71032279460271e-02,
                        2.02100733925079e-02,
                        2.36463971369542e-02,
                        2.74026915728357e-02,
                        3.14652506132444e-02,
                        3.58174829282429e-02,
                        4.04411060163317e-02,
                        4.53171343973807e-02,
                        5.04266347950069e-02,
                        5.57512610084339e-02,
                        6.12736019390519e-02,
                        6.69773829498255e-02,
                        7.28475589986517e-02,
                        7.88703322292727e-02,
                        8.50331197801452e-02,
                        9.13244910278679e-02,
                        9.77340879783772e-02,
                        1.04252538208639e-01,
                        1.10871366547237e-01,
                        1.17582909372878e-01,
                        1.24380233801599e-01,
                        1.31257062945031e-01,
                        1.38207707707289e-01,
                        1.45227005135651e-01,
                        1.52310263065985e-01,
                        1.59453210652156e-01,
                        1.66651954293987e-01,
                        1.73902938455578e-01,
                        1.81202910873333e-01,
                        1.88548891679097e-01,
                        1.95938145999193e-01,
                        2.03368159629765e-01,
                        2.10836617429103e-01,
                        2.18341384106561e-01,
                        2.25880487124065e-01,
                        2.33452101459503e-01,
                        2.41054536011681e-01,
                        2.48686221452762e-01,
                        2.56345699358723e-01,
                        2.64031612468684e-01,
                        2.71742695942783e-01,
                        2.79477769504149e-01,
                        2.87235730364833e-01,
                        2.95015546847664e-01,
                        3.02816252626866e-01,
                        3.10636941519503e-01,
                        3.18476762768082e-01,
                        3.26334916761672e-01,
                        3.34210651149156e-01,
                        3.42103257303627e-01,
                        3.50012067101685e-01,
                        3.57936449985571e-01,
                        3.65875810279643e-01,
                        3.73829584735962e-01,
                        3.81797240286494e-01,
                        3.89778271981947e-01,
                        3.97772201099286e-01,
                        4.05778573402340e-01,
                        4.13796957540671e-01,
                        4.21826943574548e-01,
                        4.29868141614175e-01,
                        4.37920180563205e-01,
                        4.45982706956990e-01,
                        4.54055383887562e-01,
                        4.62137890007651e-01,
                        4.70229918607142e-01,
                        4.78331176755675e-01,
                        4.86441384506059e-01,
                        4.94560274153348e-01,
                        5.02687589545177e-01,
                        5.10823085439086e-01,
                        5.18966526903235e-01,
                        5.27117688756998e-01,
                        5.35276355048428e-01,
                        5.43442318565661e-01,
                        5.51615380379768e-01,
                        5.59795349416641e-01,
                        5.67982042055800e-01,
                        5.76175281754088e-01,
                        5.84374898692498e-01,
                        5.92580729444440e-01,
                        6.00792616663950e-01,
                        6.09010408792398e-01,
                        6.17233959782450e-01,
                        6.25463128838069e-01,
                        6.33697780169485e-01,
                        6.41937782762089e-01,
                        6.50183010158361e-01,
                        6.58433340251944e-01,
                        6.66688655093089e-01,
                        6.74948840704708e-01,
                        6.83213786908386e-01,
                        6.91483387159697e-01,
                        6.99757538392251e-01,
                        7.08036140869916e-01,
                        7.16319098046733e-01,
                        7.24606316434025e-01,
                        7.32897705474271e-01,
                        7.41193177421404e-01,
                        7.49492647227008e-01,
                        7.57796032432224e-01,
                        7.66103253064927e-01,
                        7.74414231541921e-01,
                        7.82728892575836e-01,
                        7.91047163086478e-01,
                        7.99368972116378e-01,
                        8.07694250750291e-01,
                        8.16022932038457e-01,
                        8.24354950923382e-01,
                        8.32690244169987e-01,
                        8.41028750298844e-01,
                        8.49370409522600e-01,
                        8.57715163684985e-01,
                        8.66062956202683e-01,
                        8.74413732009721e-01,
                        8.82767437504206e-01,
                        8.91124020497459e-01,
                        8.99483430165226e-01,
                        9.07845617001021e-01,
                        9.16210532771399e-01,
                        9.24578130473112e-01,
                        9.32948364292029e-01,
                        9.41321189563734e-01,
                        9.49696562735755e-01,
                        9.58074441331298e-01,
                        9.66454783914439e-01,
                        9.74837550056705e-01,
                        9.83222700304978e-01,
                        9.91610196150662e-01,
                        1.00000000000000e00,
                    ]
                ).reshape(-1, 1)
            else:
                raise ValueError(
                    "x_grid_ini and x_grid_points must be specified either both or none of them"
                )
        elif args.x_grid_points is None:
            raise ValueError(
                "x_grid_ini and x_grid_points must be specified either both or none of them"
            )
        else:
            x_grid = np.geomspace(args.x_grid_ini, 1.0, args.x_grid_points)
        tcard, opcard = evolve.construct_eko_cards(
            args.theoryID, op_card_info, t_card_info, args.q_fin, args.q_points, x_grid
        )
        eko_op = evolve.construct_eko_for_fit(tcard, opcard, args.dump)


def cli_evolven3fit(
    configuration_folder, q_fin, q_points, op_card_info, t_card_info, dump, load
):
    """Evolves the fitted PDFs.

    The q_grid starts at the Q0 given by the theory but
    the last point is q_fin and its number of
    points can be specified by q_points. If just one of the
    two is not specified by the user, the default grid
    will be used.

    If a path is given for the dump option, the eko
    will be dumped in that path after the computation.

    If a path is given for the load option, the eko
    to be used for the evolution will be loaded from that
    path.

    The two options are mutually exclusive.
    """
    return evolve.evolve_fit(
        configuration_folder, q_fin, q_points, op_card_info, t_card_info, load, dump
    )


if __name__ == "__main__":
    main()
