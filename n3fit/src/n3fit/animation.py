#!/usr/bin/env python3

import matplotlib
matplotlib.use("TkAgg")
import glob, os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as mpanimation

from reportengine.compat import yaml

import logging

log = logging.getLogger(__name__)

ALL_PARTICLES = [
        (2, "c", ),
        (3, "s"),
        (4, "u"),
        (5, "d"),
        (6, "g"),
        (7, "anti_d"),
        (8, "anti_u"),
        (9, "anti_s"),
        (10, "anti_c"),
        ]

POS_PASS = "POS_PASS"


def parse_arguments():
    """ Wrapper around ArgumentParser """
    from argparse import ArgumentParser

    parser = ArgumentParser()

    # Main arguments
    parser.add_argument("name_of_fit", help="Folder where to find the history_step_i members", type=str)

    # Optional arguments
    parser.add_argument("-s", "--steps", type=int, default=200)
    parser.add_argument("-b", "--begin_point", type=int, default=1)
    parser.add_argument("-r", "--replicas", type=int, default=20)
    parser.add_argument("-e", "--epochs_per_step", type=int, default=200)
    parser.add_argument("-f", "--frame_separation", type=int, default=200, help="Separation in ms between frames")
    parser.add_argument("-m", "--how_many", type=int, default=25)

    args = parser.parse_args()

    return args


class AnimationObject:
    """
    Class containing all information to generate the animations

    This class provides a `.function` method which can be passed to
    `matplotlib.animation.FuncAnimation`
    """
    def __init__(self, xgrid, step_data, how_many, frame_separation, begin_point=1, alpha=0.4, color="indianred"):
        self.begin_point = begin_point - 1
        self.xgrid = xgrid
        self.alpha = alpha
        self.color = color
        self.step_data = step_data
        self.how_many = how_many
        self.frame_separation = frame_separation

    def define_plot(self, line, ax, particle_tuple, log_scale=False):
        """ 
        """
        self.ax = ax
        self.line = line
        self.particle_index = particle_tuple[0]
        self.particle = particle_tuple[1]
        self.log = log_scale
        self.set_name()

    def set_name(self):
        plt.xlabel("$x$")
        if "anti" in self.particle:
            parname = self.particle.replace("anti_","")
            plt.ylabel(r"$x\bar{" + parname + "}(x)$")
        else:
            plt.ylabel(r"$x" + self.particle + "(x)$")

    def set_title(self, epoch):
        if "anti" in self.particle:
            parname = self.particle.replace("anti_","")
            self.ax.set_title(r"$\bar{{{0}}}$ epoch: {1}".format(parname, epoch))
        else:
            self.ax.set_title(r"${0}$ epoch: {1}".format(self.particle, epoch))

    def initialization(self):
        self.line.set_data([], [])
        return (self.line,)

    def function(self, j):
        i = j + self.begin_point
        step_dict = self.step_data[i]
        epochs = step_dict["epochs"]
        central_val = step_dict["central"][:, self.particle_index]
        miny = step_dict["miny"][:, self.particle_index]
        maxy = step_dict["maxy"][:, self.particle_index]
        tr_chi2 = step_dict["tr_chi2"]
        vl_chi2 = step_dict["vl_chi2"]

        self.ax.collections.clear()
        if self.log:
            self.ax.set_xscale("log")

        self.ax.fill_between(self.xgrid, miny, maxy, alpha=self.alpha, color=self.color)
        self.set_title(epochs)
        self.line.set_data(self.xgrid, central_val)
        self.line.set_color(self.color)

        lab = mpatches.Patch(
            color=self.color, label="$\chi^2_{{vl}}$ = {0:.3f}\n$\chi^2_{{tr}}$ = {1:.3f}".format(vl_chi2, tr_chi2)
        )
        plt.legend(handles=[lab], loc="upper right")
        return (self.line,)


def plot_this(particle_tuple, animator, log_scale=False):
    # First create the figure
    fig = plt.figure()
    ax = plt.axes()
    plot_line, = ax.plot([], [])
    particle = particle_tuple[0]
    particle_name = particle_tuple[1]
    # TODO: Define the axis depending on the particle we are targetting
    if log_scale:
        min_x = 1e-5
        max_x = 1.0
        min_y = -0.05
        max_y = 1.0
    else:
        min_x = 0.0
        max_x = 1.0
        min_y = -0.05
        max_y = 1.0
    ax.axis( [min_x, max_x, min_y, max_y] )

    animator.define_plot(plot_line, ax, particle_tuple, log_scale=log_scale)

    anim = mpanimation.FuncAnimation(
        fig,
        animator.function,
        frames=animator.how_many,
        init_func=animator.initialization,
        interval=animator.frame_separation,
    )

    if log_scale:
        log_string = "log"
    else:
        log_string = ""

    title = "{0}_{1}_{2}.gif".format("animation", particle_name, log_string)
    anim.save(title, writer="pillow")
    print("\n > {0} saved\n".format(title))


def fake_postfit(chi2s, keep=50):
    ave_chi2 = np.average(chi2s)
    order_chi2 = sorted(chi2s, key=lambda x: np.abs(x - ave_chi2))
    # Keep the `keep` best
    throw = order_chi2[keep:]
    # For each one we throw, find what is its index number
    # and it to the kill list
    remove_list = [chi2s.index(t) for t in throw]
    return remove_list


def main(arguments):
    directory = arguments.name_of_fit
    replicas = arguments.replicas
    steps = arguments.steps

    basename = directory + "/history_step_{0}/replica_{1}/" + directory
    # define the names of the interesting files
    info_file_template = basename + ".fitinfo"
    grid_file_template = basename + ".exportgrid"

    # If steps or replicas are not set, try to find out how many there are by yourself
    if steps == -1:
        steps = len(glob.glob(basename.format("*", 1) + "*"))
    if replicas:
        replicas = len(glob.glob(basename.format(0, "*") + "*"))

    remove_indices = []
    indices_removed = False
    xgrid = None
    all_steps_info = []
    # Step 0. Loop over all steps from end to start
    #   The reason we start from the end is because we will use the last step to
    #   see which replicas are useless and remove them, so that all subsequents
    #   iterations of the loop are faster
    for step in range(steps, 0, -1):
        log.info(f"Entering step {step}")
        all_pdfs = []  # List of all PDFs for this step
        all_tr_chi2 = []
        all_vl_chi2 = []
        all_ex_chi2 = []

        for replica in range(1, replicas + 1):  # Loop over replicas
            if replica in remove_indices:
                continue

            # Read the information of the replica
            info_filename = info_file_template.format(step, replica)
            # If the replica is not there, skip
            # this can happen for some steps and not others
            if not os.path.isfile(info_filename):
                continue
            with open(info_filename, "r") as info_file:
                # Step 1. Check whether positivity passes
                first_line = next(info_file).split()
                pos_check = first_line[4]
                if pos_check != POS_PASS:
                    # If positivity does not pass, all steps before this one
                    # wont pass either so add to the ignore list
                    remove_indices.append(replica)
                    continue
                all_vl_chi2.append(float(first_line[1]))
                all_tr_chi2.append(float(first_line[2]))
                all_ex_chi2.append(float(first_line[3]))

            # Read the grid file
            grid_filename = grid_file_template.format(step, replica)
            if not os.path.isfile(grid_filename):
                continue
            with open(grid_filename, "r") as grid_file:
                data = yaml.safe_load(grid_file)
            # We only need to save the grid once because it is always the same
            if not xgrid:
                xgrid = data["xgrid"]
            # Save the PDF data
            pdfgrid = np.array(data["pdfgrid"])
            all_pdfs.append(pdfgrid)

        # Did we get ANY pdf for this step?
        if not all_pdfs:
            log.info(f"No data found for step {step}")
            continue

        # If this is the first time we get to this point, remove indices
        # for doing this we are going to perform a fake postfit
        if not indices_removed:
            indices_removed = True
            to_remove = fake_postfit(all_vl_chi2, keep=arguments.how_many)
            # Now remove the information for the ignored replicas
            # as it wil mess with the average for this one
            for ind in to_remove:
                all_ex_chi2.pop(ind)
                all_vl_chi2.pop(ind)
                all_tr_chi2.pop(ind)
                all_pdfs.pop(ind)
                remove_indices.append(ind)
            log.info(" > Skipping replicas: {0}".format(remove_indices))

        # Now, with the information we have, find the central replica
        np_all_pdfs = np.array(all_pdfs)
        central_replica = np.average(np_all_pdfs, axis=0)

        # Get the avg chi2
        ave_ex_c2 = np.average(all_ex_chi2)
        ave_tr_c2 = np.average(all_tr_chi2)
        ave_vl_c2 = np.average(all_vl_chi2)

        # Get the envelope
        min_values = np.min(np_all_pdfs, axis=0)
        max_values = np.max(np_all_pdfs, axis=0)

        data_dict = {
            "epochs": step * arguments.epochs_per_step,
            "central": central_replica,
            "miny": min_values,
            "maxy": max_values,
            "tr_chi2": ave_tr_c2,
            "vl_chi2": ave_vl_c2,
            "ex_chi2": ave_ex_c2,
            "good_replicas": len(all_tr_chi2),
        }
        all_steps_info.append(data_dict)

    how_many = len(all_steps_info)
    all_steps_info.reverse()  # so that the first one corresponds to step 1

    # Now generate an Animation Object which can be passed down to the plotting function
    animator = AnimationObject(xgrid, all_steps_info, how_many, arguments.frame_separation)

    for particle_tuple in ALL_PARTICLES:
        plot_this(particle_tuple, animator, log_scale = False)
        plot_this(particle_tuple, animator, log_scale = True)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
