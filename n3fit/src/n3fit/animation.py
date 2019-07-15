#!/usr/bin/env python3
import lhapdf
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'
import matplotlib.animation as animation
import sys
import os
import numpy as np

import matplotlib.patches as mpatches

from reportengine.compat import yaml

from ipdb import set_trace



# Monte Carlo Particle Numbering Scheme
particle_names = {1: 'u', 2: 'd', 3: 's', 4: 'c', 5: 'b', 6: 't', 21: 'g',
        -1: 'anti_u', -2: 'anti_d', -3: 'anti_s', -4: 'anti_c', -5: 'anti_b', -6: 'anti_t'}

def get_particle_index(par_number):
    if par_number < -6 or (par_number > 6 and par_number != 21) or par_number == 0:
        print("""Error: wrong number for the particle
        {0}""".format(particles))
        sys.exit()
    gluon = 6
    if par_number == 21:
        return gluon
    else:
        return par_number + gluon

def aparse():
    from argparse import ArgumentParser
    parser = ArgumentParser()

    # Main arguments
    parser.add_argument('name_of_fit', type = str)

    parser.add_argument('-s', '--steps', type = int, default = 200)
    parser.add_argument('-b', '--begin_point', type = int, default = 1)
    parser.add_argument('-r', '--replicas', type = int, default = 20)
    parser.add_argument('-p', '--particle', type = int, default = [21], nargs = "+")
    parser.add_argument('-e', '--epochs_per_step', type = int, default = 200)
    parser.add_argument('-f', '--frame_separation', type = int, default = 200, help = "Separation in ms between frames")
    parser.add_argument('-m', '--how_many', type = int, default = 25)

    args = parser.parse_args()

    return args

class AnimationObject:
    def __init__(self, xgrid, step_data,
            line, ax, particle,
            begin_point = 1, log_scale = False,
            alpha = 0.4, color = 'indianred'):
        self.xgrid = xgrid
        self.step_data = step_data

        self.begin_point = begin_point - 1
        self.ax = ax
        self.line = line
        self.particle = particle
        self.particle_index = get_particle_index(particle)

        self.alpha = alpha
        self.color = color
        self.log = log_scale


    def set_title(self, epoch):
        if self.particle > 0:
            self.ax.set_title(
                    r"${0}$ epoch: {1}".format(particle_names[self.particle], epoch)
                    )
        else:
            self.ax.set_title(
                    r"$\bar{{{0}}}$ epoch: {1}".format(particle_names[-self.particle], epoch)
                    )

    def initialization(self):
        self.line.set_data([],[])
        return self.line,

    def function(self, j):
        i = j + self.begin_point
        step_dict = self.step_data[i]
        epochs = step_dict['epochs']
        central_val = step_dict['central'][:, self.particle_index]
        miny = step_dict['miny'][:, self.particle_index]
        maxy = step_dict['maxy'][:, self.particle_index]
        tr_chi2 = step_dict['tr_chi2']
        vl_chi2 = step_dict['vl_chi2']

        self.ax.collections.clear()
        if self.log:
            self.ax.set_xscale('log')

        self.ax.fill_between(
                self.xgrid, miny, maxy, alpha = self.alpha, color = self.color
                )
        self.set_title(epochs)
        self.line.set_data(self.xgrid, central_val)
        self.line.set_color(self.color)

        lab = mpatches.Patch(
                color = self.color,
                label = "$\chi^2_{{vl}}$ = {0:.3f}\n$\chi^2_{{tr}}$ = {1:.3f}".format(vl_chi2, tr_chi2),
                )
        plt.legend(handles = [lab], loc = 'upper right')
        return self.line,

def main(arguments):

    directory = arguments.name_of_fit

    steps = arguments.steps
    begin_point = arguments.begin_point
    replicas = arguments.replicas
    particles = arguments.particle
    epochs_per_step = arguments.epochs_per_step

    for p in particles:
        # Safety check
        _ = get_particle_index(p)

    # Now instead of using LHAPDF let us read the replicas by ourselves
    basename = directory + "/history_step_{0}/replica_{1}/" + directory
    # define the names of the interesting files\
    info_file = basename + ".fitinfo"
    grid_file = basename + ".exportgrid"

    xgrid = []
    all_steps_info = []
    remove_indices = []
    for step in range(steps, 0, -1): # Loop over steps
        print("Entering on step {0}".format(step))
        all_pdfs = []
        all_tr_chi2 = []
        all_vl_chi2 = []
        all_ex_chi2 = []

        for replica in range(1, replicas+1): # Loop over replicas
            if replica in remove_indices:
                continue
            # Read the information of the replica
            set_trace()
            info = info_file.format(step, replica)
            # If the replica is not there, skip
            if not os.path.isfile(info):
                continue
            with open(info, 'r') as filein:
                line_str = next(filein)
                line = line_str.split()
                pos_pass = line[4]
                if pos_pass != 'POS_PASS':
                    continue
                all_vl_chi2.append( float(line[1]) )
                all_tr_chi2.append( float(line[2]) )
                all_ex_chi2.append( float(line[3]) )
            # Read the exportgrid file
            grid = grid_file.format(step, replica)
            if not os.path.isfile(grid):
                continue
            with open(grid, 'r') as filegrid:
                data = yaml.safe_load(filegrid)
            if not xgrid: # We only need to load the xgrid once
                xgrid = data['xgrid']
            # Save the complete PDF for now for each replica
#             pdfgrid = np.array(data['pdfgrid'])[:,get_particle_index(particle)]
            pdfgrid = np.array(data['pdfgrid'])
            all_pdfs.append( pdfgrid )
        # Check whether any of the replicas contained data for this run
        if not all_pdfs:
            print("No data found for step {0}".format(step))
            continue
        if not remove_indices:
            use_chi2 = all_vl_chi2
            # Now do a fake postfit looking at the chi2
            ave_c2 = np.average(use_chi2)
            order = sorted(use_chi2, key = lambda x: np.abs(x-ave_c2))
            # Keep only the 25 best
            throw = order[arguments.how_many:]
            for t in throw:
                ind = use_chi2.index(t)
                remove_indices.append( ind )
                all_ex_chi2.pop(ind)
                all_vl_chi2.pop(ind)
                all_tr_chi2.pop(ind)
                all_pdfs.pop(ind)
            print(" > Skipping replicas {0}".format(remove_indices))

        # Find the central replica
        np_all_pdf = np.array(all_pdfs)
        central_replica = np.average(np_all_pdf, axis = 0)

        # Get the average chi2
        ave_ex_c2 = np.average(all_ex_chi2)
        ave_tr_c2 = np.average(all_tr_chi2)
        ave_vl_c2 = np.average(all_vl_chi2)

        # Now save the information for this step
        min_values = np.min(np_all_pdf, axis = 0)
        max_values = np.max(np_all_pdf, axis = 0)

        data_dict = {
                'epochs' : step*epochs_per_step,
                'central' : central_replica,
                'miny' : min_values,
                'maxy': max_values,
                'tr_chi2' : ave_tr_c2,
                'vl_chi2' : ave_vl_c2,
                'ex_chi2' : ave_ex_c2,
                'good_replicas' : len(all_tr_chi2),
                }
        all_steps_info.append( data_dict )

    # Now we have all information, we have to plot
    how_many = len(all_steps_info)
    all_steps_info.reverse()

    # Now create a subfunction with all plot routines
    def plot_this(particle, log_scale = False):
            # Create the figure
            fig = plt.figure()
            ax = plt.axes()
            plot_line, = ax.plot([], [])

            # Define the axis depending on the particle we are targetting
            if log_scale:
                if particle == -1 or particle == -2:
                    ax.axis([1e-5,1,-0.05,0.65])
                elif particle == -3:
                    ax.axis([1e-5,1,-0.3,0.4])
                elif particle == 1:
                    ax.axis([1e-5,1,-0.05,0.65])
                elif particle == 2:
                    ax.axis([1e-5,1,-0.01,0.80])
                elif particle == 3:
                    ax.axis([1e-5,1,-0.3,0.4])
                elif particle == 4:
                    ax.axis([1e-5,1,-0.005,0.04])
                else:
                    ax.axis([1e-5,1,-0.05,3.2])
            else:
                if particle == -1 or particle == -2:
                    ax.axis([0,1,-0.05,0.16])
                elif particle == -3:
                    ax.axis([0,1,-0.011,0.2])
                elif particle == 1:
                    ax.axis([0,1,-0.05,0.5])
                elif particle == 2:
                    ax.axis([0,1,-0.01,0.80])
                elif particle == 3:
                    ax.axis([0,1,-0.01,0.07])
                elif particle == 4:
                    ax.axis([0,1,-0.005,0.035])
                else:
                    ax.axis([0,1,-0.05,1.5])

            plt.xlabel("$x$")
            if particle > 0:
                plt.ylabel(r"$x"+particle_names[particle]+"(x)$")
            else:
                plt.ylabel(r"$x\bar{"+particle_names[-particle]+"}(x)$")

            # First generate the AnimationObject that carries all necessary information
            animation_object = AnimationObject(xgrid, all_steps_info,
                    plot_line, ax, particle,
                    log_scale = log_scale)

            anim = animation.FuncAnimation(
                    fig,
                    animation_object.function,
                    frames = how_many,
                    init_func = animation_object.initialization,
                    interval = arguments.frame_separation,
                    )
            if log_scale:
                log_string = "log"
            else:
                log_string = ""
            title = "{0}_{1}_{2}.gif".format(directory, particle_names[particle], log_string)
            anim.save(title, writer = 'pillow')
            print("\n > {0} saved\n".format(title))

    for particle in particles:
        plot_this(particle, log_scale = False)
        plot_this(particle, log_scale = True)

if __name__ == "__main__":

    args = aparse()

    main(args)
