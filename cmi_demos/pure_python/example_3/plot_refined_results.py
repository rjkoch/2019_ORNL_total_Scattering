# Tutorial 3 Refining PDFs from SrFe2As2 to investigate phase transformation
#
# This python script will help yo parse out the refined information.
# Tutorial 3 needs to be run before running this script.
#
# Import packages that we will need.
import os
import yaml
import sys

import matplotlib.pyplot as plt
import numpy as np

# This is the same as in the refinement file.
FIT_ID_BASE = "Fit_SrFe2As2_"

# This is where we will save the plots we make here.
T_SRIES_PLOT_DIR = "T_series_plots"


def main():
    """
    This will run by default when the file is executed using
    "python file.py" in the command line

    Parameters
    ----------
    None

    Returns
    ----------
    None
    """

    # We now read a dictionary from a file, if it exists
    # This .yaml file should have been created after running tutorial 3.
    if os.path.exists(FIT_ID_BASE + "refined_params.yml"):
        with open(FIT_ID_BASE + "refined_params.yml", 'r') as infile:
            refined_dict = yaml.safe_load(infile)

    # If the file doesn't exist, print a message to the terminal and exit
    elif not os.path.exists(FIT_ID_BASE + "refined_params.yml"):
        print(FIT_ID_BASE + "refined_params.yml does not exist.")

        # If we dont find the yaml file, we have nothing to do
        # so we exit.
        sys.exit(1)

    # Make a place to put the saved plots if one doesn't exist.
    if not os.path.exists(T_SRIES_PLOT_DIR):
        os.makedirs(T_SRIES_PLOT_DIR)

    # Use a standard style format.
    plt.style.use(os.path.join(os.pardir,
                               "utils",
                               "billinge.mplstyle"))

    # Close any plots if they're open already.
    plt.clf()
    plt.close('all')

    # Create a subplot axis.
    fig, ax = plt.subplots(1, 1)

    # Let's look at the c lattice parameter.
    var_name = "c"

    # We loop over all the structures we considered
    # In this case it's just two.
    for structure in refined_dict.keys():

        # We get the x values we want to plot against
        # Here, it's the temperature.
        xs = [temp for temp in refined_dict[structure].keys()]

        # We get the y values we want to plot
        # Here, it's defined by the string in "var_name"
        # above, or the c lattice parameter.
        ys = [refined_dict[structure][temp][var_name] for
               temp in refined_dict[structure].keys()]

        # In case they aren't sorted, we sort by temperature (x).
        xs, ys= zip(*sorted(zip(xs, ys)))

        # Plot the data, using the structure type as the label.
        ax.plot(xs,
                ys,
                '-o',
                label=structure)

    # Outside the loop, we label the axes.
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r"Lattice parameter, c ($\rm \AA$)")

    # Outside the loop, we create a legend.
    plt.legend()

    # And we show the plot.
    plt.show()

    fig.savefig(os.path.join(T_SRIES_PLOT_DIR, 'full_' + var_name + '_T.pdf'),
                format='pdf')

    # Close any plots if they're open already
    plt.clf()
    plt.close('all')


    # Create a subplot axis.
    fig, ax = plt.subplots(1, 1)

    # Let's look at the fit quality rw.
    var_name = "rw"

    # We loop over all the structures we considered
    # In this case it's just two.
    for structure in refined_dict.keys():

        # We get the x values we want to plot against
        # Here, it's the temperature.
        xs = [temp for temp in refined_dict[structure].keys()]

        # We get the y values we want to plot
        # Here, it's defined by the string in "var_name"
        # above, or rw.
        ys = [refined_dict[structure][temp][var_name] for
               temp in refined_dict[structure].keys()]

        # In case they aren't sorted, we sort by temperature (x).
        xs, ys= zip(*sorted(zip(xs, ys)))

        # Plot the data, using the structure type as the label.
        ax.plot(xs,
                ys,
                '-o',
                label=structure)

    # Outside the loop, we label the axes.
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r"$\rm R_w$")

    # Outside the loop, we create a legend.
    plt.legend()

    # And we show the plot.
    plt.show()

    fig.savefig(os.path.join(T_SRIES_PLOT_DIR, 'full_' + var_name + '_T.pdf'),
                format='pdf')

    # Close any plots if they're open already
    plt.clf()
    plt.close('all')

    # Create a subplot axis.
    fig, ax = plt.subplots(1, 1)

    # Let's look at all the some thermal parameters.


    # We loop over all the structures we considered
    # In this case it's just two.
    for structure in refined_dict.keys():
        first_key = next(iter(refined_dict["orthorhombic"].keys()))
        thermal_names = [key for key in refined_dict[structure][first_key].keys() if "U" in key]

        for thermal in thermal_names:

            # We get the x values we want to plot against
            # Here, it's the temperature.
            xs = [temp for temp in refined_dict[structure].keys()]

            # We get the y values we want to plot
            # Here, it's defined by the string in "var_name"
            # above, or the thermal parameter.
            ys = [refined_dict[structure][temp][thermal] for
                  temp in refined_dict[structure].keys()]

            # In case they aren't sorted, we sort by temperature (x).
            xs, ys = zip(*sorted(zip(xs, ys)))

            # Plot the data, using the structure type as the label.
            ax.plot(xs,
                    ys,
                    '-o',
                    label=structure + "_" + thermal)




    # Outside the loop, we label the axes.
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r"Thermal paramete ($\rm \AA^2$)")

    # Outside the loop, we create a legend.
    plt.legend(fontsize=10, ncol=2)

    # And we show the plot.
    plt.show()


    fig.savefig(os.path.join(T_SRIES_PLOT_DIR, 'full_thermals_T.pdf'),
                format='pdf')

    # Close any plots if they're open already
    plt.clf()
    plt.close('all')

    # Create a subplot axis
    fig, ax = plt.subplots(1, 1)

    structure = "orthorhombic"

    var_name = "a"
    a_params = [refined_dict[structure][temp][var_name]
                for temp in refined_dict[structure].keys()]

    var_name = "b"
    b_params = [refined_dict[structure][temp][var_name]
                for temp in refined_dict[structure].keys()]


    xs = [temp for temp in refined_dict[structure].keys()]

    ys = [(a-b)/(0.5*(a+b)) for
          a, b in zip(a_params, b_params)]

    xs, ys = zip(*sorted(zip(xs, ys)))

    ax.plot(xs,
            ys,
            '-o')

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r"Orthorhombicity ($\rm \AA$)")

    plt.legend()
    plt.show()

    fig.savefig(os.path.join(T_SRIES_PLOT_DIR, 'orthorhombicity_T.pdf'),
                format='pdf')


    # End of function


if __name__ == "__main__":
    main()

# End of file