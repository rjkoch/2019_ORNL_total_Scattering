#!/usr/bin/env python
###################################
#                                 #
# File coded by: Robert J. Koch   #
#                                 #
###################################
#
# Example 1, Refining PDF from Ni to characterize instrument
#
# This Diffpy-CMI notebook will carry out a structural refinement of a measured
# PDF from nickel.  Nickel is a standard, and allows us to understand how the
# instrument affects our data, so that we can correct for it when we
# look at an unknown sample (in example 2).
#
# Import packages that we will need.
import sys
from pathlib import Path

import yaml

sys.path.append(str(Path().absolute().parent.parent.parent))

from diffpy.srfit.fitbase import FitResults
from diffpy.structure.parsers import getParser
from scipy.optimize import least_squares
from cmi_demos.utils.helpers import makerecipe
from cmi_demos.utils.helpers import plotresults

############### Config ##############################
# Give a file path to where your pdf (.gr) and (.cif) files are located.
# In this case it is in a folder called "data," in the same directory as
# this script.
DPATH = Path("data")

# Give an identifying name for the refinement, similar
# to what you would name a fit tree in PDFGui.
FIT_ID = "Fit_Ni_Bulk"

# Specify the names of the input PDF and cif files.
GR_NAME = "Ni.gr"
CIF_NAME = "Ni.cif"

######## Experimental PDF Config ######################
# Specify the min, max, and step r-values of the PDF (that we want to fit over)
# also, specify the Q_max and Q_min values used to reduce the PDF.
PDF_RMIN = 1.5
PDF_RMAX = 50
PDF_RSTEP = 0.01
QMAX = 25
QMIN = 0.1

######## PDF Initialize refinable variables #############
# We explicitly specify the initial values for the lattice
# parameter, scale, and isotropic thermal parameters, as well
# as a correlated motion parameter, in this case delta_2.
CUBICLAT_I = 3.52
SCALE_I = 0.4
UISO_I = 0.005
DELTA2_I = 2

# We give initial values for the instrumental parameters, but because
# this is a calibrant, we will also refine these values.
QDAMP_I = 0.04
QBROAD_I = 0.02


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

    # Make some folders to store our output files.
    resdir = Path("res")
    fitdir = Path("fit")
    figdir = Path("fig")

    folders = [resdir, fitdir, figdir]

    # Loop over all folders
    for folder in folders:

        # If the folder does not exist...
        if not folder.exists():
            # ...then we create it.
            folder.mkdir()

    # Let the user know what fit we are running by printing to terminal.
    basename = FIT_ID
    print(f"\n{basename}\n")

    # Establish the full location of the data.
    data = DPATH / GR_NAME

    # Establish the location of the cif file with the structure of interest
    # and load it into a diffpy structure object.
    strudir = DPATH
    cif_file = strudir / CIF_NAME

    p_cif = getParser('cif')
    structure = p_cif.parseFile(str(cif_file))
    space_group = p_cif.spacegroup.short_name

    # Initialize the Fit Recipe by giving it this diffpy structure
    # as well as the path to the data file.
    recipe = makerecipe(cif_file, data)

    # We first want to add a scale parameter to our fit recipe!
    recipe.addVar(recipe.crystal.s1, SCALE_I, tag="scale")

    # Let's set the calculation range!
    recipe.crystal.profile.setCalculationRange(xmin=PDF_RMIN, xmax=PDF_RMAX, dx=PDF_RSTEP)

    # assign Q_max and Q_min, all part of the PDF Generator object.
    # It's possible that the PDFParse function we used above
    # already parsed out ths information, but in case it didn't, we set it
    # explicitly again here.
    recipe.crystal.G1.setQmax(QMAX)
    recipe.crystal.G1.setQmin(QMIN)

    # Initialize and add the instrument parameters, Q_damp and Q_broad, and
    # delta and instrumental parameters to Fit Recipe.
    # We give them unique names, and tag them with our choice of relevant strings.
    recipe.addVar(recipe.crystal.G1.delta2,
                  name="Ni_Delta2",
                  value=DELTA2_I,
                  tag="d2")

    recipe.addVar(recipe.crystal.G1.qdamp,
                  name="Calib_Qdamp",
                  value=QDAMP_I,
                  tag="inst")

    recipe.addVar(recipe.crystal.G1.qbroad,
                  name="Calib_Qbroad",
                  value=QBROAD_I,
                  tag="inst")

    # Configure some additional fit variables pertaining to symmetry.
    # We can use the srfit function constrainAsSpaceGroup to constrain
    # the lattice and ADP parameters according to the Fm-3m space group.
    # First we establish the relevant parameters, then we cycle through
    # the parameters and activate and tag them.
    # We must explicitly set the ADP parameters, because in this case, the
    # CIF had no ADP data.
    from diffpy.srfit.structure import constrainAsSpaceGroup

    # Create the symmetry distinct parameter sets, and constrain them
    # in the generator.
    spacegroupparams = constrainAsSpaceGroup(recipe.crystal.G1.phase,
                                             space_group)

    # Loop over all the symmetry distinct lattice parameters and add
    # them to the recipe.
    # We give them unique names, and tag them with our choice of a relevant string.
    for par in spacegroupparams.latpars:
        recipe.addVar(par,
                      value=CUBICLAT_I,
                      name="fcc_Lat",
                      tag="lat")

    # Loop over all the symmetry distinct ADPs and add
    # them to the recipe.
    # We give them unique names, and tag them with our choice of a relevant string.
    for par in spacegroupparams.adppars:
        recipe.addVar(par,
                      value=UISO_I,
                      name="fcc_ADP",
                      tag="adp")

    # Tell the Fit Recipe we want to write the maximum amount of
    # information to the terminal during fitting.
    recipe.fithooks[0].verbose = 3

    # During the optimization, we fix and free parameters sequentially
    # as you would in PDFgui. This leads to more stability in the refinement.
    # We first fix all variables. "all" is a tag which incorporates
    # every parameter.
    recipe.fix("all")

    # We then run a fit using the SciPy function "least_squares" which
    # takes as its arguments the function to be optimized, here recipe.residual,
    # as well as initial values for the fitted parameters, provided by
    # recipe.values. The x_scale="jac" argument is an optional argument
    # that provides for a bit more stability in the refinement.
    # "least_squares" is a bit more robust than "leastsq,"
    # which is another optimization function provided by SciPy.
    # "least_squares" supports bounds on refined parameters,
    #  while "leastsq" does not.

    refine_params = ["scale", "lat", "adp", "d2", "all"]

    for params in refine_params:
        recipe.free(params)
        print(f"\n****\nFitting {recipe.getNames()} against "
              f"{GR_NAME} with {CIF_NAME}\n")
        least_squares(recipe.residual, recipe.values, x_scale="jac")

    # We use the savetxt method of the profile to write a text file
    # containing the measured and fitted PDF to disk.
    # The file is named based on the basename we created earlier, and
    # written to the fitdir directory.
    profile = recipe.crystal.profile
    profile.savetxt(fitdir / (basename + ".fit"))

    # We use the FitResults function to parse out the results from
    # the optimized Fit Recipe.
    res = FitResults(recipe)

    # We print these results to the terminal.
    res.printResults()

    # We grab the fit Rw
    rw = res.rw

    # We use the saveResults method of FitResults to write a text file
    # containing the fitted parameters and fit quality indices to disk.
    # The file is named based on the basename we created earlier, and
    # written to the resdir directory.
    header = "crystal_HF.\n"
    res.saveResults(resdir / (basename + ".res"), header=header)

    # We use the plotresults function we created earlier to make a plot of
    # the measured, calculated, and difference curves. We show this
    # as an interactive window and then write a pdf file to disk.
    # The file is named based on the basename we created earlier, and
    # written to the figdir directory.
    plotresults(recipe, figdir / basename)

    # Let make a dictionary to hold our results. This way make reloading the
    # fit parameters easier later
    refined_dict = dict()

    refined_dict['rw'] = rw.item()

    recipe.free("all")
    # We loop over the variable names, the variable values, and the variable uncertainties (esd)
    for name, val, unc in zip(res.varnames, res.varvals, res.varunc):
        # We store the refined value for this variable using the "value" key.
        # We use the ".item()" method because "res.varvals" exist as
        # numpy.float64 objects, and we want them as regular python floats.
        if name not in refined_dict:
            refined_dict[name] = dict()
        refined_dict[name]["value"] = val.item()
        refined_dict[name]["uncert"] = unc.item()

    # Finally, let's write our dictionary to a yaml file!
    with open(basename + ".yml", 'w') as outfile:
        yaml.safe_dump(refined_dict, outfile)

    # End of function


if __name__ == "__main__":
    main()

# End of file
