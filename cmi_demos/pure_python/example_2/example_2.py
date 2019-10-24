#!/usr/bin/env python
###################################
#                                 #
# File coded by: Robert J. Koch   #
#                                 #
###################################
#
# Example 2, Refining PDF from nanocrystalline platinum to obtain
# an estimate for the nanoparticle (NP) size, as well as the atomic structure.
#
# This Diffpy-CMI script will carry out a structural refinement of a measured
# PDF from nanocrystalline platinum.  It is explicitly uses the instrumental parameters
# found in tutorial 1, and requires that the output files of this tutorial
# are present in a specific location.
#
# Comments in this Example will be less verbose than in Example 1.
#
# Import packages that we will need
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
DPATH = Path("data")

# Give an identifying name for the refinement, similar
# to what you would name a fit tree in PDFGui.
FIT_ID = "Fit_Pt_NP"

# Specify the names of the input PDF and cif files.
GR_NAME = "Pt-nanoparticles.gr"
CIF_NAME = "Pt.cif"

######## Experimental PDF Config ######################
# Specify the min, max, and step r-values of the PDF (that we want to fit over)
# also, specify the Q_max and Q_min values used to reduce the PDF.
PDF_RMIN = 1.5
PDF_RMAX = 50
PDF_RSTEP = 0.01
QMAX = 25
QMIN = 0.1

########PDF initialize refinable variables #############
# We explicitly specify the lattice parameters, scale,
# isotropic thermal parameters, and a correlated motion parameter.
CUBICLAT_I = 3.9
SCALE_I = 0.6
UISO_I = 0.01
DELTA2_I = 4

# For the NP case, we also provide an initial guess
# for the average crystallite size, in Angstroms.
PSIZE_I = 40

# First, let's read the fit results from the Ni fit.
# We parse out the refined values of Q_damp and Q_broad,
# instrumental parameters which will be fixed in this fit.

# We specify where to look for the fit results files, as well as
# a name of the file.
STANDARD_DIR = Path().absolute().parent / "example_1"
NIBASENAME = "Fit_Ni_Bulk"

STANDARD_RES_FILE = (STANDARD_DIR / (NIBASENAME + ".yml"))

# First, check of the file even exists
if STANDARD_RES_FILE.exists():

    # If it exists, let's open it and load the dictionary
    with open(STANDARD_RES_FILE, 'r') as infile:
        ni_dict = yaml.safe_load(infile)

    # If we find either of these strings in the dictionary, we load it
    if "Calib_Qbroad" in ni_dict:
        # ... we grab the value.
        QBROAD_I = ni_dict["Calib_Qbroad"]["value"]

    if "Calib_Qdamp" in ni_dict:
        # ... we grab the value.
        QDAMP_I = ni_dict["Calib_Qdamp"]["value"]

    # If we don't find these strings, something is wrong...            
    if "Calib_Qbroad" not in ni_dict or "Calib_Qdamp" not in ni_dict:
        # ... so we print a warning to the terminal...
        print(f"\nWe did not find instrumental parameters in"
              f"{STANDARD_RES_FILE}. Please make sure tutorial 1 has been")

        # ...and we exit the script with code "1," meaning abnormal termination
        sys.exit(1)

# If the fit result file does not exist, something is wrong.
else:

    # We print a warning to the terminal...
    print("\nTutorial 1, which refines instrumental parameters, "
          "does not appear to have been run.\nThese instrumental "
          "parameters are necessary to run this tutorial.\n"
          "Please run tutorial 1 first")

    # ...and we exit the script with code "1," meaning abnormal termination
    sys.exit(1)


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

    # Initialize the Fit Recipe by giving it this diffpy structure
    # as well as the path to the data file.

    p_cif = getParser('cif')
    structure = p_cif.parseFile(str(cif_file))
    space_group = p_cif.spacegroup.short_name

    # Initialize the Fit Recipe by giving it this diffpy structure
    # as well as the path to the data file.
    recipe = makerecipe(cif_file, data)

    # Let's set the calculation range!
    recipe.crystal.profile.setCalculationRange(xmin=PDF_RMIN, xmax=PDF_RMAX, dx=PDF_RSTEP)

    # Add, initialize, and tag variables in the Fit Recipe object.
    # In this case we also add psize, which is the NP size.
    recipe.addVar(recipe.crystal.s1, SCALE_I, tag="scale")

    # Set an equation, based on your PDF generators. Here we add an extra layer
    # of complexity, incorporating "f" int our equation. This new term
    # incorporates damping to our PDF to model the effect of finite crystallite size.
    # In this case we use a function which models a spherical NP.
    from diffpy.srfit.pdf.characteristicfunctions import sphericalCF
    recipe.crystal.registerFunction(sphericalCF, name="f")
    recipe.crystal.setEquation("s1*G1*f")

    recipe.addVar(recipe.crystal.psize, PSIZE_I, tag="psize")

    # Initialize the instrument parameters, Q_damp and Q_broad, and
    # assign Q_max and Q_min.
    # Note, here we do not add the qdamp and qbroad parameters to the fit!!!
    # They are fixed here, because we refined them in the Ni standard fit!
    recipe.crystal.G1.qdamp.value = QDAMP_I
    recipe.crystal.G1.qbroad.value = QBROAD_I
    recipe.crystal.G1.setQmax(QMAX)
    recipe.crystal.G1.setQmin(QMIN)

    # Use the srfit function constrainAsSpaceGroup to constrain
    # the lattice and ADP parameters according to the Fm-3m space group.
    from diffpy.srfit.structure import constrainAsSpaceGroup
    spacegroupparams = constrainAsSpaceGroup(recipe.crystal.G1.phase,
                                             space_group)

    # Add and initialize delta, the lattice parameter, and a thermal parameter,
    # but not instrumental parameters to Fit Recipe.
    # The instrumental parameters will remain fixed at values obtained from
    # the Ni calibrant in our previous example. As we have not added them through
    # recipe.addVar, they cannot be refined.
    for par in spacegroupparams.latpars:
        recipe.addVar(par,
                      value=CUBICLAT_I,
                      name="fcc_Lat",
                      tag="lat")

    for par in spacegroupparams.adppars:
        recipe.addVar(par,
                      value=UISO_I,
                      name="fcc_ADP",
                      tag="adp")

    recipe.addVar(recipe.crystal.G1.delta2,
                  name="Pt_Delta2",
                  value=DELTA2_I,
                  tag="d2")

    # Tell the Fit Recipe we want to write the maximum amount of
    # information to the terminal during fitting.
    recipe.fithooks[0].verbose = 0

    refine_params = ["scale", "lat", "psize", "adp", "d2", "all"]

    recipe.fix("all")

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

    # We loop over the variable names, the variable values, and the variable uncertainties (esd)
    for name, val, unc in zip(res.varnames, res.varvals, res.varunc):
        # We store the refined value for this variable using the "value" key.
        # We use the ".item()" method because "res.varvals" exist as
        # numpy.float64 objects, and we want them as regular python floats.
        if name not in refined_dict:
            refined_dict[name] = dict()
        refined_dict[name]["value"] = val.item()
        refined_dict[name]["uncert"] = unc.item()

    with open(basename + ".yml", 'w') as outfile:
        yaml.safe_dump(refined_dict, outfile)

    # End of function


if __name__ == "__main__":
    main()

# End of file
