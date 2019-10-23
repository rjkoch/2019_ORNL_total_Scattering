# Tutorial 2, Refining PDF from nanocrystalline platinum to obtain
# an estimate for the nanoparticle (NP) size, as well as the atomic structure.
#
# This Diffpy-CMI script will carry out a structural refinement of a measured
# PDF from nanocrystalline platinum.  It is explicitly uses the instrumental parameters
# found in tutorial 1, and requires that the output files of this tutorial
# are present in a specific location.
#
# Comments in this tutorial will be less verbose than in tutorial 1. 
#
# Import packages that we will need
from pathlib import Path
import sys
import yaml

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.pdf import PDFParser, PDFGenerator
from diffpy.structure.parsers import getParser
from scipy.optimize import least_squares

############### Config ##############################
# Give a file path to where your pdf (.gr) and (.cif) files are located.
DPATH = Path("data")

# Give an identifying name for the refinement, similar
# to what you would name a fit tree in PDFGui.
FIT_ID = "Fit_Pt_NP"

# Specify the names of the input PDF and cif files.
GR_NAME = "Pt-nanoparticles.gr"
CIF_NAME = "Pt.cif"

# Specify the space-group.
SPACEGROUP = "Fm-3m"

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

# If you'd like to run this in parallel (multi-CPUs):
# Ensure that the "psutil" python package is installed in your environment
# Switch this to "True"
RUN_PARALLEL = False

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
        refined_dict = yaml.safe_load(infile)

    # If we find either of these strings in the dictionary, we load it
    if "Calib_Qbroad" in refined_dict:
        # ... we grab the value.
        QBROAD_I = refined_dict["Calib_Qbroad"]["value"]

    if "Calib_Qdamp" in refined_dict:
        # ... we grab the value.
        QDAMP_I = refined_dict["Calib_Qdamp"]["value"]

    # If we don't find these strings, something is wrong...            
    if "Calib_Qbroad" not in refined_dict or "Calib_Qdamp" not in refined_dict:
        
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


######## Functions that will carry out the refinement ##################
# Make the Fit Recipe object, similar to the Ni example.
def makerecipe(structure_file, data_file):
    """
    Basic function for creating and properly constraining a fit recipe.

    Parameters
    ----------
    structure_file : Path object or str
        Path to *.cif file, containing a structural model to use to fit the PDF data.
    data_file : Path object or str
        Path to data file containing PDF data to fit against.

    Returns
    -------
    recipe : FitRecipe object
        An initialized fit recipe object, ready for fitting.
    """
    ######## Profile Section ##################
    # Create a Profile object for the experimental dataset.
    # This handles all details about the dataset.
    # We also tell this profile the range and mesh of points in r-space.
    profile = Profile()
    parser = PDFParser()
    parser.parseFile(data_file)
    profile.loadParsedData(parser)
    profile.setCalculationRange(xmin=PDF_RMIN, xmax=PDF_RMAX, dx=PDF_RSTEP)



    p_cif = getParser('cif')
    structure = p_cif.parseFile(str(structure_file))
    space_group = p_cif.spacegroup.short_name

    # Make sure we allow for anisotropy
    structure.anisotropy = True

    ######## PDF Generator Section ##################
    # Create a PDF Generator object for a periodic structure model.
    generator_crystal1 = PDFGenerator("G1")
    generator_crystal1.setStructure(structure, periodic=True)

    # Initialize the instrument parameters, Q_damp and Q_broad, and
    # assign Q_max and Q_min.
    generator_crystal1.qdamp.value = QDAMP_I
    generator_crystal1.qbroad.value = QBROAD_I
    generator_crystal1.setQmax(QMAX)
    generator_crystal1.setQmin(QMIN)

    # If you have a multi-core computer (you probably do),
    # you can run your refinement in parallel!
    # This requires that you set "RUN_PARALLEL" to "True" above.
    # The psutil python package is also required for the bit of
    # code below, where we make sure not to overload your CPUs.
    if RUN_PARALLEL:
        import psutil
        import multiprocessing
        syst_cores = multiprocessing.cpu_count()
        cpu_percent = psutil.cpu_percent()
        avail_cores = np.floor((100 - cpu_percent) / (100.0 / syst_cores))
        ncpu = int(np.max([1, avail_cores]))
        generator_crystal1.parallel(ncpu)

    ######## Fit Contribution Section ##################
    # Create a Fit Contribution object.
    contribution = FitContribution("crystal")
    contribution.addProfileGenerator(generator_crystal1)

    # Set an equation, based on your PDF generators. Here we add an extra layer
    # of complexity, incorporating "f" int our equation. This new term
    # incorporates damping to our PDF to model the effect of finite crystallite size.
    # In this case we use a function which models a spherical NP.
    from diffpy.srfit.pdf.characteristicfunctions import sphericalCF
    contribution.registerFunction(sphericalCF, name="f")
    contribution.setEquation("s1*G1*f")

    # Set the Fit Contribution profile to the Profile object.
    contribution.setProfile(profile, xname="r")

    ######## Recipe Section ##################
    # Create the Fit Recipe object that holds all the details of the fit.
    recipe = FitRecipe()
    recipe.addContribution(contribution)


    # Add, initialize, and tag variables in the Fit Recipe object.
    # In this case we also add psize, which is the NP size.
    recipe.addVar(contribution.s1, SCALE_I, tag="scale")
    recipe.addVar(contribution.psize, PSIZE_I, tag="psize")

    # Use the srfit function constrainAsSpaceGroup to constrain
    # the lattice and ADP parameters according to the Fm-3m space group.
    from diffpy.srfit.structure import constrainAsSpaceGroup
    spacegroupparams = constrainAsSpaceGroup(generator_crystal1.phase,
                                             space_group)
    for par in spacegroupparams.latpars:
        recipe.addVar(par,
                      value=CUBICLAT_I,
                      fixed=False,
                      name="fcc_Lat",
                      tag="lat")

    for par in spacegroupparams.adppars:
        recipe.addVar(par,
                      value=UISO_I,
                      fixed=False,
                      name="fcc_ADP",
                      tag="adp")

    # Add delta, but not instrumental parameters to Fit Recipe.
    # The instrumental parameters will remain fixed at values obtained from
    # the Ni calibrant in our previous example. As we have not added them through
    # recipe.addVar, they cannot be refined.
    recipe.addVar(generator_crystal1.delta2,
                  name="Pt_Delta2",
                  value=DELTA2_I,
                  tag="d2")

    return recipe

    # End of function


def plotresults(recipe, figname):
    """
    Creates plots of the fitted PDF and residual, and writes them to disk
    as *.pdf files.

    Parameters
    ----------
    recipe :    The optimized Fit Recipe object containing the PDF data
                we wish to plot
    figname :   string, the location and name of the figure file to create

    Returns
    ----------
    None
    """
    r = recipe.crystal.profile.x

    g = recipe.crystal.profile.y
    gcalc = recipe.crystal.profile.ycalc
    diffzero = -0.65 * max(g) * np.ones_like(g)
    diff = g - gcalc + diffzero

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.style.use(str(Path().absolute().parent / "utils" / "billinge.mplstyle"))

    fig, ax1 = plt.subplots(1, 1)

    ax1.plot(r, g, ls="None",
             marker="o", ms=5, mew=0.2,
             mfc="None", label="G(r) Data")

    ax1.plot(r, gcalc, lw=1.3, label="G(r) Fit")

    ax1.plot(r, diff, lw=1.2, label="G(r) diff")

    ax1.plot(r, diffzero, lw=1.0, ls="--", c="black")

    ax1.set_xlabel("r($\mathrm{\AA}$)")
    ax1.set_ylabel("G($\mathrm{\AA}$$^{-2}$)")
    ax1.tick_params(axis="both",
                    which="major",
                    top=True,
                    right=True)

    ax1.set_xlim(PDF_RMIN, PDF_RMAX)
    ax1.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(figname.parent / (figname.name + ".pdf"),
                format="pdf")

    # End of function


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
    print(basename)

    # Establish the full location of the data.
    data = DPATH / GR_NAME

    # Establish the location of the cif file with the structure of interest
    # and load it into a diffpy structure object.
    strudir = DPATH
    cif_file = strudir / CIF_NAME

    # Initialize the Fit Recipe by giving it this diffpy structure
    # as well as the path to the data file.
    recipe = makerecipe(cif_file, data)
    recipe.fix("all")

    # Tell the Fit Recipe we want to write the maximum amount of
    # information to the terminal during fitting.
    recipe.fithooks[0].verbose = 3

    refine_params = ["scale", "lat", "psize", "adp", "d2", "all"]

    for params in refine_params:
        recipe.free(params)
        print(f"\n****\nFitting {recipe.getNames()} against "
              f"{GR_NAME} with {CIF_NAME}\n")
        input("\nPress enter to continue...")
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
    # plt.ion()

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

    with open(basename + ".yml", 'w') as outfile:
        yaml.safe_dump(refined_dict, outfile)


    # End of function


if __name__ == "__main__":
    main()

# End of file
