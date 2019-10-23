# Tutorial 1, Refining PDF from Ni to characterize instrument
#
# This Diffpy-CMI script will carry out a structural refinement of a measured
# PDF from nickel.  Nickel is a standard, and allows us to understand how the
# instrument affects our data, so that we can correct for it when we
# look at an unknown sample (in tutorial 2).
#
# Import packages that we will need.
from pathlib import Path
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

# If you'd like to run this in parallel (multi-CPUs):
# Ensure that the "psutil" package is installed in your environment
# Switch this to "True"
RUN_PARALLEL = False


######## Functions that will carry out the refinement ##################
# Make the recipe that the fit will follow.
# This Fit recipe object contains the PDF data, information on all the
# structures and all relevant details necessary to run the fit.
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

    ######## PDF Generator Section ##################
    # Create a PDF Generator object for a periodic structure model.
    # Here we name it "G1" and we give it the structure object.
    # This Generator will later compute the model PDF for the structure
    # object we provide it here.
    generator_crystal1 = PDFGenerator("G1")
    generator_crystal1.setStructure(structure, periodic=True)

    # Initialize the instrument parameters, Q_damp and Q_broad, and
    # assign Q_max and Q_min, all part of the PDF Generator object.
    # It's possible that the PDFParse function we used above
    # already parsed out ths information, but in case it didn't, we set it
    # explicitly again here.
    # All parameter objects can have their value assigned using the
    # below ".value = " syntax.
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
    # Create a Fit Contribution object, and name it "crystal."
    # We then give the PDF Generator object we created just above
    # to this Fit Contribution object. The Fit Contribution holds
    # the equation used to fit the PDF.
    contribution = FitContribution("crystal")
    contribution.addProfileGenerator(generator_crystal1)

    # Set an equation, within the Fit Contribution, based on your PDF
    # Generators. Here we simply have one Generator, G1, and a scale variable,
    # s1. Using this structure is a very flexible way of adding additional
    # Generators (ie. multiple structural phases), experimental Profiles,
    # PDF characteristic functions (ie. shape envelopes), and more.
    contribution.setEquation("s1*G1")


    # Set the experimental profile, within the Fit Contribution object,
    # to the Profile object we created earlier.
    contribution.setProfile(profile, xname="r")

    ######## Recipe Section ##################
    # Create the Fit Recipe object that holds all the details of the fit,
    # defined in the lines above. We give the Fit Recipe the Fit
    # Contribution we created earlier.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    # Add a variable to the Fit Recipe object, initialize the variables
    # with some value, and tag it with a string. Here we add the scale
    # parameter from the Fit Contribution. The ".addVar" method can be
    # used generally to add variables to the Fit Recipe.
    recipe.addVar(contribution.s1, SCALE_I, tag="scale")

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
    spacegroupparams = constrainAsSpaceGroup(generator_crystal1.phase,
                                             space_group)

    # Loop over all the symmetry distinct lattice parameters and add
    # them to the recipe.
    # We give them unique names, and tag them with our choice of a relevant string.
    for par in spacegroupparams.latpars:
        recipe.addVar(par,
                      value=CUBICLAT_I,
                      fixed=False,
                      name="fcc_Lat",
                      tag="lat")

    # Loop over all the symmetry distinct ADPs and add
    # them to the recipe.
    # We give them unique names, and tag them with our choice of a relevant string.
    for par in spacegroupparams.adppars:
        recipe.addVar(par,
                      value=UISO_I,
                      fixed=False,
                      name="fcc_ADP",
                      tag="adp")

    # Add delta and instrumental parameters to Fit Recipe.
    # These parameters are contained as part of the PDF Generator object
    # and initialized with values as defined in the opening of the script.
    # We give them unique names, and tag them with our choice of relevant strings.
    recipe.addVar(generator_crystal1.delta2,
                  name="Ni_Delta2",
                  value=DELTA2_I,
                  tag="d2")

    recipe.addVar(generator_crystal1.qdamp,
                  fixed=False,
                  name="Calib_Qdamp",
                  tag="inst")

    recipe.addVar(generator_crystal1.qbroad,
                  fixed=False,
                  name="Calib_Qbroad",
                  tag="inst")

    # Return the Fit Recipe object to be optimized
    return recipe

    # End of function

def plotresults(recipe, figname):
    """
    Creates plots of the fitted PDF and residual, displays them, and
    then writes them to disk as *.pdf files.

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
