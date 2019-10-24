#!/usr/bin/env python
###################################
#                                 #
# File coded by: Robert J. Koch   #
#                                 #
###################################
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import Profile
from diffpy.srfit.pdf import PDFParser, PDFGenerator
from diffpy.structure.parsers import getParser

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
    figname = Path(figname)
    r = recipe.crystal.profile.x

    g = recipe.crystal.profile.y
    gcalc = recipe.crystal.profile.ycalc
    diffzero = -0.65 * max(g) * np.ones_like(g)
    diff = g - gcalc + diffzero

    fig, ax1 = plt.subplots(1, 1)

    ax1.plot(r, g, ls="None",
             marker="o", ms=5, mew=0.2,
             mfc="None", label="G(r) Data", markeredgecolor="blue")

    ax1.plot(r, gcalc, lw=1.3, label="G(r) Fit",color="red")

    ax1.plot(r, diff, lw=1.2, label="G(r) diff",c="green")

    ax1.plot(r, diffzero, lw=1.0, ls="--", c="black")

    ax1.set_xlabel("r ($\mathrm{\AA}$)")
    ax1.set_ylabel("G ($\mathrm{\AA}$$^{-2}$)")
    ax1.tick_params(axis="both",
                    which="major",
                    top=True,
                    right=True)

    ax1.set_xlim(r.min(), r.max())
    ax1.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(figname.parent / (figname.name + ".pdf"),
                format="pdf")

    # End of function


######## Functions that will carry out the refinement ##################
# Make the recipe that the fit will follow.
# This Fit recipe object contains the PDF data, information on all the
# structures and all relevant details necessary to run the fit.
def makerecipe(structure_file,
               data_file):
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


    # Return the Fit Recipe object to be optimized
    return recipe

    # End of function
