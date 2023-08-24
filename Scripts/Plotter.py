from coffea.util import load
import matplotlib.pyplot as plt
import mplhep as hep
import os

# Output Directory
outputDir = '/home/mukund2/Projects/CoffeaTrial/outputs'


# load coffea
output = load("output.coffea")
plt.style.use(hep.style.ROOT)

#################################################
# get plot from .coffea for Electron_pt
fig1, ax = plt.subplots()
hep.histplot(
    output["Run2016B_ver1_el"]["electron_pt"], label="electron pT", ax=ax
)
plt.savefig(os.path.join(outputDir,"Electron_pT.pdf"))
# stack the MC and overlay data with errorbar
fig2, ax2 = plt.subplots()
hep.histplot(
    [output[s]["electron_pt"] for s in output.keys() if "Run" not in s],
    histtype="fill",
    stack=True,
    label=[s for s in output.keys() if "Run" not in s],
    ax=ax2,
)

hep.histplot(
    [output[s]["electron_pt"] for s in output.keys() if "Run" in s],
    histtype="errorbar",
    label=[s for s in output.keys() if "Run" in s],
    ax=ax2,
    color="k",
)
plt.legend()
plt.savefig(os.path.join(outputDir,"Electron_pT_stacked.pdf"), dpi=300, bbox_inches="tight")

#################################################
# get plot from .coffea for Electron_eta
fig3, ax = plt.subplots()
hep.histplot(
    output["Run2016B_ver1_el"]["electron_eta"], label="electron eta", ax=ax
)
plt.savefig(os.path.join(outputDir,"Electron_eta.pdf"))
# stack the MC and overlay data with errorbar
fig4, ax2 = plt.subplots()
hep.histplot(
    [output[s]["electron_eta"] for s in output.keys() if "Run" not in s],
    histtype="fill",
    stack=True,
    label=[s for s in output.keys() if "Run" not in s],
    ax=ax2,
)

hep.histplot(
    [output[s]["electron_eta"] for s in output.keys() if "Run" in s],
    histtype="errorbar",
    label=[s for s in output.keys() if "Run" in s],
    ax=ax2,
    color="k",
)
plt.legend()
plt.savefig(os.path.join(outputDir,"Electron_eta_stacked.pdf"), dpi=300, bbox_inches="tight")


from matplotlib.backends.backend_pdf import PdfPages
figures = [fig1, fig2, fig3, fig4]
output_file = os.path.join(outputDir,"combined_plots.pdf")
with PdfPages(output_file) as pdf:
    for fig in figures:
        pdf.savefig(fig)  # Save the figure as a page in the PDF
