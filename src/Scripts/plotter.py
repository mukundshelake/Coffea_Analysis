from coffea.util import load
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

out = load("coffeaOutput.coffea")
ytProjection = out["DoubleMuon"]["DoubleMuon"]["yMatrix"].project("y_t")
ytbarProjection = out["UL2016preVFP_ttbarSemi_Sample"]["UL2016preVFP_ttbarSemi_Sample"]["yMatrix"].project("y_tbar")

with PdfPages('output_plots.pdf') as pdf:
    fig, ax = plt.subplots(figsize=(10, 6))
    ytProjection.plot(ax=ax)
    ytbarProjection.plot(ax=ax)
    # ax.set_xscale("log")
    ax.legend(title="Dimuon charge")
    pdf.savefig() 



