import argparse
import os
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from coffea.util import load

# Set CMS style
plt.style.use(hep.style.CMS)

def plot_histogram(hist_data, title, xlabel, output_path):
    """Plot a single histogram and save to file"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot the histogram
    hep.histplot(hist_data, ax=ax, label="ttbar SemiLeptonic")

    # Compute statistics
    counts = hist_data.values()
    edges = hist_data.axes[0].edges
    centers = 0.5 * (edges[:-1] + edges[1:])
    total_entries = int(np.sum(counts))
    
    # Weighted mean and std
    if total_entries > 0:
        mean = np.sum(centers * counts) / total_entries
        std = np.sqrt(np.sum(counts * (centers - mean)**2) / total_entries)
    else:
        mean, std = 0, 0

    # Try to fit with ROOT (Gaussian) and overlay
    try:
        import ROOT
        from array import array
        # prepare arrays for TGraphErrors
        n = len(centers)
        xa = array('d', centers.tolist())
        ya = array('d', counts.tolist())
        ex = array('d', [0.0] * n)
        ey = array('d', [np.sqrt(c) if c>=0 else 0.0 for c in counts.tolist()])
        g = ROOT.TGraphErrors(n, xa, ya, ex, ey)
        f1 = ROOT.TF1('f1', 'gaus', centers[0], centers[-1])
        g.Fit(f1, 'Q')
        # evaluate fit for overlay
        xs = np.linspace(centers[0], centers[-1], 300)
        amp = f1.GetParameter(0)
        fit_mean = f1.GetParameter(1)
        sigma = f1.GetParameter(2)
        fit_vals = amp * np.exp(-0.5 * ((xs - fit_mean) / sigma) ** 2)
        ax.plot(xs, fit_vals, color='r', lw=2, label=f'Gauss fit: $\\mu$={fit_mean:.1f}, $\\sigma$={sigma:.1f}')
    except Exception as e:
        # If ROOT not available or fit failed, skip fitting
        print('Fit skipped (ROOT not available or fit failed):', e)
    
    # Add statistics box
    stats_text = f"Entries: {total_entries}\\nMean: {mean:.2f}\\nStd Dev: {std:.2f}"
    ax.text(0.65, 0.95, stats_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Events", fontsize=14)
    ax.set_title(title, fontsize=16, pad=20)
    ax.legend(fontsize=12, loc='upper left')
    ax.grid(True, alpha=0.3)
    
    # Add CMS label
    hep.cms.label(loc=0, data=False, llabel="Simulation", ax=ax)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot ttbar reconstruction test histograms")
    parser.add_argument('-i', '--input', type=str, required=True, 
                        help='Input .coffea file path')
    parser.add_argument('-o', '--outdir', type=str, default='plots', 
                        help='Output directory for plots')
    args = parser.parse_args()
    
    # Create output directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    # Load the coffea file
    print(f"Loading data from {args.input}")
    data = load(args.input)
    
    # Extract the dataset name (should be single dataset)
    dataset_names = list(data.keys())
    if len(dataset_names) == 0:
        print("Error: No datasets found in file")
        return
    
    dataset_name = dataset_names[0]
    print(f"Processing dataset: {dataset_name}")
    
    histos = data[dataset_name]['histos']
    
    # Define histogram configurations (name, title, xlabel)
    hist_configs = [
        ("gen_top_mass", "Generator-Level Top Quark Mass", "Gen Top Mass [GeV]"),
        ("gen_antitop_mass", "Generator-Level Anti-Top Quark Mass", "Gen Anti-Top Mass [GeV]"),
        ("gen_ttbar_mass", r"Generator-Level $t\bar{t}$ System Mass", r"Gen $t\bar{t}$ Mass [GeV]"),
        ("reco_top_mass", "Reconstructed Top Quark Mass", "Reco Top Mass [GeV]"),
        ("reco_antitop_mass", "Reconstructed Anti-Top Quark Mass", "Reco Anti-Top Mass [GeV]"),
        ("reco_ttbar_mass", r"Reconstructed $t\bar{t}$ System Mass", r"Reco $t\bar{t}$ Mass [GeV]"),
        ("delta_top_mass", r"$\Delta$ Top Mass (Gen - Reco)", r"$\Delta$ Top Mass [GeV]"),
        ("delta_antitop_mass", r"$\Delta$ Anti-Top Mass (Gen - Reco)", r"$\Delta$ Anti-Top Mass [GeV]"),
        ("delta_ttbar_mass", r"$\Delta$ $t\bar{t}$ Mass (Gen - Reco)", r"$\Delta$ $t\bar{t}$ Mass [GeV]"),
    ]
    
    # Plot each histogram
    for hist_name, title, xlabel in hist_configs:
        if hist_name in histos:
            output_path = os.path.join(args.outdir, f"{hist_name}.png")
            plot_histogram(histos[hist_name], title, xlabel, output_path)
        else:
            print(f"Warning: Histogram '{hist_name}' not found in data")
    
    print(f"\nAll plots saved to {args.outdir}/")


if __name__ == '__main__':
    main()
