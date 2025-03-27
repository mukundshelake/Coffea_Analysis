import ROOT

# Create histograms for different categories
hist_data = ROOT.TH1F("data", "Data;Tag e p_{T} (GeV);Events", 50, 0, 400)
hist_tt_prompt = ROOT.TH1F("tt_prompt", "t#bar{t} prompt #mu", 50, 0, 400)
hist_tt_tau = ROOT.TH1F("tt_tau", "t#bar{t} #tau #rightarrow #mu", 50, 0, 400)
hist_other_top = ROOT.TH1F("other_top", "Other Top", 50, 0, 400)

# Fill histograms with random data
for _ in range(1000):
    hist_data.Fill(ROOT.gRandom.Gaus(200, 50))
    hist_tt_prompt.Fill(ROOT.gRandom.Gaus(150, 40))
    hist_tt_tau.Fill(ROOT.gRandom.Gaus(100, 30))
    hist_other_top.Fill(ROOT.gRandom.Gaus(50, 20))

# Style histograms
hist_data.SetMarkerStyle(20)
hist_data.SetMarkerSize(1)
hist_data.SetLineColor(ROOT.kBlack)

hist_tt_prompt.SetFillColor(ROOT.kCyan)
hist_tt_tau.SetFillColor(ROOT.kMagenta)
hist_other_top.SetFillColor(ROOT.kYellow)

# Disable the statistics box for all histograms
hist_data.SetStats(False)
hist_tt_prompt.SetStats(False)
hist_tt_tau.SetStats(False)
hist_other_top.SetStats(False)

# Create a stack for the backgrounds
stack = ROOT.THStack("stack", "Signal region;Tag e p_{T} (GeV);Events")
stack.Add(hist_tt_prompt)
stack.Add(hist_tt_tau)
stack.Add(hist_other_top)

# Create a canvas
canvas = ROOT.TCanvas("canvas", "Histogram Canvas", 800, 800)

# Upper pad for the main plot
pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad1.SetBottomMargin(0)  # Remove space between pads
pad1.Draw()
pad1.cd()
# Set y-axis to log scale
pad1.SetLogy()
# Draw the stack and data
stack.Draw("HIST")
hist_data.Draw("E SAME")

# Add a legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.SetBorderSize(0)
# legend.SetNColumns(2)
legend.AddEntry(hist_data, "Data", "lep")
legend.AddEntry(hist_tt_prompt, "t#bar{t} prompt #mu", "f")
legend.AddEntry(hist_tt_tau, "t#bar{t} #tau #rightarrow #mu", "f")
legend.AddEntry(hist_other_top, "Other Top", "f")
legend.Draw()

# Add "CMS Preliminary" text
cms_text = ROOT.TLatex()
cms_text.SetNDC()  # Use normalized device coordinates
cms_text.SetTextSize(0.04)
cms_text.SetTextFont(42)
cms_text.DrawLatex(0.1, 0.92, "CMS Preliminary")  # Position (x=0.1, y=0.92)

# Add luminosity text
lumi_text = ROOT.TLatex()
lumi_text.SetNDC()
lumi_text.SetTextSize(0.04)
lumi_text.SetTextFont(42)
lumi_text.DrawLatex(0.7, 0.92, "36.3 fb^{-1} (13 TeV, 2016)")  # Position (x=0.7, y=0.92)

# Go back to the main canvas and create a lower pad for the ratio plot
canvas.cd()
pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.3)
pad2.SetTopMargin(0)  # Remove space between pads
pad2.SetBottomMargin(0.3)
pad2.Draw()
pad2.cd()

# Create a ratio plot
ratio = hist_data.Clone("ratio")
stack_total = stack.GetStack().Last()
ratio.Divide(stack_total)
ratio.SetTitle(";Tag e p_{T} (GeV);Data/MC")
# ratio.GetYaxis().SetRangeUser(0.7, 1.3)
ratio.GetYaxis().SetNdivisions(505)
ratio.GetYaxis().SetTitleSize(0.1)
ratio.GetYaxis().SetTitleOffset(0.5)
ratio.GetXaxis().SetTitleSize(0.1)
ratio.GetXaxis().SetTitleOffset(1.0)
ratio.GetXaxis().SetLabelSize(0.1)
ratio.GetYaxis().SetLabelSize(0.1)
ratio.Draw("E")

# Create a graph for the uncertainty band
uncertainty_band = ROOT.TGraphAsymmErrors(stack_total)

# Loop over bins to set the uncertainty
for i in range(1, stack_total.GetNbinsX() + 1):
    x = stack_total.GetBinCenter(i)
    y = 1  # Center the band at ratio = 1
    error = stack_total.GetBinError(i) / stack_total.GetBinContent(i) if stack_total.GetBinContent(i) > 0 else 0
    uncertainty_band.SetPoint(i - 1, x, y)
    uncertainty_band.SetPointError(i - 1, 0, 0, error, error)

# Style the uncertainty band
uncertainty_band.SetFillColor(ROOT.kGray + 2)
uncertainty_band.SetFillStyle(3001)  # Semi-transparent fill
uncertainty_band.Draw("E2 SAME")

# Loop over all bins in the histogram
for i in range(1, hist_data.GetNbinsX() + 1):  # Bin index starts from 1
    bin_content = hist_data.GetBinContent(i)
    bin_error = hist_data.GetBinError(i)
    print(f"Bin {i}: Content = {bin_content}, Error = {bin_error}, Total Error Bar Length = {2 * bin_error}")

# Save the canvas
canvas.SaveAs("histogram.png")
