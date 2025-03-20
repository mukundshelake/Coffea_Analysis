import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load CSV file
csvFile = "/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/Scripts/BDTImplementation/Binary/QQ_vs_nonQQ/outputs/Feb10_3/scores_binary.csv"
df = pd.read_csv(csvFile)  # Change to your actual file path

# Separate BDT scores based on class
qq_scores = df[df["y_actual"] == 0]["score"]
nonqq_scores = df[df["y_actual"] == 1]["score"]

# Print lengths of qq_scores and nonqq_scores
print(f"Length of qq_scores: {len(qq_scores)}")
print(f"Length of nonqq_scores: {len(nonqq_scores)}")

# Plot the histograms
plt.figure(figsize=(8, 6))
plt.hist(qq_scores, bins=30, histtype='step', label="qq", color="blue")
plt.hist(nonqq_scores, bins=30, histtype='step', label="non-qq", color="red")

# Labels and legend
plt.xlabel("BDT Score")
plt.ylabel("Event Count")
plt.title("BDT Score Distribution")
plt.legend()
plt.grid(True)

# Save plot to file
output_file = csvFile.replace("scores_binary.csv", "scoreHist.png")
plt.savefig(output_file)
print(f"Plot saved to {output_file}")

