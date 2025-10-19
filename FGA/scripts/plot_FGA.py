import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

script_dir = Path(__file__).resolve().parent
results_directory = Path(script_dir.parent, "results").resolve()

# File paths
ogm_file = Path(results_directory, "ogm_fga.csv")
wgs_file = Path(results_directory, "wgs_fga.csv")

# Load data
df_ogm = pd.read_csv(ogm_file)
df_wgs = pd.read_csv(wgs_file)

# Merge on SampleID
df = pd.merge(df_ogm, df_wgs, on="SampleID", suffixes=("_OGM", "_WGS"))

# Compute Pearson correlation coefficient
r, p_value = pearsonr(df["FGA_OGM"], df["FGA_WGS"])
print(f"Pearson correlation: r = {r:.3f}, p-value = {p_value:.3e}")

# Scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(df["FGA_OGM"], df["FGA_WGS"], color="dodgerblue", edgecolor="k", alpha=0.7)
plt.plot([0, 1], [0, 1], color="red", linestyle="--", label="y = x")
plt.xlabel("FGA (OGM)")
plt.ylabel("FGA (WGS)")
plt.title("Comparison of FGA: OGM vs WGS")
plt.grid(True)

# Add correlation text to the plot
plt.text(
    0.05,
    0.9,
    f"r = {r:.2f}",
    transform=plt.gca().transAxes,
    fontsize=12,
    bbox=dict(facecolor="white", alpha=0.5),
)

plt.legend()
plt.tight_layout()

plot_file = Path(results_directory, "fga_comparison.png")
plt.savefig(plot_file, dpi=300)

plt.show()
