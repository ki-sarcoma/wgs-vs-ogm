import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

# Paths
script_dir = Path(__file__).resolve().parent
results_directory = Path(script_dir.parent, "results").resolve()

ogm_file = Path(results_directory, "ogm_fga.csv")
wgs_file = Path(results_directory, "wgs_fga.csv")

# Load data
df_ogm = pd.read_csv(ogm_file)
df_wgs = pd.read_csv(wgs_file)

# Merge
df = pd.merge(df_ogm, df_wgs, on="SampleID", suffixes=("_OGM", "_WGS"))

# Pearson correlation
r, p_value = pearsonr(df["FGA_OGM"], df["FGA_WGS"])
n_samples = len(df)

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(
    df["FGA_OGM"],
    df["FGA_WGS"],
    color="dodgerblue",
    edgecolor="k",
    alpha=0.6,
    s=70,
)
line_range = [0, 1]
plt.plot(line_range, line_range, color="red", linestyle="dotted", label="y = x")

# Plot labels and grid
plt.xlabel("FGA (OGM)", fontsize=12)
plt.ylabel("FGA (WGS)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.4)

# Add correlation
plt.text(
    -0.03,
    1,
    f"n= {n_samples}",
    fontsize=10,
    bbox=dict(
        facecolor="white",
        boxstyle="round,rounding_size=0.2, pad=0.3",
        alpha=0.3,
        edgecolor="gray",
    ),
)

# Add sample size
plt.text(
    -0.03,
    0.94,
    f"Pearson r = {r:.2f}",
    fontsize=10,
    bbox=dict(
        facecolor="white",
        boxstyle="round,rounding_size=0.2, pad=0.3",
        alpha=0.3,
        edgecolor="gray",
    ),
)

plt.legend()
plt.tight_layout()

plot_file = Path(results_directory, "fga_comparison.png")
plt.savefig(plot_file, dpi=300)

plt.show()
