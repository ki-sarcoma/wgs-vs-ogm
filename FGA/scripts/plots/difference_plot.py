import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# Paths
script_dir = Path(__file__).resolve().parent
results_directory = Path(script_dir.parent.parent, "results").resolve()

ogm_file = Path(results_directory, "ogm_fga.csv")
wgs_file = Path(results_directory, "wgs_fga.csv")

# Load data
df_ogm = pd.read_csv(ogm_file)
df_wgs = pd.read_csv(wgs_file)

# Merge data on SampleID
df = pd.merge(df_ogm, df_wgs, on="SampleID", suffixes=("_OGM", "_WGS"))
n_samples = len(df)

# Bland-Altman calculations
mean_values = df[["FGA_OGM", "FGA_WGS"]].mean(axis=1)
diff_values = df["FGA_OGM"] - df["FGA_WGS"]
mean_diff = np.mean(diff_values)
std_diff = np.std(diff_values, ddof=1)

# Limits of agreement
loa_upper = mean_diff + 1.96 * std_diff
loa_lower = mean_diff - 1.96 * std_diff

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(
    mean_values, diff_values, color="dodgerblue", alpha=0.6, edgecolor="k", s=70
)

# Horizontal lines
plt.axhline(mean_diff, color="red", linestyle="dotted")
plt.axhline(loa_upper, color="green", linestyle="dotted")
plt.axhline(loa_lower, color="green", linestyle="dotted")

# Text labels at the right end of each line
x_text = max(mean_values) - 0.045
x_text_main = max(mean_values) - 0.02
y_text_offset = 0.003

# Mean difference
plt.text(
    x_text_main,
    mean_diff + y_text_offset,
    "Mean",
    color="red",
    va="bottom",
    fontsize=10,
)
plt.text(
    x_text_main,
    mean_diff - y_text_offset,
    f"{mean_diff:.3f}",
    color="red",
    va="top",
    fontsize=10,
)

# Upper limit
plt.text(
    x_text - 0.005,
    loa_upper + y_text_offset,
    "+1.96 SD",
    color="green",
    va="bottom",
    fontsize=10,
)
plt.text(
    x_text - 0.005,
    loa_upper - y_text_offset,
    f"{loa_upper:.3f}",
    color="green",
    va="top",
    fontsize=10,
)

# Lower limit
plt.text(
    x_text,
    loa_lower + y_text_offset,
    "-1.96 SD",
    color="green",
    va="bottom",
    fontsize=10,
)
plt.text(
    x_text,
    loa_lower - y_text_offset,
    f"{loa_lower:.3f}",
    color="green",
    va="top",
    fontsize=10,
)


# Add sample size
plt.text(
    0,
    0.157,
    f"n= {n_samples}",
    fontsize=10,
    bbox=dict(
        facecolor="white",
        boxstyle="round,rounding_size=0.2, pad=0.3",
        alpha=0.9,
        edgecolor="gray",
    ),
)

plt.xlabel("Mean FGA (OGM & WGS)", fontsize=12)
plt.ylabel("Difference in FGA (OGM minus WGS)", fontsize=12)
plt.tight_layout()

# Save plot
plot_file = Path(results_directory, "difference_plot_FGA.png")
plt.savefig(plot_file, dpi=300)

plt.show()
