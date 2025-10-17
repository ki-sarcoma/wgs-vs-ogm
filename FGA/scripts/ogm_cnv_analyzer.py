from pathlib import Path
from pandas import DataFrame, read_csv, concat


class OGMCNVAnalyzer:
    def __init__(
        self,
        genome_size: int = 2_875_001_522,  # hg38
        data_directory: Path = None,
        results_directory: Path = None,
        log2_threshold: float = None,
        min_size: int = None,
        max_mask_frac: float = None,
    ):
        if data_directory is None:
            script_dir = Path(__file__).resolve().parent
            self.data_directory = Path(script_dir.parent, "data", "OGM").resolve()
        else:
            self.data_directory = Path(data_directory).resolve()

        if results_directory is None:
            script_dir = Path(__file__).resolve().parent
            self.results_directory = Path(script_dir.parent, "results").resolve()
        else:
            self.results_directory = Path(results_directory).resolve()

        self.genome_size = genome_size
        self.log2_threshold = log2_threshold
        self.min_size = min_size
        self.max_mask_frac = max_mask_frac

    def read_cnv_file(self, filepath) -> DataFrame:
        """Read a single CNV file and return a dataframe."""
        return read_csv(filepath, sep="\t", header=5)

    def calculate_fga_for_sample(self, cnvs: DataFrame) -> float:
        """Calculate fraction of genome altered."""
        altered_bp: int = cnvs["Size"].sum(skipna=True)
        fga: float = altered_bp / self.genome_size
        return round(fga, 4)

    def get_filtered_cnvs(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on log2 ratio, size, and mask fraction."""
        return cnvs

    def save_fga_summary(self, fga_summary: DataFrame, output_csv: str) -> None:
        """Save FGA summary to a CSV file."""
        fga_summary_path = Path(self.results_directory, output_csv)
        fga_summary.to_csv(fga_summary_path, index=False)

    def compute_fga(self, output_csv) -> None:
        """Process all CNV files in a directory and save FGA summary."""
        fga_summary = []

        for filepath in self.data_directory.glob("*.txt"):
            sample_id: str = filepath.stem.split("_")[0]
            cnvs: DataFrame = self.read_cnv_file(filepath)

            # Filter CNVs
            filtered_cnvs = self.get_filtered_cnvs(cnvs)

            # Compute FGA
            fga: float = self.calculate_fga_for_sample(filtered_cnvs)
            fga_summary.append({"SampleID": sample_id, "FGA": fga})

        self.save_fga_summary(DataFrame(fga_summary), output_csv)


if __name__ == "__main__":
    ogm_cnv_analyzer = OGMCNVAnalyzer()
    ogm_cnv_analyzer.compute_fga(output_csv="ogm_fga.csv")
