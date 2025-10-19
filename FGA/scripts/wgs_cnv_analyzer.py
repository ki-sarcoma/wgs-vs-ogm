from pathlib import Path
from pandas import DataFrame, read_csv


class WGSCNVAnalyzer:
    def __init__(
        self,
        data_directory: Path = None,
        results_directory: Path = None,
        genome_size: int = 2_881_033_286,  # Autosomal hg19
        keep_autosomes_only: bool = True,
        min_size: int = None,
        log2_threshold: float = None,
    ):
        if data_directory is None:
            script_dir = Path(__file__).resolve().parent
            self.data_directory = Path(script_dir.parent, "data", "WGS").resolve()
        else:
            self.data_directory = Path(data_directory).resolve()

        if results_directory is None:
            script_dir = Path(__file__).resolve().parent
            self.results_directory = Path(script_dir.parent, "results").resolve()
        else:
            self.results_directory = Path(results_directory).resolve()

        self.genome_size = genome_size
        self.min_size = min_size
        self.log2_threshold = log2_threshold
        self.keep_autosomes_only = keep_autosomes_only

    def read_cnv_file(self, filepath) -> DataFrame:
        """Read a single CNV file and return a dataframe."""
        return read_csv(filepath, sep="\t", header=0)

    def get_filtered_cnvs_by_autosomes(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs to keep only autosomes (chromosomes 1-22)."""
        cnvs_filtered = cnvs[cnvs["chromosome"].astype(str).str.isdigit()]
        return cnvs_filtered[cnvs_filtered["chromosome"].astype(int).between(1, 22)]

    def get_filtered_cnvs_by_size(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on minimum size threshold."""
        cnvs["size"] = cnvs["end"] - cnvs["start"] + 1
        return cnvs[cnvs["size"] >= self.min_size]

    def get_filtered_cnvs_by_log2_ratio(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on log2 ratio threshold."""
        return cnvs[cnvs["log2"].abs() >= self.log2_threshold]

    def get_filtered_cnvs(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on log2 ratio, size, and mask fraction."""
        if self.keep_autosomes_only:
            cnvs: DataFrame = self.get_filtered_cnvs_by_autosomes(cnvs)
        if self.min_size is not None:
            cnvs: DataFrame = self.get_filtered_cnvs_by_size(cnvs)
        if self.log2_threshold is not None:
            cnvs: DataFrame = self.get_filtered_cnvs_by_log2_ratio(cnvs)

        return cnvs

    def calculate_fga_for_sample(self, cnvs: DataFrame) -> float:
        """Calculate fraction of genome altered."""
        altered_bp: int = cnvs["size"].sum(skipna=True)
        fga: float = altered_bp / self.genome_size
        return round(fga, 4)

    def save_fga_summary(self, fga_summary: DataFrame, output_csv: str) -> None:
        """Save FGA summary to a CSV file."""
        fga_summary_path = Path(self.results_directory, output_csv)
        fga_summary.to_csv(fga_summary_path, index=False)

    def compute_fga(self, output_csv) -> None:
        """Process all CNV files in a directory and save FGA summary."""
        fga_summary = []

        for filepath in self.data_directory.glob("*.cns"):
            sample_id: str = filepath.stem.split("-")[2][1:-2]
            cnvs: DataFrame = self.read_cnv_file(filepath)

            # Filter CNVs
            filtered_cnvs = self.get_filtered_cnvs(cnvs)

            # Compute FGA
            fga: float = self.calculate_fga_for_sample(filtered_cnvs)
            fga_summary.append({"SampleID": sample_id, "FGA": fga})

        self.save_fga_summary(DataFrame(fga_summary), output_csv)


if __name__ == "__main__":
    ogm_cnv_analyzer = WGSCNVAnalyzer()
    ogm_cnv_analyzer.compute_fga(output_csv="wgs_fga.csv")
