from pathlib import Path
from pandas import DataFrame, Series, read_csv
from numpy import log2, maximum


class CNVAnalyzer:
    def __init__(
        self,
        data_directory: Path = None,
        results_directory: Path = None,
        genome_size: int = None,
        keep_autosomes_only: bool = True,
        min_size: int = None,
        log2_threshold: float = None,
    ):
        if data_directory is None:
            self.data_directory = self.get_data_directory()
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

    def get_data_directory(self) -> Path:
        """Get the data directory path."""
        pass

    def get_sample_id_from_path(self, filepath: Path) -> str:
        """Extract sample ID from file path."""
        pass

    def read_cnv_file(self, filepath) -> DataFrame:
        """Read a single CNV file and return a dataframe."""
        return read_csv(filepath, sep="\t", header=0)

    def get_sample_purity(self, sample_id: str) -> float:
        """Retrieve sample purity from a metadata file."""
        tumor_fraction_path = Path(self.data_directory.parent, "tumor_fraction.csv")
        tumor_fraction: DataFrame = read_csv(tumor_fraction_path, sep=",")
        result: Series = tumor_fraction.loc[
            tumor_fraction["SampleID"].astype(str) == sample_id, "Tumor fraction"
        ]

        if result.empty:
            raise ValueError(f"SampleID {sample_id} not found in tumor_fraction.csv")

        return result.values[0]

    def get_purity_corrected_log2_ratios_cnvs(
        self, sample_id: str, cnvs: DataFrame
    ) -> DataFrame:
        """Calculate purity-corrected log2 ratios for CNVs."""
        sample_purity: float = self.get_sample_purity(sample_id)
        copy_number_obs = 2 * (2 ** cnvs["log2"])
        copy_number_tumor = maximum(
            (copy_number_obs - 2 * (1 - sample_purity)) / sample_purity, 0.01
        )
        cnvs["log2_corrected"] = round(log2(copy_number_tumor / 2), 2)
        return cnvs

    def get_filtered_cnvs_by_autosomes(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs to keep only autosomes (chromosomes 1-22)."""
        cnvs_filtered = cnvs[cnvs["chromosome"].astype(str).str.isdigit()]
        return cnvs_filtered[cnvs_filtered["chromosome"].astype(int).between(1, 22)]

    def get_filtered_cnvs_by_size(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on minimum size threshold."""
        return cnvs[cnvs["size"] >= self.min_size]

    def get_filtered_cnvs_by_log2_ratio(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on log2 ratio threshold."""
        return cnvs[cnvs["log2_corrected"].abs() >= self.log2_threshold]

    def get_filtered_cnvs(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on log2 ratio, size, and mask fraction."""
        if self.keep_autosomes_only:
            cnvs: DataFrame = self.get_filtered_cnvs_by_autosomes(cnvs)
        if self.min_size is not None:
            cnvs: DataFrame = self.get_filtered_cnvs_by_size(cnvs)
        if self.log2_threshold is not None:
            cnvs: DataFrame = self.get_filtered_cnvs_by_log2_ratio(cnvs)

        return cnvs

    def save_fga_summary(self, fga_summary: DataFrame, output_csv: str) -> None:
        """Save FGA summary to a CSV file."""
        fga_summary_path = Path(self.results_directory, output_csv)
        fga_summary.to_csv(fga_summary_path, index=False)

    def compute_fga(self, output_csv) -> None:
        """Process all CNV files in a directory and save FGA summary."""
        fga_summary = []

        for filepath in self.data_directory.glob("*.txt"):
            sample_id: str = self.get_sample_id_from_path(filepath)
            cnvs: DataFrame = self.read_cnv_file(filepath)

            # Calculate purity-corrected log2 ratios
            cnvs: DataFrame = self.get_purity_corrected_log2_ratios_cnvs(
                sample_id, cnvs
            )

            # Filter CNVs
            filtered_cnvs = self.get_filtered_cnvs(cnvs)

            # Compute FGA
            fga: float = self.calculate_fga_for_sample(filtered_cnvs)
            fga_summary.append({"SampleID": sample_id, "FGA": fga})

        fga_summary.sort(key=lambda x: int(x["SampleID"]))
        self.save_fga_summary(DataFrame(fga_summary), output_csv)
