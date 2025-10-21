from pathlib import Path
from numpy import log2, maximum
from pandas import DataFrame, Series, read_csv

from cnv_analyzer import CNVAnalyzer


class WGSCNVAnalyzer(CNVAnalyzer):
    def __init__(
        self,
        data_directory: Path = None,
        results_directory: Path = None,
        genome_size: int = 2_881_033_286,  # Autosomal hg19
        keep_autosomes_only: bool = True,
        min_size: int = None,
        log2_threshold: float = None,
    ):
        super().__init__(
            data_directory=data_directory,
            results_directory=results_directory,
            genome_size=genome_size,
            keep_autosomes_only=keep_autosomes_only,
            min_size=min_size,
            log2_threshold=log2_threshold,
        )

    def get_data_directory(self) -> Path:
        """Get the data directory path."""
        script_dir = Path(__file__).resolve().parent
        return Path(script_dir.parent, "data", "WGS").resolve()

    def get_sample_id_from_path(self, filepath: Path) -> str:
        """Extract sample ID from file path."""
        return filepath.stem.split("-")[2][1:-2]

    def get_size_cnvs(self, cnvs: DataFrame) -> DataFrame:
        """Calculate size column for CNVs."""
        cnvs["size"] = cnvs["end"] - cnvs["start"] + 1
        return cnvs

    def calculate_fga_for_sample(self, cnvs: DataFrame) -> float:
        """Calculate fraction of genome altered."""
        altered_bp: int = cnvs["size"].sum(skipna=True)
        fga: float = altered_bp / self.genome_size
        return round(fga, 4)

    def compute_fga(self, output_csv) -> None:
        """Process all CNV files in a directory and save FGA summary."""
        fga_summary = []

        for filepath in self.data_directory.glob("*.cns"):
            sample_id: str = self.get_sample_id_from_path(filepath)
            cnvs: DataFrame = self.read_cnv_file(filepath)

            # Calculate size column
            cnvs: DataFrame = self.get_size_cnvs(cnvs)

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


if __name__ == "__main__":
    ogm_cnv_analyzer = WGSCNVAnalyzer(min_size=50_000, log2_threshold=0.2)
    ogm_cnv_analyzer.compute_fga(output_csv="wgs_fga.csv")
