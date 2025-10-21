from pathlib import Path
from pandas import DataFrame, Series, read_csv
from numpy import log2, maximum

from cnv_analyzer import CNVAnalyzer


class OGMCNVAnalyzer(CNVAnalyzer):
    def __init__(
        self,
        data_directory: Path = None,
        results_directory: Path = None,
        genome_size: int = 2_875_001_522,  # Autosomal hg38
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
        return Path(script_dir.parent, "data", "OGM").resolve()

    def get_sample_id_from_path(self, filepath: Path) -> str:
        """Extract sample ID from file path."""
        return filepath.stem.split("_")[0]

    def read_cnv_file(self, filepath) -> DataFrame:
        """Read a single CNV file and return a dataframe."""
        return read_csv(filepath, sep="\t", header=5)

    def get_purity_corrected_log2_ratios_cnvs(
        self, sample_id: str, cnvs: DataFrame
    ) -> DataFrame:
        """Calculate purity-corrected log2 ratios for CNVs."""
        sample_purity: float = self.get_sample_purity(sample_id)
        copy_number_obs = cnvs["fractionalCopyNumber"]
        copy_number_tumor = maximum(
            (copy_number_obs - 2 * (1 - sample_purity)) / sample_purity, 0.01
        )
        cnvs["log2_corrected"] = round(log2(copy_number_tumor / 2), 2)
        return cnvs

    def get_filtered_cnvs_by_autosomes(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs to keep only autosomes (chromosomes 1-22)."""
        cnvs_filtered = cnvs[cnvs["Chromosome"].astype(str).str.isdigit()]
        return cnvs_filtered[cnvs_filtered["Chromosome"].astype(int).between(1, 22)]

    def get_filtered_cnvs_by_size(self, cnvs: DataFrame) -> DataFrame:
        """Filter CNVs based on minimum size threshold."""
        return cnvs[cnvs["Size"] >= self.min_size]

    def calculate_fga_for_sample(self, cnvs: DataFrame) -> float:
        """Calculate fraction of genome altered."""
        altered_bp: int = cnvs["Size"].sum(skipna=True)
        fga: float = altered_bp / self.genome_size
        return round(fga, 4)


if __name__ == "__main__":
    ogm_cnv_analyzer = OGMCNVAnalyzer(min_size=50_000, log2_threshold=0.2)
    ogm_cnv_analyzer.compute_fga(output_csv="ogm_fga.csv")
