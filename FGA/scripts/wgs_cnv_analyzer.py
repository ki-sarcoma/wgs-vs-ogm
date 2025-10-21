from pathlib import Path
from pandas import DataFrame, Series

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
        chr_col_name: str = "chromosome",
        size_col_name: str = "size",
    ):
        super().__init__(
            data_directory=data_directory,
            results_directory=results_directory,
            genome_size=genome_size,
            keep_autosomes_only=keep_autosomes_only,
            min_size=min_size,
            log2_threshold=log2_threshold,
            chr_col_name=chr_col_name,
            size_col_name=size_col_name,
        )

    def get_data_directory(self) -> Path:
        """Get the data directory path."""
        script_dir = Path(__file__).resolve().parent
        return Path(script_dir.parent, "data", "WGS").resolve()

    def get_sample_id_from_path(self, filepath: Path) -> str:
        """Extract sample ID from file path."""
        return filepath.stem.split("-")[2][1:-2]

    def get_copy_number_obs(self, cnvs: DataFrame) -> Series:
        """Get observed copy number from CNVs."""
        return 2 * (2 ** cnvs["log2"])

    def get_size_cnvs(self, cnvs: DataFrame) -> DataFrame:
        """Calculate size column for CNVs."""
        cnvs["size"] = cnvs["end"] - cnvs["start"] + 1
        return cnvs


if __name__ == "__main__":
    ogm_cnv_analyzer = WGSCNVAnalyzer(min_size=50_000, log2_threshold=0.2)
    ogm_cnv_analyzer.compute_fga(output_csv="wgs_fga.csv")
