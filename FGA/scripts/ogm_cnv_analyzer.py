from pathlib import Path
from pandas import DataFrame, Series, read_csv

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
        chr_col_name: str = "Chromosome",
        size_col_name: str = "Size",
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
        return Path(script_dir.parent, "data", "OGM").resolve()

    def get_sample_id_from_path(self, filepath: Path) -> str:
        """Extract sample ID from file path."""
        return filepath.stem.split("_")[0]

    def read_cnv_file(self, filepath) -> DataFrame:
        """Read a single CNV file and return a dataframe."""
        return read_csv(filepath, sep="\t", header=5)

    def get_copy_number_obs(self, cnvs: DataFrame) -> Series:
        """Get observed copy number from CNVs."""
        return cnvs["fractionalCopyNumber"]


if __name__ == "__main__":
    ogm_cnv_analyzer = OGMCNVAnalyzer(min_size=50_000, log2_threshold=0.2)
    ogm_cnv_analyzer.compute_fga(output_csv="ogm_fga.csv")
