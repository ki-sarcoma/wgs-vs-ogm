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
