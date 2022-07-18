from pathlib import Path
import shutil
import os

WORK_DIRECTORY = Path(os.path.dirname(__file__))

DATA_DIRECTORY = WORK_DIRECTORY.parents[2].joinpath("data")
os.makedirs(DATA_DIRECTORY, exist_ok=True)

CHROMOSOME_13_REFERENCE = DATA_DIRECTORY.joinpath("chromosome13.fasta")
CHROMOSOME_13_READS = DATA_DIRECTORY.joinpath("reads_13.fastq.gz")

CHROMOSOME_19_REFERENCE = DATA_DIRECTORY.joinpath("chromosome19.fasta")
CHROMOSOME_19_READS = DATA_DIRECTORY.joinpath("reads_19.fastq.gz")
