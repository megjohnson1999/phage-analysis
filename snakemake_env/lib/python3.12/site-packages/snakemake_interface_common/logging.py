__author__ = "Johannes Köster"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


def get_logger():
    """Retrieve the logger singleton from snakemake."""
    from snakemake.logging import logger

    return logger
