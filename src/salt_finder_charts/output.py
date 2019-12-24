from abc import ABC
import enum
from io import BytesIO
from typing import BinaryIO, Optional

from salt_finder_charts.finder_charts import FinderChart, FinderChartMetadata


class OutputFormat(enum.Enum):
    """
    Available output formats.

    """

    PDF = "PDF"
    PNG = "PNG"


def output_pdf(finder_chart: FinderChart, metadata: FinderChartMetadata) -> BinaryIO:
    """
    Generate a binary stream with a PDF document containing the given finder chart.

    Parameters
    ----------
    finder_chart : FinderChart
        Finding chart.
    metadata : FinderChartMetadata
        Finding chart metadata.

    Returns
    -------
    BytesIO
        A binary stream containing a PDF document with the finder chart.

    """

    out = BytesIO()
    finder_chart.plot.save(out, format="pdf")
    return out


def output_png(finder_chart: FinderChart, metadata: FinderChartMetadata) -> BinaryIO:
    """
    Generate a binary stream with a PDF document containing the given finder chart.

    Parameters
    ----------
    finder_chart : FinderChart
        Finding chart.
    metadata : FinderChartMetadata
        Finding chart metadata.

    Returns
    -------
    BytesIO
        A binary stream containing a PDF document with the finder chart.

    """

    out = BytesIO()
    finder_chart.plot.save(out, format="png")
    return out
