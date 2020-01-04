import enum
from io import BytesIO
from typing import BinaryIO

from salt_finder_charts.finder_charts import FinderChart, FinderChartMetadata


class OutputFormat(enum.Enum):
    """
    Available output formats.

    """

    PDF = "PDF"
    PNG = "PNG"
    SVG = "SVG"

    def extension(self):
        """
        The extension to use in filenames of files in this format.

        Returns
        -------
        str
            The file extension.

        """

        if self == OutputFormat.PDF:
            return "pdf"
        elif self == OutputFormat.PNG:
            return "png"
        elif self == OutputFormat.SVG:
            return "svg"
        else:
            raise ValueError(
                f"No file extension defined for output format {self.value}."
            )

    def mime_type(self) -> str:
        """
        The MIME type for this output format.

        Returns
        -------
        str
            MIME type

        """

        if self == OutputFormat.PDF:
            return "application/pdf"
        elif self == OutputFormat.PNG:
            return "image/png"
        elif self == OutputFormat.SVG:
            return "image/svg+xml"
        else:
            raise ValueError(f"No MIME type defined for output format {self.value}")


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


def output_svg(finder_chart: FinderChart, metadata: FinderChartMetadata) -> BinaryIO:
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
    finder_chart.plot.save(out, format="svg")
    return out
