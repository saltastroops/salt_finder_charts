from abc import ABC
import enum
import io
from typing import TextIO
import urllib.parse
import urllib.request

from astropy.coordinates.angles import Angle
import astropy.io.fits as pyfits
from astropy.units import Quantity
from astropy import units as u


class Survey(enum.Enum):
    """
    Image survey.

    """

    POSS1_BLUE = "POSS1 Blue"
    POSS1_RED = "POSS1 Red"
    POSS2UKSTU_BLUE = "POSS2/UKSTU Blue"
    POSS2UKSTU_IR = "POSS2/UKSTU IR"
    POSS2UKSTU_RED = "POSS2/UKSTU Red"
    TWO_MASS_H = "2MASS-H"
    TWO_MASS_J = "2MASS-J"
    TWO_MASS_K = "2MASS-K"


class ImageService(ABC):
    """
    A generator of FITS HDU list objects to use for the image in a finder chart.

    """

    def image(self, ra: Quantity, dec: Quantity) -> pyfits.HDUList:
        """
        Generate a HDU list to use for the image in a finder chart.

        Parameters
        ----------
        ra : Quantity
            Right ascension of the image centre.
        dec : Quantity
            Declination of the image centre.

        Returns
        -------
        HDUList
            The HDU list object for the image.

        """

        raise NotImplementedError

    def source(self) -> str:
        """
        The name of the image source (such as a survey).

        Returns
        -------
        str
            The image source.

        """

        raise NotImplementedError


class SurveyImageService(ImageService):
    """
    Image service for getting an image using one of a set of image surveys.

    Parameters
    ----------
    survey : Survey
        Image survey.

    """

    STSCI_SURVEYS = [
        Survey.POSS2UKSTU_RED,
        Survey.POSS2UKSTU_BLUE,
        Survey.POSS2UKSTU_IR,
        Survey.POSS1_RED,
        Survey.POSS1_BLUE,
    ]

    SKY_VIEW_SURVEYS = [Survey.TWO_MASS_J, Survey.TWO_MASS_H, Survey.TWO_MASS_K]

    def __init__(self, survey: Survey):
        self.survey = survey

    def image(self, ra: Quantity, dec: Quantity) -> pyfits.HDUList:
        # grab 10' x 10' image from server and pull it into pyfits
        if self.survey in SurveyImageService.STSCI_SURVEYS:
            survey_identifiers = {
                Survey.POSS2UKSTU_RED: "poss2ukstu_red",
                Survey.POSS2UKSTU_BLUE: "poss2ukstu_blue",
                Survey.POSS2UKSTU_IR: "poss2ukstu_ir",
                Survey.POSS1_RED: "poss1_red",
                Survey.POSS1_BLUE: "poss1_blue",
            }
            url = "http://archive.stsci.edu/cgi-bin/dss_search"
            params = urllib.parse.urlencode(
                {
                    "v": survey_identifiers[self.survey],
                    "r": "%f" % ra.to_value(u.deg),
                    "d": "%f" % dec.to_value(u.deg),
                    "e": "J2000",
                    "h": 10.0,
                    "w": 10.0,
                    "f": "fits",
                    "c": "none",
                }
            ).encode("utf-8")
        elif self.survey in SurveyImageService.SKY_VIEW_SURVEYS:
            survey_identifiers = {
                Survey.TWO_MASS_J: "2mass-j",
                Survey.TWO_MASS_H: "2mass-h",
                Survey.TWO_MASS_K: "2mass-k",
            }
            ra = Angle(ra)
            dec = Angle(dec)
            url = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
            params = urllib.parse.urlencode(
                {
                    "Position": "'%d %d %f, %d %d %f'"
                    % (
                        round(ra.hms[0]),
                        ra.hms[1],
                        ra.hms[2],
                        round(dec.dms[0]),
                        abs(dec.dms[1]),
                        abs(dec.dms[2]),
                    ),
                    "Survey": survey_identifiers[self.survey],
                    "Coordinates": "J2000",
                    "Return": "FITS",
                    "Pixels": 700,
                    "Size": 0.1667,
                }
            ).encode("utf-8")
        else:
            raise Exception(f"Unsupported survey: {self.survey}")
        fits_data = io.BytesIO()
        data = urllib.request.urlopen(url, params).read()
        fits_data.write(data)
        fits_data.seek(0)
        return pyfits.open(fits_data)

    def source(self) -> str:
        return str(self.survey.value)


class FITSImageService(ImageService):
    """
    Image service for an image from a FITS file.

    This class should be used with care, as it does not check whether the FITS file is
    consistent with the right ascension abd declination passed to the image method. So
    you may easily end up with completely nonsensical finder charts.

    Parameters
    ----------
    fits : text stream
         Text stream containing the FITS file.

    """

    def __init__(self, fits: TextIO):
        self.hdu = pyfits.open(fits)

    def image(self, ra: Quantity, dec: Quantity) -> pyfits.HDUList:
        """
        Return the HDU list object from the FITS file.

        Note that this method does not use the passed right ascension and declination
        and that it does not check whether they are consistent with the image.

        It is thus the user's responsibility to ensure that the FITS file and the given
        right ascension and declination are consistent.

        Parameters
        ----------
        ra : Quantity
            Right ascension.
        dec : Quantity
            Declination.

        Returns
        -------
        HDUList
            The HDU list object for the image.

        """

        return self.hdu

    def source(self) -> str:
        return "User-supplied FITS"
