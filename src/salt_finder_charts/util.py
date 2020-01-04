from datetime import datetime, timedelta
from typing import List, NamedTuple, BinaryIO
import zipfile

import astropy.units as u
from astropy.units import Quantity
from defusedxml.minidom import parseString
import pytz


class MagnitudeRange(NamedTuple):
    """
    A magnitude range.

    Properties
    ----------
    bandpass : str
        Bandpass (such as V) for which the magnitudes are given.
    max_magnitude : float
        Maximum (faintest) magnitude.
    min_magnitude : float
        Minimum (brightest) magnitude.

    """

    bandpass: str
    max_magnitude: float
    min_magnitude: float


class MOSMaskSlit(NamedTuple):
    """
    A MOS mask slit.

    Attributes
    ----------
    ra : Quantity
        Right ascension of the slit center.
    dec : Quantity
        Declination of the slit center.
    tilt : Quantity
        Tilt of the slit.
    width: Quantity
        Slit width, as an angle.
    height : Quantity
        Slit height, as an angle.

    """

    ra: Quantity
    dec: Quantity
    tilt: Quantity
    width: Quantity
    height: Quantity


class MOSMaskReferenceStar(NamedTuple):
    """
    A MOS mask reference star.

    Attributes
    ----------
    ra : Quantity
        Right ascension of the star.
    dec : Quantity
        Declination of the star.

    """

    ra: Quantity
    dec: Quantity


class MOSMask:
    """
    Access to properties of a MOS mask definition.

    Parameters
    ----------
    rsmt : BinaryIO
        RSMT file for the MOS mask.

    Attributes
    ----------
    dec : Quantity
        Declination of the mask center.
    reference_stars : list of MOSMaskReferenceStar
        List of reference stars defined for the mask.
    right_ascension : Quantity
        Right ascensio of the mask center.
    rotation_angle : Quantity
        Rotation angle of the mask.
    slits : list of MOSMaskSlit
        Right ascension of the slits defined for the mask.
    """

    def __init__(self, rsmt: BinaryIO):
        # extract the MOS mask definition from the RSMT file
        mos_zip = zipfile.ZipFile(rsmt)
        mos_xml = mos_zip.read("Slitmask.xml").decode("UTF-8")
        mos_zip.close()
        doc = parseString(mos_xml)

        # extract the mask position and rotation angle
        parameter_elenents = doc.getElementsByTagName("parameter")
        parameters = {}
        for par in parameter_elenents:
            name = par.getAttribute("name")
            val = par.getAttribute("value")
            parameters[name] = val
        self._ra = float(parameters["CENTERRA"]) * u.deg
        self._dec = float(parameters["CENTERDEC"]) * u.deg
        self._pa = float(parameters["ROTANGLE"]) * u.deg

        # extract the slit details
        slit_elements = doc.getElementsByTagName("slit")
        self._slits = []
        for slit in slit_elements:
            tilt = 0.0 * u.deg
            if "tilt" in slit.attributes.keys():
                tilt = float(slit.attributes["tilt"].value) * u.deg
            self._slits.append(
                MOSMaskSlit(
                    ra=float(slit.attributes["xce"].value) * u.deg,
                    dec=float(slit.attributes["yce"].value) * u.deg,
                    width=float(slit.attributes["width"].value) * u.arcsec,
                    height=float(slit.attributes["length"].value) * u.arcsec,
                    tilt=tilt,
                )
            )

        # extract the reference star details
        reference_star_elements = doc.getElementsByTagName("refstar")
        self._reference_stars = []
        for ref in reference_star_elements:
            self._reference_stars.append(
                MOSMaskReferenceStar(
                    ra=float(ref.attributes["xce"].value) * u.deg,
                    dec=float(ref.attributes["yce"].value) * u.deg,
                )
            )

    @property
    def declination(self) -> Quantity:
        return self._dec

    @property
    def reference_stars(self) -> List[MOSMaskReferenceStar]:
        return self._reference_stars

    @property
    def right_ascension(self) -> Quantity:
        return self._ra

    @property
    def position_angle(self) -> Quantity:
        return self._pa

    @property
    def slits(self) -> List[MOSMaskSlit]:
        return self._slits


def julian_day_start(t: datetime) -> datetime:
    """
    Return the start of a Julian day containing a datetime.

    If the given datetime is noon (UTC) or later this is noon (UTC) of the given date,
    otherwise it is noon (UTC) of the previous day.

    Parameters
    ----------
    t : datetime
        Datetime (must be timezone-aware).

    Returns
    -------
    datetime
        Default start time.

    """

    if not t.tzinfo:
        raise ValueError("The datetime must be timezone-aware.")

    t = t.astimezone(pytz.utc)
    noon = t.replace(hour=12, minute=0, second=0, microsecond=0)
    if t < noon:
        return noon - timedelta(days=1)
    else:
        return noon


def julian_day_end(t: datetime) -> datetime:
    """
    Return the end of a Julian day containing a datetime.

    If the given datetime is noon (UTC) or earlier this is noon (UTC) of the given date,
    otherwise it is noon (UTC) of the following day.

    Parameters
    ----------
    t : datetime
        Datetime (must be timezone-aware).

    Returns
    -------
    datetime
        Default end time.

    """

    if not t.tzinfo:
        raise ValueError("The datetime must be timezone-aware.")

    t = t.astimezone(pytz.utc)
    noon = t.replace(hour=12, minute=0, second=0, microsecond=0)
    if t < noon:
        return noon
    else:
        return noon + timedelta(days=1)
