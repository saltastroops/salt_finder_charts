from abc import ABC
from datetime import datetime, timedelta
from typing import Any, BinaryIO, Dict, List, NamedTuple, Optional, Tuple, Callable
import zipfile

from asteria import asteria
from asteria.instruments import BaseInstrument, Star
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from defusedxml.minidom import parseString
import pytz


"""
Metadata.

Any Metadata instance must be JSON-serializable.

"""

Metadata = Dict[str, Any]


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
        self.xml = mos_zip.read("Slitmask.xml").decode("UTF-8")
        mos_zip.close()
        doc = parseString(self.xml)

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


def estimated_position_angle(ra: Quantity, dec: Quantity, radius_range: Tuple[Quantity, Quantity]=(1 * u.arcmin, 3 * u.arcmin), mag_range: Tuple[float, float]=(15, 18), min_star_separation: Quantity=10 * u.arcsec) -> Optional[Quantity]:
    """
    Find a suitable position angle.

    The GAIA star catalog is used to find a suitable star with which a slit can be
    properly positioned, and the position angle of that star relative to the target is
    returned.

    Parameters
    ----------
    ra : Quantity
        Right ascension of the target, as an angle.
    dec : Quantity
        Declination of the target, as an angle.
    radius_range : pair of Quantity
        The inner and outer radius (as an angle) of the annulus in which a suitable may
        be located.
    mag_range : pair of float
        The minimum (brightest) and maximum (faintest) magnitude a suitable star may
        have.
    min_star_separation : Quantity
        The minimum angular distance a suitable star must have from neighbouring stars.

    """

    # set up the con ditions for the Gaia star catalog search
    instr = _build_position_angle_instrument(
        radius_range=radius_range,
        mag_range=mag_range,
        min_star_separation=min_star_separation
    )
    center = SkyCoord(ra, dec)
    instr.target = center

    # search for stars matching the conditions
    stars = asteria.view_sky(instr)
    matching_stars = asteria.find_best_stars(instr, stars)

    # don't calculate a position angle if there is no suitable star
    if len(matching_stars) == 0:
        return None

    # Find the best of the matching stars.
    best_star = sorted(matching_stars, key=_sorting_key(radius_range, mag_range))[0]
    best_star_coord = SkyCoord(best_star.ra, best_star.dec)

    return center.position_angle(best_star_coord)


def _build_position_angle_instrument(radius_range: Tuple[Quantity, Quantity], mag_range: Tuple[float, float], min_star_separation: Quantity) -> Callable[[Star], float]:
    """
    Create an "instrument" for the conditions stars must match for use in positioning a
    slit.

    Parameters
    ----------
    radius_range : pair of Quantity
        The inner and outer radius (as an angle) of the annulus in which a suitable may
        be located.
    mag_range : pair of float
        The minimum (brightest) and maximum (faintest) magnitude a suitable star may
        have.
    min_star_separation : Quantity
        The minimum angular distance a suitable star must have from neighbouring stars.

    Returns
    -------
    BaseInstrument
        The "instrument" with the search conditions.

    """

    class _PositionAngleInstrument(BaseInstrument):
        def best_stars(self, stars):
            return [s for s in stars if s.merit >= 4]

    instr = _PositionAngleInstrument(
        instr_name='PositionAngle',
        instr_fov=radius_range[1],
        inner_excl_distance=radius_range[0],
        nearby_limit=min_star_separation,
        bright_limit=mag_range[0],
        faint_limit=mag_range[1])

    return instr


def _sorting_key(radius_range: Tuple[Quantity, Quantity], mag_range: Tuple[float, float]):
    target_radius = radius_range[0]
    target_magnitude = 0.5 * sum(mag_range)

    return lambda star: abs((star.radius - target_radius).to_value(u.arcmin)) + 0.2 * abs(star.g_mag - target_magnitude)
