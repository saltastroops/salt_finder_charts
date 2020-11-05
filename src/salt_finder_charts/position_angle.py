from typing import Optional, Tuple

from gilmenel import gilmenel
from gilmenel.instruments import BaseInstrument, Star
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

MIN_RADIUS = 1 * u.arcmin
MAX_RADIUS = 3 * u.arcmin

MIN_MAG = 15
MAX_MAG = 18

MIN_STAR_SEPARATION = 10 * u.arcsec


def estimated_position_angle(
    ra: Quantity,
    dec: Quantity,
    radius_range: Tuple[Quantity, Quantity] = (MIN_RADIUS, MAX_RADIUS),
    mag_range: Tuple[float, float] = (MIN_MAG, MAX_MAG),
    min_star_separation: Quantity = MIN_STAR_SEPARATION,
) -> Optional[Quantity]:
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
        min_star_separation=min_star_separation,
    )
    center = SkyCoord(ra, dec)
    instr.point_to(center)

    # search for stars matching the conditions
    gilmenel.init()
    stars = gilmenel.view_sky(instr)
    matching_stars = gilmenel.find_best_stars(instr, stars)

    # don't calculate a position angle if there is no suitable star
    if len(matching_stars) == 0:
        return None

    # Find the best of the matching stars.
    best_star = sorted(matching_stars, key=_sorting_key(radius_range, mag_range))[0]
    best_star_coord = SkyCoord(best_star.ra, best_star.dec)

    return center.position_angle(best_star_coord)


def _build_position_angle_instrument(
    radius_range: Tuple[Quantity, Quantity],
    mag_range: Tuple[float, float],
    min_star_separation: Quantity,
) -> BaseInstrument:
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
        name="PositionAngle",
        instr_fov=radius_range[1],
        inner_excl_distance=radius_range[0],
        nearby_limit=min_star_separation,
        bright_limit=mag_range[0],
        faint_limit=mag_range[1],
    )

    return instr


def _sorting_key(
    radius_range: Tuple[Quantity, Quantity], mag_range: Tuple[float, float]
):
    """
    A key function for sorting stars.

    The most suitable star comes first in a sorted list.

    Parameters
    ----------
    radius_range : pair of Quantity
        The inner and outer radius (as an angle) of the annulus in which a suitable may
        be located.
    mag_range : pair of float
        The minimum (brightest) and maximum (faintest) magnitude a suitable star may
        have.

    Returns
    -------
    float
        A function which can be used for sorting.

    """

    target_radius = radius_range[0]
    target_magnitude = 0.5 * sum(mag_range)

    return lambda star: abs(
        (star.radius - target_radius).to_value(u.arcmin)
    ) + 0.2 * abs(star.g_mag - target_magnitude)
