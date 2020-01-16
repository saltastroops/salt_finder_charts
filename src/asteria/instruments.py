import abc
import re

from typing import (
    List,
    Optional,
)

from astropy import units as u
from astropy.coordinates import SkyCoord


class Star(object):
    def __init__(
            self,
            ra: u.deg,
            dec: u.deg,
            radius: u.arcmin,
            g_mag: float,
            merit: Optional[int] = None,
    ):
        self.sky_coord = SkyCoord(ra, dec, frame='icrs')
        self.ra = self.sky_coord.ra
        self.dec = self.sky_coord.dec

        self.radius = radius  # arcmin
        self.g_mag = g_mag

        # the stars coordinates in the instrument's reference frame
        self.instr_coord = None

        # at first, all stars have a merit of 0
        # as the pass each filter, the merit value increases
        # this will allow for operations in marginal conditions
        if merit:
            self.merit = merit
        else:
            self.merit = 0

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Star):
            return (u.isclose(self.ra, other.ra, rtol=1e-6)
                    and u.isclose(self.dec, other.dec, rtol=1e-6)
                    and u.isclose(self.radius, other.radius, rtol=1e-2)
                    and self.g_mag == other.g_mag
                    and self.merit == other.merit)
        return False  # pragma: no cover

    def __str__(self):
        return (f"ra={self.ra:.6f}, "
                f"dec={self.dec:.6f}, "
                f"radius={self.radius:.3f}, "
                f"g_mag={self.g_mag:.2f}, "
                f"merit={self.merit}")  # pragma: no cover

    def __repr__(self):
        return re.sub(
            r'0+ ', '0 ', f"\nStar({self.ra.to_value(u.deg):.6f} * d,"
            f" {self.dec.to_value(u.deg):.6f} * d,"
            f" {self.radius.to_value(u.arcmin):.3f} * m,"
            f" {self.g_mag}, merit={self.merit})")  # pragma: no cover


class BaseInstrument(metaclass=abc.ABCMeta):
    def __init__(
            self,
            instr_name: str,
            instr_fov: u.arcmin,
            inner_excl_distance: u.arcmin,
            nearby_limit: u.arcsec,
            bright_limit: float,
            faint_limit: float,
    ):
        self.instr_name = instr_name
        self.instr_fov = instr_fov  # arcminutes radius
        self.inner_excl_distance = inner_excl_distance  # arcminutes radius
        self.nearby_limit = nearby_limit  # arcseconds diameter
        self.bright_limit = bright_limit
        self.faint_limit = faint_limit
        self.inner_excl_shape = 'circle'

        self.target = None  # after init > instr.target = SkyCoord(ra=x, dec=y)

    def star_available(self, star: Star) -> bool:
        '''Check if star location falls within allowable geometry'''

        # Check if star falls within FOV
        if star.radius > self.instr_fov:
            return False

        # Check if star falls within exclusion zone
        if self.inner_excl_shape == 'circle':
            if star.radius <= self.inner_excl_distance:
                return False

        # TODO: fix the errors in this case, the result does not seem right
        elif self.inner_excl_shape == 'square':
            delta_ra = abs(star.sky_coord.ra - self.target.ra)
            delta_dec = abs(star.sky_coord.dec - self.target.dec)

            if (delta_ra <= self.inner_excl_distance
                    and delta_dec <= self.inner_excl_distance):
                return False
        else:
            raise NotImplementedError

        return True

    def filter_geometry(self, stars: List[Star]) -> List[Star]:
        for s in stars:
            if self.star_available(s):
                s.merit = 1

        return stars

    def filter_nearby_pairs(self, stars: List[Star]) -> List[Star]:
        for s in stars:

            def s_is_close_to(t):
                return (t.g_mag < s.g_mag + 2
                        and  # t mag brighter than s mag -2
                        abs(t.ra - s.ra) <= self.nearby_limit and
                        abs(t.dec - s.dec) <= self.nearby_limit)

            nearby_stars = any(  # any returns True at first test star found
                True for t in stars
                if s_is_close_to(t) and s is not t)  # t -> test star

            if nearby_stars is False:
                s.merit = 2

        return stars

    def filter_magnitudes(self, stars: List[Star]) -> List[Star]:
        for s in stars:
            # remove faint stars
            if s.g_mag > self.faint_limit:
                pass  # leave faint stars at their current merit value

            # remove bright stars
            elif s.g_mag < self.bright_limit:
                s.merit = 3

            # just right
            else:
                s.merit = 4

        return stars

    def filter(self, stars: List[Star]) -> List[Star]:
        # TODO: improve algorithm to check for stars in layers of magnitude
        self.filter_geometry(stars)
        self.filter_nearby_pairs([s for s in stars if s.merit == 1])
        self.filter_magnitudes([s for s in stars if s.merit == 2])

        return stars

    @abc.abstractmethod
    def best_stars(self, stars: List[Star]) -> List[Star]:
        '''Return the best guide stars from the given selection'''


class GapInstrument(BaseInstrument):
    def __init__(
            self,
            instr_name: str,
            instr_fov: u.arcmin,
            inner_excl_distance: u.arcmin,
            nearby_limit: u.arcsec,
            bright_limit: float,
            faint_limit: float,
            slit_gap_radius: u.arcmin,
            slit_gap_angle: u.deg,
    ):
        super().__init__(
            instr_name,
            instr_fov,
            inner_excl_distance,
            nearby_limit,
            bright_limit,
            faint_limit,
        )
        self.slit_gap_radius = slit_gap_radius  # arcmin
        self.slit_gap_angle = slit_gap_angle  # degrees

    def star_available(self, star: Star) -> bool:
        '''Check if star location falls within allowable geometry'''

        result = super().star_available(star)
        if result is False:
            return result

        # Check if star falls within a vertical gap of distance slit_gap_radius
        star.instr_coord = star.sky_coord.copy()

        star.instr_coord = star.instr_coord.transform_to(self.instr_frame)
        ra = (star.instr_coord.lon.to_value(u.deg) + 180) % 360 - 180

        if abs(ra) <= self.slit_gap_radius.to_value(u.deg):
            return False

        return True

    def filter_geometry(self, stars: List[Star]) -> List[Star]:
        self.instr_frame = self.target.skyoffset_frame(
            rotation=self.slit_gap_angle)

        for s in stars:
            if self.star_available(s):
                s.merit = 1

        return stars

    def filter(self, stars: List[Star]) -> List[Star]:
        self.filter_geometry(stars)
        super().filter_nearby_pairs([s for s in stars if s.merit == 1])
        super().filter_magnitudes([s for s in stars if s.merit == 2])

        return stars

    def best_stars(self, stars: List[Star]) -> List[Star]:
        '''Return the best guide stars from the given selection'''
        raise NotImplementedError
