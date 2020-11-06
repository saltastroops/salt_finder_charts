from abc import ABC
import bisect
from datetime import datetime, timedelta
from typing import List, NamedTuple, Optional, Tuple

from astropy.units import Quantity
from astropy import units as u
from astroquery.jplhorizons import Horizons
from dateutil.parser import parse
from dateutil.tz import tzutc

from salt_finder_charts.util import MagnitudeRange, Metadata


class Ephemeris(NamedTuple):
    """
    An ephemeris.

    Properties
    ----------
    dec : Quantity
        Declination, as an angle between -90 and 90 degrees.
    dec_rate : Quantity
        Proper motion (rate of change) of the declination.
    epoch : datetime
        The epoch, i.e. a datetime for which this ephemeris is valid.
    magnitude_range : Optional[MagnitudeRange]
        The range of magnitudes the target may have, plus the bandpass for which the
         magnitudes are given.
    ra : Quantity
        Right ascension, as an angle between 0 and 360 degrees.
    ra_rate : Quantity
        Proper motion (rate of change) of the right ascension.
    """

    dec: Quantity
    dec_rate: Quantity
    epoch: datetime
    magnitude_range: Optional[MagnitudeRange]
    ra: Quantity
    ra_rate: Quantity


class EphemerisService(ABC):
    """
    An abstract base class for services providing ephemerides.

    """

    def ephemerides(self, start_time: datetime, end_time: datetime) -> List[Ephemeris]:
        """
        A list of ephemeris values that covers at least a given time interval.

        The ephemeris values are sorted by time (from earliest to latest).

        The first returned ephemeris may be for a time earlier than the start time and
        the last returned ephemeris may be for a time later than the end time.

        Parameters
        ----------
        start_time : datetime
            Start time.
        end_time : datetime
            End time.

        Returns
        -------
        list of Ephemeris
            List of ephemeris values at least covering the interval.

        """

        raise NotImplementedError

    def is_sidereal_target(self) -> bool:
        """
        Whether the ephemerides are for a sidereal target.

        Returns
        -------
        bool
            Whether the ephemerides are for a sidereal target.

        """

        raise NotImplementedError

    @staticmethod
    def center_position(ephemerides: List[Ephemeris]) -> Tuple[Quantity, Quantity]:
        """
        The center of the positions defined in a list of ephemerides.

        The center position is returned as a tuple of the right ascension and the
        declination.

        If the maximum and minimum right ascension differ by more than 180 degrees, it
        is assumed that the target right ascension crosses 360 degrees mark.

        If the center right ascension is greater than 360 degrees, the equivalent angle
        netween 0 and 360 degrees is used instead.

        Parameters
        ----------
        ephemerides : list of Ephemeris
            Ephemerides.

        Returns
        -------
        Tuple[Quantity, Quantity]
            The center position.

        """

        ra_key = lambda e: e.ra
        ra_min = min(ephemerides, key=ra_key).ra
        ra_max = max(ephemerides, key=ra_key).ra

        if ra_max - ra_min > 180 * u.deg:
            ra_min += 360 * u.deg

        dec_key = lambda e: e.dec
        dec_min = min(ephemerides, key=dec_key).dec
        dec_max = max(ephemerides, key=dec_key).dec

        ra_center = (ra_min + ra_max) / 2.0
        dec_center = (dec_min + dec_max) / 2.0

        if ra_center > 360 * u.deg:
            ra_center -= 360 * u.deg

        return ra_center, dec_center

    @staticmethod
    def find_magnitude_range(ephemerides: List[Ephemeris]) -> Optional[MagnitudeRange]:
        """
        The overall magnitude range for a list of ephemerides.

        If the magnitude range of any ephemeris is None, or if the minimum or maximum
        magnitude, ot if the bandpass is not the same for all ephemerides, or if the
        given list of ephemerides is empty, None is returned.

        Parameters
        ----------
        ephemerides : list of Ephemeris
            Ephemerides.

        Returns
        -------
        Optional[MagnitudeRange]
            Magnitude range.

        """

        if len(ephemerides) == 0:
            return None

        bandpass = None
        if ephemerides[0].magnitude_range:
            bandpass = ephemerides[0].magnitude_range.bandpass
        if bandpass is None:
            return None

        min_magnitude = 1e100
        max_magnitude = 1e-100
        for ephemeris in ephemerides:
            mr = ephemeris.magnitude_range
            if (
                mr is None
                or mr.min_magnitude is None
                or mr.max_magnitude is None
                or mr.bandpass != bandpass
            ):
                return None
            if mr.min_magnitude < min_magnitude:
                min_magnitude = mr.min_magnitude
            if mr.max_magnitude > max_magnitude:
                max_magnitude = mr.max_magnitude

        return MagnitudeRange(
            min_magnitude=min_magnitude, max_magnitude=max_magnitude, bandpass=bandpass
        )

    def metadata(self) -> Metadata:
        """
        Metadata characterising this ephemeris service.

        Returns
        -------
        Metadata
            Metadata for the ephemeris service.

        """

        raise NotImplementedError


def _cover_time_interval(
    ephemerides: List[Ephemeris], start_time: datetime, end_time: datetime
) -> List[Ephemeris]:
    """
    The smallest possible list of ephemerides from a given list that completely covers
    a given time interval.

    Parameters
    ----------
    ephemerides : list of Ephemeris
    start_time : timezone-aware start time
    end_time : timezone-aware end time

    Returns
    -------
    list of Ephemeris
        The list of ephemerides covering the time interval.

    """

    if start_time.tzinfo is None or end_time.tzinfo is None:
        raise ValueError("The start and end time must be timezone-aware.")

    if start_time >= end_time:
        raise ValueError("The start time must be earlier than the end time")

    all_times = [e.epoch for e in ephemerides]
    start_index = bisect.bisect_right(all_times, start_time)
    if not start_index and all_times[0] != start_index:
        raise ValueError("The start time isn't covered by the ephemerides.")
    end_index = bisect.bisect(all_times, end_time)
    if end_index == len(all_times) and all_times[-1] != end_time:
        raise ValueError("The end time isn't covered by the ephemerides.")
    if end_index > 0 and all_times[end_index - 1] == end_time:
        end_index -= 1

    return ephemerides[start_index - 1 : end_index + 1]


class ConstantEphemerisService(EphemerisService):
    """
    An ephemeris generator for constant ephemerides.

    Parameters
    ----------
    ra : Quantity
        Right ascension, as an angle.
    dec : Quantity
        Declination, as an angle.
    magnitude_range : MagnitudeRange
        Magnitude range (optional).

    """

    def __init__(
        self, ra: Quantity, dec: Quantity, magnitude_range: Optional[MagnitudeRange]
    ) -> None:
        self.ra = ra
        self.dec = dec
        self.magnitude_range = magnitude_range

    def ephemerides(self, start_time: datetime, end_time: datetime) -> List[Ephemeris]:
        # enforce timezones
        if start_time.tzinfo is None or end_time.tzinfo is None:
            raise ValueError("The start and end time must be timezone-aware.")

        return [
            Ephemeris(
                dec=self.dec,
                dec_rate=0 * u.deg / u.second,
                epoch=start_time,
                magnitude_range=self.magnitude_range,
                ra=self.ra,
                ra_rate=0 * u.deg / u.second,
            ),
            Ephemeris(
                dec=self.dec,
                dec_rate=0 * u.deg / u.second,
                epoch=end_time,
                magnitude_range=self.magnitude_range,
                ra=self.ra,
                ra_rate=0 * u.deg / u.second,
            ),
        ]

    def is_sidereal_target(self) -> bool:
        return True

    def metadata(self) -> Metadata:
        return dict()


class HorizonsEphemerisService(EphemerisService):
    """
    An ephemeris generator using the JPL Horizons service.

    In order to avoid missing ephemerides, you should choose a start time at least two
    days earlier than the start time from which you generate finder charts, and you
    should choose an end time at least two days later than the end time until which
    you generate finder charts. This is necessary as the calculated visibility windows
    are not strictly confined to the time interval for which finder charts shall be
    generated.

    Parameters
    ----------
    object_id : str
        Identifier of the object whose ephemerides are generated.
    start_time : datetime
        Time of the first ephemeris to get from Horizons (must be timezone-aware).
    end_time : datetime
        Time of the last ephemeris to get from Horizons (must be timezone-aware).
    stepsize : Quantity
        Time between ephemerides queried from Horizons (must be at least 5 minutes).

    """

    def __init__(
        self,
        object_id: str,
        start_time: datetime,
        end_time: datetime,
        stepsize: Quantity,
    ):
        SALT_OBSERVATORY_ID = "B31"

        # enforce timezones
        if start_time.tzinfo is None or end_time.tzinfo is None:
            raise ValueError("The start and end time must be timezone-aware.")

        # avoid overly excessive queries
        self.stepsize = stepsize
        if self.stepsize < 5 * u.minute:
            raise ValueError("The sampling interval must be at least 5 minutes.")

        # query Horizons
        self.object_id = object_id
        start = start_time.astimezone(tzutc()).strftime("%Y-%m-%d %H:%M:%S")
        # Make sure the whole time interval is covered by the queried ephemerides
        end_time_with_margin = end_time + timedelta(seconds=stepsize.to_value(u.second))
        stop = end_time_with_margin.astimezone(tzutc()).strftime("%Y-%m-%d %H:%M:%S")
        # Horizons requires an int for the step size. As round() might call NumPy's
        # round method and thus produce a float, we have to round "manually" using
        # the int function.
        step = f"{int(0.5 + stepsize.to_value(u.minute))}m"
        obj = Horizons(
            id=self.object_id,
            location=SALT_OBSERVATORY_ID,
            epochs={"start": start, "stop": stop, "step": step},
        )
        ephemerides = obj.ephemerides()

        # store the ephemerides in the format we need
        self._ephemerides = []
        for row in range(len(ephemerides)):
            epoch = parse(ephemerides["datetime_str"][row]).replace(tzinfo=tzutc())
            ra = float(ephemerides["RA"][row]) * u.deg
            dec = float(ephemerides["DEC"][row]) * u.deg
            ra_rate = float(ephemerides["RA_rate"][row]) * u.arcsec / u.hour
            dec_rate = ephemerides["DEC_rate"][row] * u.arcsec / u.hour
            magnitude = ephemerides["V"][row] if "V" in ephemerides.keys() else 0
            magnitude_range = MagnitudeRange(
                min_magnitude=magnitude, max_magnitude=magnitude, bandpass="V"
            )
            self._ephemerides.append(
                Ephemeris(
                    ra=ra,
                    dec=dec,
                    ra_rate=ra_rate,
                    dec_rate=dec_rate,
                    magnitude_range=magnitude_range,
                    epoch=epoch,
                )
            )

    def ephemerides(self, start_time: datetime, end_time: datetime) -> List[Ephemeris]:
        if start_time.tzinfo is None or end_time.tzinfo is None:
            raise ValueError("The start and end time must be timezone-aware.")

        return _cover_time_interval(self._ephemerides, start_time, end_time)

    def is_sidereal_target(self) -> bool:
        return False

    def metadata(self) -> Metadata:
        return {
            "horizons_id": self.object_id,
            "horizons_stepsize": f"{self.stepsize.to_value(u.min)} min",
        }
