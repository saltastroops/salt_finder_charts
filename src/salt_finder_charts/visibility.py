import bisect
from datetime import datetime, timedelta
from typing import cast, List, Optional, Tuple

from astropy.units import Quantity
from astropy import units as u
from dateutil.tz import tzutc
import ephem
import numpy as np

from salt_finder_charts.ephemerides import EphemerisService


SUTH_LONGITUDE = 20.8108 * u.deg

SUTH_LATITUDE = -32.3755556 * u.deg

MIN_ALTITUDE_SINE = np.sin(46.18 * u.deg)

MAX_ALTITUDE_SINE = np.sin(59.36 * u.deg)


def visibility_windows(
    ephemeris_service: EphemerisService, start_time: datetime, end_time: datetime
) -> List[Tuple[datetime, datetime]]:
    """
    Returns the visibility windows for all the nights in the
    specified time interval.

    Both the start and end time should be outside a night.

    Parameters
    ----------
    ephemeris_service : EphemerisService
        Service for getting the required ephemerides.
    start_time : datetime
        Start time (must be timezone-aware).
    end_time:  datetime
        End time (must be timezone-aware).

    Returns
    -------
    list of intervals
        The visibility windows.

    """

    # enforce timezones
    if start_time.tzinfo is None or end_time.tzinfo is None:
        raise ValueError("The start and end time must be timezone-aware.")

    if start_time >= end_time:
        raise Exception("The start time must be earlier than the end time.")

    t = start_time
    dt = timedelta(days=1)
    windows = []
    while t < end_time:
        windows_for_night = _visibility_windows_next_night(ephemeris_service, t)
        for window in windows_for_night:
            window_start = window[0]
            if window_start < end_time:
                window_end = min(window[1], end_time)
                windows.append((window_start, window_end))
        t += dt

    return windows


def fov_fitting_intervals(
    intervals: List[Tuple[datetime, datetime]],
    ephemeris_generator: EphemerisService,
    fov_radius: Quantity,
) -> List[Tuple[datetime, datetime]]:

    """
    Split time intervals up so that the positions within every interval fit into the
    field of view.

    Parameters
    ----------
    intervals : list of intervals
        Time intervals to split up.
    ephemeris_generator : EphemerisService
        Generator for the required ephemerides.
    fov_radius : Quantity
        The radius of the field of view, as an angle.

    Returns
    -------
    list of tuple
        The list of time intervals.

    """

    fitting_intervals: List[Tuple[datetime, datetime]] = []

    for interval in intervals:
        start = interval[0]
        end = interval[1]
        ephemerides = ephemeris_generator.ephemerides(start, end)
        epochs = [e.epoch for e in ephemerides]
        right_ascensions = [e.ra for e in ephemerides]
        declinations = [e.dec for e in ephemerides]

        # handle the transition from 360 to 0 degrees
        _transform_right_ascensions(right_ascensions, True)

        # find the RA, dec center
        center_ra, center_dec = EphemerisService.center_position(ephemerides)

        # are all positions located in the FOV around this center?
        all_in_fov = True
        for i in range(0, len(epochs)):
            if not is_in_fov(
                center_ra, center_dec, right_ascensions[i], declinations[i], fov_radius
            ):
                all_in_fov = False
                break
        if all_in_fov:
            fitting_intervals.append(interval)
        else:
            # split the interval in half (more or less)
            if len(epochs) <= 2:
                raise ValueError("The ephemerides are too sparse to cover the FOV.")
            td = epochs[-1] - epochs[0]
            dt = (
                td.microseconds + (td.seconds + td.days * 24 * 3600) * 10 ** 6
            ) / 10 ** 6
            dt_center = int(round(0.5 * dt))
            center_time = epochs[0] + timedelta(seconds=dt_center)
            index = bisect.bisect(epochs, center_time)
            if index == len(epochs) - 1:
                index -= 1
            interval1 = (epochs[0], epochs[index])
            interval2 = (epochs[index], epochs[-1])
            if interval1[0] != interval1[1]:
                fitting_intervals.extend(
                    fov_fitting_intervals([interval1], ephemeris_generator, fov_radius)
                )
            if interval2[0] != interval2[1]:
                fitting_intervals.extend(
                    fov_fitting_intervals([interval2], ephemeris_generator, fov_radius)
                )

    return fitting_intervals


def is_in_fov(
    center_ra: Quantity,
    center_dec: Quantity,
    ra: Quantity,
    dec: Quantity,
    fov_radius: Quantity,
) -> bool:
    """
    Check whether a right ascension and declination are within the field of view (FOV).

    Parameters
    ----------
    center_ra : Quantity
        Right ascension of the FOV's center, as an angle.
    center_dec : Quantity
        Declination of the FOV's center, as an angle.
    ra : Quantity
        Right ascension, as an angle.
    dec : Quantity
        Declination, as an angle
    fov_radius : Quantity
        Radius of the FOV, as an angle.

    Returns
    -------
    bool
        Whether the right ascension and declination are within the field of view.

    """

    dra = (ra - center_ra) * abs(np.cos(dec))
    ddec = dec - center_dec
    d = np.sqrt(dra ** 2 + ddec ** 2)

    return cast(bool, d <= fov_radius)


def _is_visible_with_salt(ra: Quantity, dec: Quantity, t: datetime) -> bool:
    """
    Checks whether a target is visible by SALT.

    Parameters
    ----------
    ra: Quantity
        Right ascension.
    dec: Quantity
        Declination.
    t: datetime
        Time (must be timezone-aware).

    """

    observer = _salt_observer(t)

    lst = observer.sidereal_time() * u.radian  # PyEphem uses radians
    hour_angle = lst - ra

    phi = SUTH_LATITUDE
    h = hour_angle

    sin_altitude = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(h)

    return cast(bool, MIN_ALTITUDE_SINE < sin_altitude < MAX_ALTITUDE_SINE)


def _visibility_windows_next_night(
    ephemeris_service: EphemerisService, t: datetime
) -> List[Tuple[datetime, datetime]]:
    """
    Returns the visibility windows for the night following
    the given time.

    Parameters
    ----------
    ephemeris_service : EphemerisService
        Service to use for getting the required ephemerides.
    t: datetime
        Time (must be timezone-aware).

    Returns
    -------
    list of tuple
        The visibility windows for the next night.

    """

    # get the night data
    observer = _salt_observer(t)

    sunset = observer.next_setting(ephem.Sun()).datetime().replace(tzinfo=tzutc())
    sunrise = observer.next_rising(ephem.Sun()).datetime().replace(tzinfo=tzutc())
    night_ephemerides = ephemeris_service.ephemerides(sunset, sunrise)

    # get the visibility windows
    windows: List[Tuple[datetime, datetime]] = []
    dt = timedelta(seconds=300)
    tp = sunset
    window_start: Optional[datetime] = None
    window_end: Optional[datetime] = None
    while tp <= sunrise:
        epochs = [e.epoch for e in night_ephemerides]
        index = bisect.bisect(epochs, tp)
        ra = night_ephemerides[index].ra
        dec = night_ephemerides[index].dec
        visible = _is_visible_with_salt(ra, dec, tp)
        if visible:
            if window_start is None:
                window_start = tp
            window_end = tp
        else:
            if window_start is not None and window_end is not None:
                windows.append((window_start - dt, window_end + dt))
                window_start = None
                window_end = None
        tp += dt
    if window_start is not None and window_end is not None:
        windows.append((window_start - dt, window_end + dt))

    return windows


def _salt_observer(t: datetime) -> ephem.Observer:
    """
    Returns a PyEphem Observer instance for the right ascension and declination of
    Sutherland.

    Parameters
    ----------
    t: datetime
        Datetime for the observer (must be timezone-aware).

    Returns
    -------
    Observer
        PyEphem Observer instance.

    """

    # Ensure timezone.
    if t.tzinfo is None:
        raise ValueError("The time must be timezone-aware.")

    observer = ephem.Observer()
    observer.lat = SUTH_LATITUDE.to_value(u.radian)
    observer.lon = SUTH_LONGITUDE.to_value(u.radian)
    observer.date = ephem.Date(t)

    return observer


def _transform_right_ascensions(
    right_ascensions: List[Quantity], continuous_at_360: bool
) -> None:
    """
    Transform the given right ascensions. If the continuousAt360 flag is true, values
    between 0 and 1 degree are increased by 360 degrees if there are values between 359
    and 360 degrees. Otherwise values greater than or equal to 360 degrees are reduced
    by 360 degrees.

    The transformation is done in place.

    This function assumes that the difference between subsequent right ascensions does
    not exceed 1 degree and that the overall range of right ascensions does not exceed a
    few degrees.

    Parameters
    ----------
    right_ascensions : list of Quantity
        List of right ascensions to transform.
    continuous_at_360 : bool
        Whether the transformed right ascensions may extend beyond 360 degrees.

    """

    if continuous_at_360:
        just_before_360 = False
        just_after_0 = False
        for r in right_ascensions:
            if r >= 359 * u.deg:
                just_before_360 = True
            if r <= 1 * u.deg:
                just_after_0 = True
        if just_before_360 and just_after_0:
            for i in range(0, len(right_ascensions)):
                if right_ascensions[i] <= 1 * u.deg:
                    right_ascensions[i] += 360 * u.deg
    else:
        for i in range(0, len(right_ascensions)):
            if right_ascensions[i] >= 360 * u.deg:
                right_ascensions[i] -= 360 * u.deg
