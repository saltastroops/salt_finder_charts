from datetime import datetime, timedelta
from typing import BinaryIO, Generator, Optional, Tuple

import astropy.units as u
from astropy.units import Quantity
from dateutil.tz import tzutc

from salt_finder_charts.image import Survey, SurveyImageService
from salt_finder_charts.mode import Mode, ModeDetails, ImagingModeDetails, \
    LongslitModeDetails, SlotModeDetails, MOSModeDetails
from salt_finder_charts.output import output_pdf, output_png, OutputFormat
from salt_finder_charts.util import MagnitudeRange, MOSMask
from salt_finder_charts import finder_charts
from salt_finder_charts.ephemerides import HorizonsEphemerisService, \
    ConstantEphemerisService, EphemerisService

TimeInterval = Tuple[datetime, datetime]


def standard_finder_charts(
        # arguments which are always required
        mode: Mode,
        basic_annotations: bool,
        output_format: OutputFormat,
        # time interval
        start_time: Optional[datetime]=None,
        end_time: Optional[datetime]=None,
        # ephemerides
        ra: Optional[Quantity]=None,
        dec: Optional[Quantity]=None,
        min_magnitude: Optional[float]=None,
        max_magnitude: Optional[float]=None,
        bandpass: Optional[str]=None,
        horizons_id: Optional[str]=None,
        horizons_stepsize: Optional[Quantity]=None,
        # image
        survey: Optional[Survey]=Survey.POSS2UKSTU_RED,
        # instrument mode details
        position_angle: Optional[Quantity]=None,
        slitwidth: Optional[Quantity]=None,
        mos_mask_rsmt: Optional[BinaryIO]=None,
        # miscellaneous
        title: Optional[str]=None,
) -> Generator[BinaryIO, None, None]:
    """
    Create standard SALT finder charts.

    Some of the parameters are mutually exclusive. For example, it does mot make sense
    to specify a slit width if you generate finding charts for imaging mode. In some
    cases such combinations will raise an error, but in others some of the parameters
    may just be ignored.

    If no start time is given, the beginning of the current Julian day is assumed. If no
    end time is given, the end of the current Julian day is assumed.

    Parameters
    ----------
    mode : Mode
        Observation mode (such as imaging or MOS).
    basic_annotations : bool
        Whether only basic annotations should be added to the finder chart.
    output_format : OutputFormat
        Output format (such as PDF) to use for the generated finder charts.
    start_time : datetime
        Start time from which to generate finder charts.
    end_time : datetime
        End time until which to generate finder charts.
    ra : Quantity
        Right ascension of the finder chart center.
    dec : Quantity
        Declination of the finder chart center.
    min_magnitude : float
        Minimum magnitude of the target.
    max_magnitude L: float
        Maximum magnitude of the target.
    bandpass : str
        Bandpass (such as V) for the magnitudes,
    horizons_id : str
        Identifier for a target in the Horizons database.
    horizons_stepsize : Quantity
        Time between ephemerides queried from the Horizons service. The default is 5
        minutes.
    survey : Survey
        The image survey from which the finder chart image shall be taken.
    position_angle : Quantity
        The position angle.
    slitwidth : Quantity
        The width of the longslit, as an angle.
    mos_mask_rsmt : BinaryIO
        Input stream containing an RSMT file for a MOS setup.
    title : str
        Title for the finder chart.

    Returns
    -------
    Generator of BinaryIO
        The finder charts as input streams.

    """

    # time interval

    # get default start and end time if need be
    if not start_time:
        start_time = _default_start_time()
    if not end_time:
        end_time = _default_end_time()

    # ensure there are timezones
    if start_time.tzinfo is None:
        raise ValueError('The start time must be timezone-aware.')
    if end_time.tzinfo is None:
        raise ValueError('The end time must be timezone aware.')

    # ephemerides

    mos_mask: Optional[MOSMask] = None
    if mode == Mode.MOS:
        if mos_mask_rsmt is None:
            raise ValueError("A RSMT file must be supplied if a finding chart is generated for MOS mode.")
        if ra or dec or position_angle:
            raise ValueError("You must not supply a right ascension, declination or position angle in MOS mode, as they are taken from the MOS mask definition.")
        mos_mask = MOSMask(mos_mask_rsmt)
        ra = mos_mask.right_ascension
        dec = mos_mask.declination
        position_angle = mos_mask.position_angle

    if horizons_id:
        # get ephemerides from Horizons
        if ra is not None or dec is not None:
            raise ValueError("No right ascension or declination must be supplied if a Horizons identifier is supplied.")
        if horizons_stepsize is None:
            horizons_stepsize = 5 * u.minute
        ephemeris_service: EphemerisService = HorizonsEphemerisService(object_id=horizons_id, start_time=start_time - timedelta(days=2), end_time=end_time + timedelta(days=2), stepsize=horizons_stepsize)
    else:
        # use ephemerides for a non-sidereal target
        if ra is None:
            raise ValueError("The right ascension is missing.")
        if dec is None:
            raise ValueError("The declination is missing.")
        if min_magnitude is not None and (max_magnitude is None or bandpass is None):
            raise ValueError("You must supply a maximum magnitude and bandpass if you supply a minimum magnitude.")
        if max_magnitude is not None and (min_magnitude is None or bandpass is None):
            raise ValueError("You must supply a minimum magnitude and bandpass if you supply a maximum magnitude.")
        if bandpass is not None and (min_magnitude is None or max_magnitude is None):
            raise ValueError("You must supply a minimum and maximum magnitude if you supply a bandpass.")
        magnitude_range: Optional[MagnitudeRange] = None
        if min_magnitude is not None and max_magnitude is not None and bandpass is not None:
            magnitude_range = MagnitudeRange(min_magnitude=min_magnitude, max_magnitude=max_magnitude, bandpass=bandpass)
        ephemeris_service = ConstantEphemerisService(ra=ra, dec=dec, magnitude_range=magnitude_range)

    # image

    image_service = SurveyImageService(survey=survey)

    # mode details

    if mode is None:
        raise ValueError("You must specify an instrument mode.")
    if mode == Mode.IMAGING or mode == Mode.HRS:
        mode_details: ModeDetails = ImagingModeDetails(position_angle)
    elif mode == Mode.SLOT:
        mode_details = SlotModeDetails(pa=position_angle)
    elif mode == Mode.LONGSLIT:
        if slitwidth is None:
            raise ValueError("A slit width is required if a finding chart is generated for longslit mode.")
        mode_details = LongslitModeDetails(slitwidth=slitwidth, pa=position_angle)
    elif mode == Mode.MOS:
        mode_details = MOSModeDetails(mos_mask)
    else:
        raise ValueError(f"Mode unsupported: {mode.value}")

    # output

    if output_format == OutputFormat.PDF:
        output = output_pdf
    elif output_format == OutputFormat.PNG:
        output = output_png
    else:
        raise ValueError(f'Output format unsupported: {output_format.value}')

    # generate the finder charts
    return finder_charts(mode_details=mode_details,
                         start_time=start_time,
                         end_time=end_time,
                         ephemeris_service=ephemeris_service,
                         image_service=image_service,
                         title=title,
                         basic_annotations=basic_annotations,
                         output = output)


def _default_start_time() -> datetime:
    """
    Return the default start time.

    If the current time is noon (UTC) or later this is noon (UTC) of the current date,
    otherwise it is noon (UTC) of the previous day.

    In other words, the start time is the beginning of the current Julian day.

    Returns
    -------
    datetime
        Default start time.

    """

    now = datetime.now(tz=tzutc())
    noon = now.replace(hour=12, minute=0, second=0, microsecond=0)
    if now < noon:
        return noon - timedelta(days=1)
    else:
        return noon


def _default_end_time() -> datetime:
    """
    Return the default end time.

    If the current time is noon (UTC) or earlier this is noon (UTC) of the current date,
    otherwise it is noon (UTC) of the following day.

    In other words, the end time is the end of the current Julian day.

    Returns
    -------
    datetime
        Default end time.

    """

    now = datetime.now(tz=tzutc())
    noon = now.replace(hour=12, minute=0, second=0, microsecond=0)
    if now < noon:
        return noon
    else:
        return noon + timedelta(days=1)
