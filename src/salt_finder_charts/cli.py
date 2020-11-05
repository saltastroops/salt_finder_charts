from typing import Optional, BinaryIO

import astropy.units as u
import click
from datetime import datetime
import os
import pytz
from salt_finder_charts import standard_finder_charts, __version__
from salt_finder_charts.image import Survey
from salt_finder_charts.mode import Mode
from salt_finder_charts.output import OutputFormat
from salt_finder_charts.position_angle import (
    estimated_position_angle,
    MAX_MAG,
    MAX_RADIUS,
    MIN_MAG,
    MIN_RADIUS,
    MIN_STAR_SEPARATION,
)
from salt_finder_charts.util import julian_day_start, julian_day_end


@click.command()
@click.option("--bandpass", type=str, help="bandpass (such as V) for the magnitudes")
@click.option(
    "--basename",
    type=str,
    default="FinderChart",
    help="Basename for the saved finder chart files.",
)
@click.option("--basic-annotations", is_flag=True, help="add basic annotations only")
@click.option("--dec", type=float, help="declination of the finder chart center")
@click.option(
    "--end-time",
    type=click.DateTime(),
    help="emd time until which to generate finder charts",
)
@click.option("--horizons-id", type=str, help="identifier for the Horizons service")
@click.option(
    "--horizons-stepsize",
    type=int,
    default=5,
    help="minutes between ephemerides queried from the Horizoms service",
)
@click.option("--max-magnitude", type=float, help="maximum magnitude of the target")
@click.option("--min-magnitude", type=float, help="minimum magnitude of the target")
@click.option(
    "--mode",
    type=click.Choice([mode.value for mode in Mode], case_sensitive=False),
    required=True,
    help="observation mode",
)
@click.option(
    "--mos-mask-rsmt", type=click.File(mode="rb"), help="RSMT file defining a MOS mask"
)
@click.option(
    "--output-dir",
    type=click.Path(exists=True, file_okay=False, writable=True, resolve_path=True),
    required=True,
    help="directory where to save the generated finder chart files",
)
@click.option(
    "--output-format",
    type=click.Choice([of.value for of in OutputFormat], case_sensitive=False),
    default="PDF",
    help="output format of the generated finder chart files",
)
@click.option("--position-angle", type=float, help="position angle in degrees")
@click.option("--ra", type=float, help="right ascension of the finder chart center")
@click.option("--slitwidth", type=float, help="slit width in arcseconds")
@click.option(
    "--start-time",
    type=click.DateTime(formats=["%Y-%m-%d %H:%M:%S"]),
    help="start time from when to generate finder charts",
)
@click.option(
    "--survey",
    type=click.Choice([survey.value for survey in Survey], case_sensitive=False),
    default="POSS2/UKSTU Red",
    help="survey to use for the finder chart image",
)
@click.option("--title", type=str, help="title for the finder chart")
@click.version_option(__version__)
def saltfc(
    bandpass: Optional[str],
    basename: str,
    basic_annotations: bool,
    dec: Optional[float],
    end_time: Optional[datetime],
    horizons_id: Optional[str],
    horizons_stepsize: Optional[int],
    max_magnitude: Optional[float],
    min_magnitude: Optional[float],
    mode: str,
    mos_mask_rsmt: Optional[BinaryIO],
    output_dir: str,
    output_format: str,
    position_angle: Optional[float],
    ra: Optional[float],
    slitwidth: Optional[float],
    start_time: Optional[datetime],
    survey,
    title,
):
    """
    Command for generating SALT finder charts.

    By default, the finder charts are stored as files named FinderChart-1,
    FinderChart-2, ... (with the correct file suffix according to the chosen output
    format), but you can change the base name ("FindingChart") with the --name
    parameter. No running number is added if only one finder chart is generated. The
    target directory for the finder charts must be specified with the --dir parameter.

    See the README file for a discussion of the various parameters.

    """

    # start and end time
    if start_time:
        _start_time = pytz.utc.localize(start_time)
    else:
        _start_time = julian_day_start(datetime.now(pytz.utc))
    if end_time:
        _end_time = pytz.utc.localize(end_time)
    else:
        _end_time = julian_day_end(datetime.now(pytz.utc))

    # finder chart center
    _ra = ra * u.deg if ra is not None else None
    _dec = dec * u.deg if dec else None

    # position angle
    _position_angle = position_angle * u.deg if position_angle is not None else None

    # mode
    _mode = [m for m in Mode if m.value.lower() == mode.lower()][0]

    # output
    _output_format = [
        of for of in OutputFormat if of.value.lower() == output_format.lower()
    ][0]

    # survey
    _survey = [s for s in Survey if s.value.lower() == survey.lower()][0]

    # slit width
    _slitwidth = slitwidth * u.arcsec if slitwidth is not None else None

    # Horizons query stepsize
    _horizons_stepsize = (
        horizons_stepsize * u.minute if horizons_stepsize is not None else None
    )

    # generate the finder charts
    counter = 1
    for fc in standard_finder_charts(
        mode=_mode,
        output_format=_output_format,
        start_time=_start_time,
        end_time=_end_time,
        ra=_ra,
        dec=_dec,
        min_magnitude=min_magnitude,
        max_magnitude=max_magnitude,
        bandpass=bandpass,
        horizons_id=horizons_id,
        horizons_stepsize=_horizons_stepsize,
        survey=_survey,
        position_angle=_position_angle,
        slitwidth=_slitwidth,
        mos_mask_rsmt=mos_mask_rsmt,
        basic_annotations=basic_annotations,
        title=title,
    ):
        filename = f"{basename}-{counter}.{_output_format.extension()}"
        counter += 1
        filepath = os.path.join(output_dir, filename)
        with open(filepath, "wb") as f:
            f.write(fc.read())


@click.command()
@click.option(
    "--dec", type=float, required=True, help="declination of the target, in degrees"
)
@click.option("--max-mag", type=float, default=MAX_MAG, help="maximum magnitude")
@click.option(
    "--max-radius",
    type=float,
    default=MAX_RADIUS.to_value(u.arcmin),
    help="maximum radius, in arcminutes",
)
@click.option("--min-mag", type=float, default=MIN_MAG, help="minimum magnitude")
@click.option(
    "--min-radius",
    type=float,
    default=MIN_RADIUS.to_value(u.arcmin),
    help="minimum radius, in arcminutes",
)
@click.option(
    "--min-separation",
    type=float,
    default=MIN_STAR_SEPARATION.to_value(u.arcsec),
    help="minimum separation between stars, in arcseconds",
)
@click.option(
    "--ra", type=float, required=True, help="right ascension of the target, in degrees"
)
def pa(ra: float, dec: float, min_radius, max_radius, min_mag, max_mag, min_separation):
    """
    Calculate a suitable position angle for a target.

    The position angle is the position angle of a star within an annulus between a
    minimum and maximum radius. This star must have a magnitude between a minimum and
    maximum magnitude, and it must have a minimum separation from neighbouring stars.

    The position angle is output as an angle in degrees.

    """

    print(
        estimated_position_angle(
            ra=ra * u.deg,
            dec=dec * u.deg,
            radius_range=(min_radius * u.arcmin, max_radius * u.deg),
            mag_range=(min_radius, max_radius),
            min_star_separation=min_separation * u.arcsec,
        )
    )
