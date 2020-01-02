import click
from salt_finder_charts.image import Survey
from salt_finder_charts.mode import Mode
from salt_finder_charts.output import OutputFormat


@click.command()
@click.option('--bandpass', type=str, help='bandpass (such as V) for the magnitudes')
@click.option('--basename', type='str', default='FinderChart', help='Basename for the saved finder chart files.')
@click.option('--basic-annotations', is_flag=True, help='add basic annotations only')
@click.option('--dec', type=str, help='declination of the finder chart center')
@click.option('--end-time', type=click.DateTime(), help='emd time until which to generate finder charts')
@click.option('--horizons-id', type=str, help='identifier for the Horizons service')
@click.option('--horizons-stepsize', type=int, default=5, help='minutes between ephemerides queried from the Horizoms service')
@click.option('--max-magnitude', type=float, help='maximum magnitude of the target')
@click.option('--min-magnitude', type=float, help='minimum magnitude of the target')
@click.option('--mode', type=click.Choice([mode.value for mode in Mode], case_sensitive=False), required=True, help='observation mode')
@click.option('--mos-mask-rsmt', type=click.File(mode='rb'), help='RSMT file defining a MOS mask')
@click.option('--name', type=str, default='FinderChart', help='base name for the generated finder chart files')
@click.option('--output-dir', type=click.Path(exists=True, file_okay=False, writable=True, resolve_path=True), required=True, help='directory where to save the generated finder chart files')
@click.option('--output-format', type=click.Choice([of.value for of in OutputFormat],case_sensitive=False), default='PDF', help='output format of the generated finder chart files')
@click.option('--position-angle', type=float, help='position angle in degrees')
@click.option('--ra', type=str, help='right ascension of the finder chart center')
@click.option('--slitwidth', type=float, help='slit width in arcseconds')
@click.option('--start-time', type=click.DateTime(formats=['%Y-%m-%d %H:%M:%S']), help='start time from when to generate finder charts')
@click.option('--survey', type=click.Choice([survey.value for survey in Survey], case_sensitive=False), default='POSS2/UKSTU Red', help='survey to use for the findder chart image')
@click.option('--title', type=str, help='title for the finder chart')
def saltfc(bandpass, basic_annotations, dec, end_time, horizons_id, horizons_stepsize, max_magnitude, min_magnitude, mode, mos_mask_rsmt, name, output_dir, output_format, position_angle, ra, slitwidth, start_time, survey, title):
    """
    Command for generating SALT finder charts.

    By default, the finder charts are stored as files named FinderChart-1,
    FinderChart-2, ... (with the correct file suffix according to the chosen output
    format), but you can change the base name ("FindingChart") with the --name
    parameter. No running number is added if only one finder chart is generated. The
    target directory for the finder charts must be specified with the --dir parameter.

    See the README file for a discussion of the various parameters.

    """

    click.echo("HELLO")
