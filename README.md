# SALT Finder Charts

Generate finder charts for the Southern African Large Telescope (SALT).

## This package will change...

This is a preliminary and rudimentary version, and you should expect major changes in future versions.

## Installation

The package can be installed with pip.

```bash
pip install salt_finder_charts
```

## Usage

You can generate finder charts with the `standard_finder_chart` function. This returns a generator of binary streams containing the finder charts. For example, this is how you would generate finder charts for a longslit observation of an 18th magnitude target.

```python
import astropy.units as u
from salt_finder_charts import standard_finder_charts
from salt_finder_charts.mode import Mode
from salt_finder_charts.output import OutputFormat

fcs = standard_finder_charts(
    bandpass="V",
    dec=-27.63792 * u.deg,
    max_magnitude=18,
    min_magnitude=18,
    mode=Mode.LONGSLIT,
    output_format=OutputFormat.PDF,
    position_angle=0 * u.deg,
    ra=78.6568723 * u.deg,
    slitwidth=2 * u.arcsec,
    title="Test Finding Chart"
)
```

Why a generator? Isn't there a single finder chart only? Usually that's true. However, in case of an asteroid you may get multiple finder charts, in particular if the considered time interval covers multiple nights.

Saving the generated finder charts to disk is fairly straightforward.

```python
counter = 1
for fc in fcs:
    with open(f"FindingChart-{counter}.pdf", "wb") as f:
        f.write(fc.read())
        counter += 1
``` 

If you are certain there is only one finder chart, you could just use `next()` to get it.

```python
fc = next(fcs)
```

Here is how you could view finder charts for an asteroid. Note that both the pytz and the Pillow library must have been installed, and that you have to hit enter for moving to the next finder chart.

```python
from datetime import datetime
from PIL import Image

import pytz
import astropy.units as u
from salt_finder_charts import standard_finder_charts
from salt_finder_charts.mode import Mode
from salt_finder_charts.output import OutputFormat

utc = pytz.timezone("UTC")
start = utc.localize(datetime(2019, 6, 28, 11, 31, 0))
end = utc.localize(datetime(2019, 6, 29, 12, 31, 0))

for fc in standard_finder_charts(
    mode=Mode.IMAGING,
    horizons_id='54827',
    start_time=start,
    end_time=end,
    slitwidth=2 * u.arcsec,
    basic_annotations=False,
    output_format=OutputFormat.PNG
):
    image = Image.open(fc)
    image.show()
    input('Press enter to continue')
    image.close()
```

The `standard_finder_charts` function accepts the following arguments.

Argument | Description | Default
--- | --- | ---
bandpass | Bandpass for the magnitudes, such as V.
basic_annotations | Whether to add basic annotations only. | False
dec | Declination of the finder chart center, as an angle between -90 and 90 degrees.
end_time | End time until which to generate finder charts. This is only relevant for non-sidereal targets. Must be timezone-aware. | End of the current Julian day
horizons_id | Identifier for the Horizons service.
horizons_stepsize | Time between ephemerides queried from the Horizons service | 5 minutes
max_magnitude | Maximum magnitude
min_magnitude | Minimum magnitude for the target.
mode | Observation mode, such as imaging or longslit.
mos_mask_rsmt | Binary stream containing an RSMT MOS mask definition file.
output_format | Output format of the generated finder charts, such as PDF. | `OutputFormat.PDF`
position_angle | Position angle. | Generally 0&deg; or, if appropriate, chosen to allow for easier acquisition.
ra | Right ascension of the finder chart center, as an angle between 0 and 360 degrees.
slitwidth | Slit width, as an angle.
start_time | Start time from which to generate finder charts. This is only relevant for non-sidereal targets. | Start of the current Julian day. Must be timezone-aware.
survey | Survey from which to get the finder chart image.  | `Survey.POSS2UKSTU_RED`
title | Title of the finder chart.

Which of these arguments are required depends on the combination of arguments used. Also, some arguments may be ignored when used with others. For example, the slit width won;t be used if the observing mode is imaging.

While normally you wouldn't have to use it, the `salt_finder_charts` package also offers a `finder_charts` function which lets you customise some aspects of the finder chart generation.

For example, you could create the same longslit observation finder charts as above with the following code.

```python
import astropy.units as u
from salt_finder_charts import finder_charts
from salt_finder_charts.ephemerides import ConstantEphemerisService
from salt_finder_charts.image import SurveyImageService, Survey
from salt_finder_charts.mode import LongslitModeDetails
from salt_finder_charts.output import output_pdf
from salt_finder_charts.util import MagnitudeRange

magnitude_range = MagnitudeRange(bandpass="V", max_magnitude=18, min_magnitude=18)
ephemeris_service = ConstantEphemerisService(
    ra=78.6568723 * u.deg,
    dec=-27.63792 * u.deg,
    magnitude_range=magnitude_range
)
image_service = SurveyImageService(Survey.POSS2UKSTU_RED)
mode_details = LongslitModeDetails(slitwidth=2 * u.arcsec, pa=0 * u.deg)
fcs = finder_charts(
    ephemeris_service=ephemeris_service,
    image_service=image_service,
    mode_details=mode_details,
    output=output_pdf,
    title="test Finding Chart"
)
```

Like the `standard_finder_charts` function, `finder_charts` returns a generator of binary streams with the finder charts. It uses the following arguments.

Argument | Description | Default
--- | --- | ---
basic_annotations | Whether to add basic annotations only. | False
end_time | End time until which to generate finder charts. This is only relevant for non-sidereal targets. Must be timezone-aware. | End of the current Julian day
ephemeris_service | Service for getting ephemerides.
image_service | Service for getting finder chart images.
mode_details | Observing mode and its details
output | Function for converting a finding chart into an output format such as pdf.
start_time | Start time from which to generate finder charts. This is only relevant for non-sidereal targets.  Must be timezone-aware. | Start of the current Julian day.
title | Title for the finder chart.

## Command-line interface for generating finder charts

For convenience, a command `saltfc` is provided, which you can run in a terminal. It saves the generated finding charts in a directory. The filename of the generated files consists of a basename (such as `FinderChart`) followed by a dash and a running number. The basename can be customised with a command line option. Existing files are replaced without warning,

For example, you can generate finder charts for an asteroid by running the following command in a terminal.

```bash
saltfc \
    --output-dir /tmp \
    --basename "Asteroid_54827" \
    --mode imaging \
    --horizons-id 54827 \
    --start-time "2019-06-28 11:31:00" \
    --end-time "2019-06-29 12:31:00" \
    --output-format PDF
```

The command line options for `saltfc` are essentially the same as the arguments for the `standard_finder_charts` function, plus options for customising the file base name and setting the output directory.

Argument | Description | Default
--- | --- | ---
--bandpass | Bandpass for the magnitudes, such as V.
--basename | Basename for the saved finder chart files | FinderChart
--basic_annotations | Add basic annotations only. | False
--dec | Declination of the finder chart center, as an angle between -90 and 90 degrees, in degrees.
--end_time | End time until which to generate finder charts. This is only relevant for non-sidereal targets. The time is supposed to be in UTC. | End of the current Julian day
--horizons_id | Identifier for the Horizons service.
--horizons_stepsize | Minutes between ephemerides queried from the Horizons service | 5
--max_magnitude | Maximum magnitude for the target
--min_magnitude | Minimum magnitude for the target.
--mode | Observation mode, such as imaging or longslit.
--mos_mask_rsmt | RSMT MOS mask file
--output-dir | Directory where to store the generated finder chart files.
--output_format | Output format of the generated finder charts, such as PDF. | PDF
--position_angle | Position angle, in degrees. | Generally 0 or, if appropriate, chosen to allow for easier acquisition.
--ra | Right ascension of the finder chart center, as an angle between 0 and 360 degrees, in degrees.
--slitwidth | Slit width, as an angle in arcseconds.
--start_time | Start time from which to generate finder charts. This is only relevant for non-sidereal targets. | Start of the current Julian day. The time is supposed to be in UTC.
--survey | Survey from which to get the finder chart image.  | `Survey.POSS2UKSTU_RED`
--title | Title of the finder chart.

Both the start and end time use the format `yyyy-mm-dd hh:mm:ss`. For example, valid values are `2019-01-25 9:45:16'` or `2020-01-02 12:00:00`.

## Command interface for calculating position angles

A suitable position angle for a target can be calculated with the `pa` command. It tries to find a star in an annulus around the target with a magnitude in a given magnitude range. The stsr must have a given minimum separation from its neighbouring stars.

The radii, magnitudes and minimum separation all have default values, and usually you should not have to specify them. Hence a typical call of `pa` looks as follows./

```bash
pa --ra 217.524 --dec -23.97611
```

Here is a list of the available arguments for `pa`.

Argument | Description | Default
--- | --- | ---
--dec | Declination of the target, in degrees.
--max-mag | Maximum (faintest) magnitude a suitable star may have. | 18
--max-radius | Maximum radius, in arcminutes, of the annulus in which the suitable stars are located. " | 3
--min-mag | Minimum (brightest) magnitude a suitable star may have. | 15
--min-radius | Minimum radius of the annulus in which the suitable stars are located. | 1
--min-separation | Minimum angular distance, in arcseconds, a suitable star must have from its neighbouring stars. | 10
