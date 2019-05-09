# SALT Finder Charts

Generate finder charts for the Southern African Large Telescope (SALT).

## This API will change...

This is a preliminary and rudimentary version, and you should expect major changes in future versions.

## Installation

The package can be installed with pip.

```bash
pip install salt_finder_charts
```

## Usage

You can generate a finder chart with the `finder_chart` function.

```python
from salt_finder_charts import finder_chart, Mode, Survey, Target

# Replace with a file path of your choice.
# The file extension must be one of those supported by the Pillow package.
out =  "/path/to/finder_chart.pdf"

target = Target(ra=13.1867, dec=-72.828611, name="SMC")
fc = finder_chart(target=target, mode=Mode.IMAGING, pa=15, survey=Survey.POSS2UKSTU_RED)
fc.save(out)
```

The function takes a mode, target details, a position angle and (optionally) a survey as its arguments and returns an `Image` object. Refer to the [documentation for the Pillow package](https://pillow.readthedocs.io/en/stable/) for details on how to use this object.

## Input arguments

The `finder_chart` function expects target details, a mode, a position angle and (optionally) a survey as its input arguments.

### Target

The target details should be supplied as an instance of the `Target` class. This is a named tuple with the following properties.

Property | Explanation | Example
--- | --- | ---
`dec` | Declination, in degrees between -90 and 90 | `16.9235`
`name` | Target name | `"My Interesting Target"`
`ra` | Right ascension, in degrees between 0 and 360 | `-57.36229`

The target name is included in the finder chart's title.

### Mode

The mode should be one of the enumeration members of the `Mode` enumeration, as listed in the following table.

Enumeration value | Explanation
--- | ---
`Mode.IMAGING` | An imaging observation.

### Position angle

The position angle must be a valid angle in degrees.

### Survey

The survey should be one of the enumeration members of the `Survey` enumeration. The following surveys are available.

Enumeration value | Survey
--- | ---
`Survey.POSS1_BLUE` | POSS1 Blue
`Survey.POSS1_RED` | POSS1 Red
`Survey.POSS2UKSTU_BLUE` | POSS2/UKSTU Blue
`Survey.POSS2UKSTU_IR` | POSS2/UKSTU IR
`Survey.POSS2UKSTU_RED` | POSS2/UKSTU Red
`Survey.TWO_MASS_H` | 2MASS-H
`Survey.TWO_MASS_J` | 2MASS-J
`Survey.TWO_MASS_K` | 2MASS-K

The POSS2/UKSTU Red survey is used if no survey is supplied.

Not all of the surveys cover all of the declinations observable with SALT, and you may get cryptic errors if you try to use a survey with a target position not covered by it. As such it is advisable to use the default whenever possible.
