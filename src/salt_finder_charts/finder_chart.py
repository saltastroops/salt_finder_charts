import sys
from datetime import datetime

from astropy.units import Quantity
from astropy import units as u
import numpy as np
import aplpy
import pytz
from typing import Any, BinaryIO, Callable, Generator, List, Optional, Tuple, Union

from salt_finder_charts import __version__
from salt_finder_charts.ephemerides import EphemerisService, Ephemeris
from salt_finder_charts.image import ImageService
from salt_finder_charts.mode import ModeDetails, Mode
from salt_finder_charts import visibility
from salt_finder_charts.util import julian_day_start, julian_day_end, Metadata


class FinderChart:
    """
    A finder chart for SALT.

    Parameters
    ----------
    start_time : datetime
        Start time from when the finder chart should be valid (must be timezone-aware).
    end_time : datetime
        End time until which the finder chart should be valid (must be timezone-aware).
    ephemeris_service : EphemerisService
        Service to use for getting the required ephemerides.
    image_Service : ImageService
        Service to use for getting the finder chart image.

    Attributes
    ----------
    metadata : FinderChartMetadata
        Metadata for this finder chart.

    """

    # Assume B to be the box just covering the target path. Then
    # MINIMUM_PATH_BOX_WIDTH is the minimum width of B for
    # which the path should be plotted.
    MINIMUM_PATH_BOX_WIDTH = 16.0 * u.arcsec

    def __init__(
        self,
        start_time: datetime,
        end_time: datetime,
        ephemeris_service: EphemerisService,
        image_service: ImageService,
        mode_details: ModeDetails,
        title: Optional[str],
        basic_annotations: bool = False,
    ):
        # enforce timezones
        if start_time.tzinfo is None or end_time.tzinfo is None:
            raise ValueError("The start and end time must be timezone-aware.")

        self.start_time = start_time
        self.end_time = end_time
        self.ephemeris_service = ephemeris_service
        self.ephemerides = self.ephemeris_service.ephemerides(start_time, end_time)
        self.ra, self.dec = EphemerisService.center_position(self.ephemerides)
        self.pa = mode_details.position_angle()
        self.magnitude_range = EphemerisService.find_magnitude_range(self.ephemerides)
        self.image_generator = image_service
        self.mode_details = mode_details
        self.title = title
        self.basic_annotations = basic_annotations

        self._init_plot()
        if not ephemeris_service.is_sidereal_target():
            self._draw_ephemeris_info()

    @property
    def metadata(self) -> Metadata:
        return dict(
            version=__version__,
            right_ascension=f"{self.ra.to_value(u.deg)} deg",
            declination=f"{self.dec.to_value(u.deg)} deg",
            min_magnitude=self.magnitude_range.min_magnitude
            if self.magnitude_range
            else None,
            max_magnitude=self.magnitude_range.max_magnitude
            if self.magnitude_range
            else None,
            bandpass=self.magnitude_range.bandpass if self.magnitude_range else None,
            position_angle=f"{self.pa.to_value(u.deg)} deg",
            start_time=self.start_time.isoformat(),
            end_time=self.end_time.isoformat(),
            mode=self.mode_details.mode.value,
            mode_metadata=self.mode_details.metadata(),
            ephemeris_metadata=self.ephemeris_service.metadata(),
        )

    def discard(self) -> None:
        """
        Free the plot resources.

        """

        self.plot.close()

    def draw_circle(
        self, center_ra: Quantity, center_dec: Quantity, radius: Quantity, color: str
    ) -> None:
        """
        Draw a circle.

        Parameters
        ----------
        center_ra : Quantity
            Right ascension of the circle's center., as an angle
        center_dec : Quantity
            Declination of the circle's center, as an angle.
        radius : Quantity
            Radius of the circle, as an angle.
        color : str
            Line color.

        """

        self.plot.show_circles(
            [center_ra.to_value(u.deg)],
            [center_dec.to_value(u.deg)],
            [radius.to_value(u.deg)],
            edgecolor=color,
        )

    def draw_centered_line(
        self,
        theta: Quantity,
        length: Quantity,
        ra: Quantity,
        dec: Quantity,
        color: str = "b",
        linewidth: float = 1.0,
        alpha: float = 0.7,
    ) -> None:
        """
        Draw a line with a given center point and angle.

        The values for the color, line width and opacity are the same you would use
        when drawing on a FITSFigure instance.

        Parameters
        ----------
        theta : Quantity
            Angle of the line.
        length : Quantity
            Length of the line, as an angle.
        ra : Quantity
            Right ascension of the line's center, as an angle.
        dec : Quantity
            Declination of the line's center, as an angle.
        color : str
            Line color.
        linewidth : str
            Line width.
        alpha : float
            Opacity.

        """

        _length = length / 2.0
        dx = np.sin(theta) * _length / np.cos(dec)
        dy = np.cos(theta) * _length
        coords = np.array(
            [
                [(ra + dx).to_value(u.deg), (ra - dx).to_value(u.deg)],
                [(dec + dy).to_value(u.deg), (dec - dy).to_value(u.deg)],
            ]
        )
        self.plot.show_lines([coords], color=color, linewidth=linewidth, alpha=alpha)

    def draw_centered_rectangle(
        self,
        theta: Quantity,
        width: Quantity,
        height: Quantity,
        ra: Quantity,
        dec: Quantity,
        color: str,
        linewidth: float = 1,
        alpha: float = 0.7,
    ) -> None:
        """
        Draw a rectangle with a given center, width, height and rotation angle around
        the center.

        The values for the color, line width and opacity are the same you would use
        when drawing on a FITSFigure instance.

        Parameters
        ----------
        theta : Quantity
            Rotation angle.
        width : Quantity
            Width of the rectangle, as an angle.
        height : Quantity
            Height of the rectangle, as an angle.
        ra : Quantity
            Right ascension of the rectangle's center, as an angle.
        dec : Quantity
            Declination of the rectangle's center, as an angle.
        color : str
            Line color.
        linewidth : float
            Line width.
        alpha : float
            Opacity.

        """

        height /= 2.0
        width /= 2.0
        # position of line centers
        ra_l = ra + np.cos(theta) * width / np.cos(dec)
        ra_r = ra - np.cos(theta) * width / np.cos(dec)
        dec_l = dec - np.sin(theta) * width
        dec_r = dec + np.sin(theta) * width

        # "width" in x and y direction
        dx = np.sin(theta) * height / np.cos(dec)
        dy = np.cos(theta) * height

        # use floats
        ra_l = ra_l.to_value(u.deg)
        ra_r = ra_r.to_value(u.deg)
        dec_l = dec_l.to_value(u.deg)
        dec_r = dec_r.to_value(u.deg)
        dx = dx.to_value(u.deg)
        dy = dy.to_value(u.deg)

        coords = np.array(
            [
                [ra_l, ra_l + dx, ra_r + dx, ra_r - dx, ra_l - dx, ra_l],
                [dec_l, dec_l + dy, dec_r + dy, dec_r - dy, dec_l - dy, dec_l],
            ]
        )
        self.plot.show_lines([coords], color=color, linewidth=linewidth, alpha=alpha)

    def draw_label(
        self,
        x: Union[Quantity, float],
        y: Union[Quantity, float],
        text: str,
        **kwargs: Any,
    ) -> None:
        """
        Draw a text label.

        As for labels in a FITSFigure, the position of the label may be give either as
        right ascension and declination (as an angle) or in plot coordinates. A value
        is interpreted as a right ascension or declination (rather than a plot
        coordinate) if it is a Quantity instance.

        Parameters
        ----------
        x : Quantity or float
            Right ascension (as an angle) or x coordinate in plot coordinates.
        y : Quantity of float
            Declination (as an angle) or y coordinate in plot coordinates.
        text : str
            Label text
        **kwargs
            Other keyword arguments accepted by the add_label method of the FITSFigure
            class.

        """

        if isinstance(x, Quantity) and isinstance(y, Quantity):
            _x = x.to_value(u.deg)
            _y = y.to_value(u.deg)
            relative = False
        else:
            _x = x
            _y = y
            relative = True
        self.plot.add_label(_x, _y, text, relative=relative, **kwargs)

    def _init_plot(self) -> None:
        """
        Set up the basic finder chart content.

        """

        # create a grayscale plot
        out = sys.stdout
        sys.stdout = open("/dev/null", "w")
        hdu = self.image_generator.image(self.ra, self.dec)
        self.plot = aplpy.FITSFigure(hdu)
        self.plot.show_grayscale()
        self.plot.set_theme("publication")
        sys.stdout = out

        # label for the position angle
        pa_string = "PA = %.1f" % self.mode_details.position_angle().to_value(u.deg)
        if self.mode_details.automated_position_angle():
            pa_string += " (auto)"
        self.draw_label(0.95, -0.05, pa_string, style="italic", weight="bold")

        # label for the title
        if self.title:
            self.draw_label(
                0.5, 1.03, self.title, style="italic", weight="bold", size="large"
            )

        # label for the image source
        self.draw_label(
            -0.05,
            -0.05,
            "%s" % self.image_generator.source(),
            style="italic",
            weight="bold",
        )

        # grid overlay
        self.plot.add_grid()
        self.plot.grid.set_alpha(0.2)
        self.plot.grid.set_color("b")

        # indicate the RSS field of view
        self.draw_circle(self.ra, self.dec, 4.0 * u.arcmin, "g")
        self.draw_label(
            0.79,
            0.79,
            "RSS",
            style="italic",
            weight="bold",
            size="large",
            horizontalalignment="left",
            color=(0, 0, 1),
        )

        # indicate the Salticam field of view
        self.draw_circle(self.ra, self.dec, 5.0 * u.arcmin, "g")
        self.draw_label(
            0.86,
            0.86,
            "SCAM",
            style="italic",
            weight="bold",
            size="large",
            horizontalalignment="left",
            color=(0, 0, 1),
        )

        # labels for north and east direction
        self.draw_label(
            self.ra,
            self.dec + 4.8 * u.arcmin,
            "N",
            style="italic",
            weight="bold",
            size="large",
            color=(0, 0.5, 1),
        )
        self.draw_label(
            self.ra + 4.8 * u.arcmin / np.abs(np.cos(self.dec)),
            self.dec,
            "E",
            style="italic",
            weight="bold",
            size="large",
            horizontalalignment="right",
            color=(0, 0.5, 1),
        )

        # add cross hairs
        self.draw_centered_line(
            0 * u.deg,
            8 * u.arcmin,
            self.ra,
            self.dec,
            color="g",
            linewidth=0.5,
            alpha=1.0,
        )
        self.draw_centered_line(
            90 * u.deg,
            8 * u.arcmin,
            self.ra,
            self.dec,
            color="g",
            linewidth=0.5,
            alpha=1.0,
        )

        # label for the magnitude range and bandpass
        if self.magnitude_range:
            self._show_magnitudes()

        # add mode specific content
        if not self.basic_annotations:
            self.mode_details.annotate_finder_chart(self)

    def _draw_ephemeris_info(self) -> None:
        """
        Draw ephemeris details.

        The drawn details depend on the target's movement. If the whole target path for
        the finder chart finds into a square of length MINIMUM_PATH_BOX_WIDTH, the
        following is added to the plot:

        - an arrow head at the the path's center indicating the direction of motion
        - a circle around this arrow head
        - the time interval for which the finder chart is valid

        Otherwise the following is added:

        - a line indicating the path of the target
        - an arrow head at the end of this path
        - the start time from which and the ebd time until which the finder chart is
          valid, added near the start and end of the path, respectively

        Only part of the details are shown if the finder chart contains basic
        annotations only.

        This method should only be used for non-sidereal targets.

        """

        ephemerides = self.ephemerides
        basic_annotations = self.basic_annotations
        center_ra, center_dec = EphemerisService.center_position(ephemerides)

        # is the target moving much?
        ra_min = min(ephemerides, key=lambda e: e.ra).ra
        ra_max = max(ephemerides, key=lambda e: e.ra).ra
        dec_min = min(ephemerides, key=lambda e: e.dec).dec
        dec_max = max(ephemerides, key=lambda e: e.dec).dec

        ra_width = ra_max - ra_min
        dec_width = dec_max - dec_min
        if (
            ra_width > FinderChart.MINIMUM_PATH_BOX_WIDTH
            or dec_width > FinderChart.MINIMUM_PATH_BOX_WIDTH
        ):
            significant_movement = True
        else:
            significant_movement = False

        # we have to convert angles to floats, as the vstack function does not accept
        # Quantity values
        right_ascensions = [e.ra for e in ephemerides]
        declinations = [e.dec for e in ephemerides]
        epochs = [e.epoch for e in ephemerides]
        start_time = epochs[0]
        end_time = epochs[-1]

        dra_start_to_end = right_ascensions[-1] - right_ascensions[0]
        ddec_start_to_end = declinations[-1] - declinations[0]

        # plot the target's path
        if significant_movement and not basic_annotations:
            # we have to convert angles to floats, as the vstack function does not accept
            # Quantity values
            right_ascensions_deg = [e.ra.to_value(u.deg) for e in ephemerides]
            declinations_deg = [e.dec.to_value(u.deg) for e in ephemerides]
            lv = np.vstack([right_ascensions_deg, declinations_deg])
            self.plot.show_lines(
                [lv], layer="object_path_lines", color="b", linewidth=1, alpha=1
            )

        # direction at the start and end
        ddec_start = declinations[1] - declinations[0]
        dra_end = right_ascensions[-1] - right_ascensions[-2]
        ddec_end = declinations[-1] - declinations[-2]

        if not basic_annotations:
            if significant_movement:
                # plot the arrow at the end time to show the direction
                self._draw_arrow_head(
                    right_ascensions[-1], declinations[-1], dra_end, ddec_end
                )
            else:
                ra_correction = abs(np.cos(center_dec))
                v_x, v_y = dra_start_to_end * ra_correction, ddec_start_to_end
                length = np.sqrt(v_x ** 2 + v_y ** 2)
                v_x, v_y = (
                    v_x.to_value(u.deg) / length.to_value(u.deg),
                    v_y.to_value(u.deg) / length.to_value(u.deg),
                )
                self._draw_arrow_head(
                    center_ra + 0.0013 * u.deg * v_x / ra_correction,
                    center_dec + 0.0013 * u.deg * v_y,
                    dra_start_to_end,
                    ddec_start_to_end,
                )

            # the labels shouldn't overlap with the path
            abs_vertical_shift = 0.002 * u.deg
            if significant_movement:
                label_position_start = {
                    "horizontal_alignment": "center",
                    "horizontal_position": right_ascensions[0],
                    "vertical_alignment": "top" if ddec_start > 0 else "bottom",
                    "vertical_position": declinations[0],
                    "vertical_shift": (-1 if ddec_start > 0 else 1)
                    * abs_vertical_shift,
                }
                label_position_end = {
                    "horizontal_alignment": "center",
                    "horizontal_position": right_ascensions[-1],
                    "vertical_alignment": "bottom" if ddec_start > 0 else "top",
                    "vertical_position": declinations[-1],
                    "vertical_shift": (1 if ddec_end > 0 else -1) * abs_vertical_shift,
                }
            else:
                radius = 0.5 * FinderChart.MINIMUM_PATH_BOX_WIDTH
                abs_vertical_position_offset = radius
                label_position_start = {
                    "horizontal_alignment": "center",
                    "horizontal_position": center_ra,
                    "vertical_alignment": "top" if ddec_start_to_end > 0 else "bottom",
                    "vertical_position": center_dec
                    + (-1 if ddec_start_to_end > 0 else 1)
                    * abs_vertical_position_offset,
                    "vertical_shift": (-1 if ddec_start_to_end > 0 else 1)
                    * abs_vertical_shift,
                }
                label_position_end = {
                    "horizontal_alignment": "center",
                    "horizontal_position": center_ra,
                    "vertical_alignment": "bottom" if ddec_start_to_end > 0 else "top",
                    "vertical_position": center_dec
                    + (1 if ddec_start_to_end > 0 else -1)
                    * abs_vertical_position_offset,
                    "vertical_shift": (1 if ddec_start_to_end > 0 else -1)
                    * abs_vertical_shift,
                }

            # add the start time label
            self.draw_label(
                label_position_start["horizontal_position"],
                label_position_start["vertical_position"]
                + label_position_start["vertical_shift"],
                start_time.strftime("%Y-%m-%d %H:%M UT"),
                size="8",
                horizontalalignment=label_position_start["horizontal_alignment"],
                verticalalignment=label_position_start["vertical_alignment"],
                color=(0, 0, 1),
            )

            # add the end time label
            self.draw_label(
                label_position_end["horizontal_position"],
                label_position_end["vertical_position"]
                + label_position_end["vertical_shift"],
                end_time.strftime("%Y-%m-%d %H:%M UT"),
                size="8",
                horizontalalignment=label_position_end["horizontal_alignment"],
                verticalalignment=label_position_end["vertical_alignment"],
                color=(0, 0, 1),
            )

            # add a "target circle" if the movement isn't significant
            if not significant_movement:
                self.draw_circle(
                    center_ra, center_dec, FinderChart.MINIMUM_PATH_BOX_WIDTH / 2.0, "b"
                )
        else:
            # output the time range
            self.draw_label(
                center_ra,
                center_dec - 4 * u.arcmin,
                start_time.strftime("%Y-%m-%d %H:%M UT")
                + " - "
                + end_time.strftime("%Y-%m-%d %H:%M UT"),
                size="large",
                horizontalalignment="center",
                verticalalignment="bottom",
                color=(0, 0.5, 1),
            )

    def _draw_arrow_head(
        self, ra: Quantity, dec: Quantity, dra: Quantity, ddec: Quantity
    ) -> None:
        """
        Draw an arrow head.

        The arrow head consists of two lines. The arrow head's tip (i.e. the point
        common to both lines) has the given right ascension and declination. The
        direction of the arrow (which dertermined how the arrow head is oriented) is
        defined by a difference dra in right ascension and a difference ddec in
        declination.

        Parameters
        ----------
        ra : Quantity
            Right ascension of the arrow tip, as an angle.
        dec : Quantity
            Declination of the arrow tip, as an angle.
        dra : Quantity
            Difference in right ascensions, as an angle.
        ddec : Quantity
            Difference in declinations, as an angle.

        """

        h = 0.002 * u.deg
        w = 0.0013 * u.deg
        ra_correction = abs(np.cos(dec))
        v_x, v_y = dra * ra_correction, ddec
        length = np.sqrt(v_x ** 2 + v_y ** 2)
        v_x, v_y = (
            v_x.to_value(u.deg) / length.to_value(u.deg),
            v_y.to_value(u.deg) / length.to_value(u.deg),
        )  # v is normalised and points in the direction of the arrow
        u_x, u_y = -v_y, v_x  # u is normalised and orthogonal to v
        dx_1 = (-h * v_x + w * u_x) / ra_correction
        dy_1 = -h * v_y + w * u_y
        dx_2 = (-h * v_x - w * u_x) / ra_correction
        dy_2 = -h * v_y - w * u_y

        coords = np.array(
            [
                [
                    (ra + dx_1).to_value(u.deg),
                    ra.to_value(u.deg),
                    (ra + dx_2).to_value(u.deg),
                ],
                [
                    (dec + dy_1).to_value(u.deg),
                    dec.to_value(u.deg),
                    (dec + dy_2).to_value(u.deg),
                ],
            ]
        )
        self.plot.show_lines([coords], color="b", linewidth=1, alpha=1)

    def _show_magnitudes(self) -> None:
        """
        Add a label with the magnitude range and bandpass.

        """

        # create the label text
        if self.magnitude_range is not None:
            mag_min = self.magnitude_range.min_magnitude
            mag_max = self.magnitude_range.max_magnitude
            bandpass = self.magnitude_range.bandpass
            if mag_max - mag_min < 0.1:
                mag_text = bandpass + " = %.1f" % mag_max
            else:
                mag_text = bandpass + " = %.1f - %.1f" % (mag_min, mag_max)
        else:
            mag_text = "no magnitude available"

        # add the label
        self.draw_label(
            self.ra,
            self.dec - 4.8 * u.arcmin,
            mag_text,
            style="italic",
            weight="bold",
            size="large",
            horizontalalignment="center",
            verticalalignment="bottom",
            color=(0, 0.5, 1),
        )

    @staticmethod
    def _center_position(ephemerides: List[Ephemeris]) -> Tuple[Quantity, Quantity]:
        """
        Get the centre right ascension and declination for a list of ephemeris values.

        Parameters
        ----------
        ephemerides

        Returns
        -------

        """
        # find the RA, dec center
        center_ra, center_dec = EphemerisService.center_position(ephemerides)

        return center_ra, center_dec


def finder_charts(
    mode_details: ModeDetails,
    ephemeris_service: EphemerisService,
    image_service: ImageService,
    output: Callable[[FinderChart, Metadata], BinaryIO],
    title: str = None,
    start_time: datetime = None,
    end_time: datetime = None,
    basic_annotations: bool = False,
) -> Generator[BinaryIO, None, None]:
    FOV_RADIUS = 3 * u.arcmin

    # get default start and end time if need be
    now = datetime.now(pytz.utc)
    if not start_time:
        start_time = julian_day_start(now)
    if not end_time:
        end_time = julian_day_end(now)

    # get all time intervals for which a finder chart must be generated
    intervals: List[Tuple[datetime, datetime]] = []
    if ephemeris_service.is_sidereal_target():
        intervals.append((start_time, end_time))
    else:
        visibility_windows = visibility.visibility_windows(
            ephemeris_service, start_time, end_time
        )
        fov_windows = visibility.fov_fitting_intervals(
            visibility_windows, ephemeris_service, FOV_RADIUS
        )
        for window in fov_windows:
            intervals.append(window)

    for (start, end) in intervals:
        finder_chart = FinderChart(
            start_time=start,
            end_time=end,
            ephemeris_service=ephemeris_service,
            image_service=image_service,
            mode_details=mode_details,
            title=title,
            basic_annotations=basic_annotations,
        )
        out = output(finder_chart, finder_chart.metadata)
        out.seek(0)
        finder_chart.discard()
        yield out
