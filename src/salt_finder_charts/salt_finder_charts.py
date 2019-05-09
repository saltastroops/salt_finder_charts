import io
import sys
import xml
import base64
import xml.dom.minidom
from astropy import units
from astropy.coordinates.angles import Angle
import astropy.io.fits as pyfits
import numpy as np
import datetime
import aplpy
import calendar
import bisect
import enum
import math
import ephem
import os
from PIL import Image
from typing import NamedTuple
import urllib.parse
import urllib.request


class Target(NamedTuple):
    """
    Details of a target.

    Properties
    ----------
    dec: str
        Declination, in degrees between -90 and 90.
    name:
        Target name.
    ra:
        Right ascension, in degrees between 0 and 360.
    """

    dec: float
    name: str
    ra: float


class Mode(enum.Enum):
    """
    Instrument mode.

    """

    HRS = "hrs"
    IMAGING = "imaging"
    MOS = "mos"
    SLOT = "slot"


class Survey(enum.Enum):
    """
    Image survey.

    """

    POSS1_BLUE = "poss1_blue"
    POSS1_RED = "poss1_red"
    POSS2UKSTU_BLUE = "poss2ukstu_blue"
    POSS2UKSTU_IR = "poss2ukstu_ir"
    POSS2UKSTU_RED = "poss2ukstu_red"
    TWO_MASS_H = "2mass-h"
    TWO_MASS_J = "2mass-j"
    TWO_MASS_K = "2mass-k"


class NonSiderealFindingCharts:
    """
    Class for generating finding charts for a Horizons target.
    """

    FOV_RADIUS = 0.05  # 3 arcminutes

    SUTH_LONGITUDE = 20.8108

    SUTH_LATITUDE = -32.3755556

    MIN_ALTITUDE_SINE = math.sin(math.radians(46.18))

    MAX_ALTITUDE_SINE = math.sin(math.radians(59.36))

    # Assume B to the box just covering the target path. Then
    # MINIMUM_PATH_BOX_WIDTH is the minimum width of B for
    # which the path should be plotted (in degrees).
    MINIMUM_PATH_BOX_WIDTH = 16.0 / 3600.0

    def __init__(self, ephemeris_generator):
        """
        Initialises the non-sidereal finding chart.

        Params:
        ephemeris_generator: ephemeris generator
        """
        self.ephemerides = ephemeris_generator.ephemerides()

    def fov_fitting_intervals(self, intervals):
        fitting_intervals = []

        for interval in intervals:
            start = interval[0]
            end = interval[1]
            times, ra, dec, ra_rate, dec_rate, mag = self.extract_ephemeris_data(
                start, end
            )

            # handle the transition from 360 to 0 degrees
            self.transform_ra(ra, True)

            # find the RA, dec center
            center_ra = self.find_center(ra)
            center_dec = self.find_center(dec)

            # are all positions located in the FOV around this center?
            all_in_fov = True
            for i in range(0, len(times)):
                if not self.is_in_fov(center_ra, center_dec, ra[i], dec[i]):
                    all_in_fov = False
                    break
            if all_in_fov:
                fitting_intervals.append(interval)
            else:
                # split the interval in half (more or less)
                if len(times) <= 2:
                    raise ValueError("The ephemerides are too sparse to cover the FOV.")
                td = times[-1] - times[0]
                dt = (
                    td.microseconds + (td.seconds + td.days * 24 * 3600) * 10 ** 6
                ) / 10 ** 6
                dt_center = int(round(0.5 * dt))
                center_time = times[0] + datetime.timedelta(seconds=dt_center)
                index = bisect.bisect(times, center_time)
                if index == len(times) - 1:
                    index -= 1
                interval1 = [times[0], times[index]]
                interval2 = [times[index], times[-1]]
                if interval1[0] != interval1[1]:
                    fitting_intervals.extend(self.fov_fitting_intervals([interval1]))
                if interval2[0] != interval2[1]:
                    fitting_intervals.extend(self.fov_fitting_intervals([interval2]))

        return fitting_intervals

    @staticmethod
    def find_limits(values):
        min_value = 1e100
        max_value = -1e100

        for v in values:
            if v is None:
                return None, None
            if v > max_value:
                max_value = v
            if v < min_value:
                min_value = v

        return min_value, max_value

    @staticmethod
    def find_center(values):
        min_value, max_value = NonSiderealFindingCharts.find_limits(values)
        return (min_value + max_value) / 2.0

    @classmethod
    def is_in_fov(cls, center_ra, center_dec, ra, dec):
        dra = (ra - center_ra) * abs(math.cos(math.radians(dec)))
        ddec = dec - center_dec
        d = math.sqrt(dra ** 2 + ddec ** 2)

        return d <= cls.FOV_RADIUS

    def extract_ephemeris_data(self, start_time, end_time):
        if start_time >= end_time:
            raise ValueError("The start time must be earlier than the end time")

        all_times = self.ephemerides["times"]
        start_index = bisect.bisect_right(all_times, start_time)
        if not start_index and all_times[0] != start_index:
            raise ValueError("The start time isn't covered by the ephemerides.")
        end_index = bisect.bisect(all_times, end_time)
        if end_index == len(all_times) and all_times[-1] != end_time:
            raise ValueError("The end time isn't covered by the ephemerides.")
        if end_index > 0 and all_times[end_index - 1] == end_time:
            end_index -= 1

        times = all_times[start_index - 1: end_index + 1]
        ra = self.ephemerides["ra"][start_index - 1: end_index + 1]
        dec = self.ephemerides["dec"][start_index - 1: end_index + 1]
        ra_rate = self.ephemerides["ra_rate"][start_index - 1: end_index + 1]
        dec_rate = self.ephemerides["dec_rate"][start_index - 1: end_index + 1]
        mag = self.ephemerides["mag"][start_index - 1: end_index + 1]

        return times, ra, dec, ra_rate, dec_rate, mag

    # Transforms the given right ascensions. If the continuousAt360 flag
    # is true, values between 0 and 1 degree are increased by 360 degrees
    # if there are values between 359 and 360 degrees. Otherwise values
    # greater equal 360 degrees are reduced by 360 degrees.
    @staticmethod
    def transform_ra(ra, continuous_at_360):
        if continuous_at_360:
            just_before_360 = False
            just_after_0 = False
            for r in ra:
                if r >= 359:
                    just_before_360 = True
                if r <= 1:
                    just_after_0 = True
            if just_before_360 and just_after_0:
                for i in range(0, len(ra)):
                    if ra[i] <= 1:
                        ra[i] += 360
        else:
            for i in range(0, len(ra)):
                if ra[i] >= 360:
                    ra[i] -= 360

    def salt_observer(self, t):
        """
        Returns a PyEphem Observer instance for the right ascension
        and declination of Sutherland.

        Parameters:
        t: datetime for the observer

        Returns:
        the PyEphem Observer instance
        """
        observer = ephem.Observer()
        observer.lat = math.radians(self.SUTH_LATITUDE)
        observer.lon = math.radians(self.SUTH_LONGITUDE)
        observer.date = ephem.Date(t)

        return observer

    def is_visible_with_salt(self, ra, dec, t):
        """
        Checks whether a target is visible by SALT.

        Parameters:
        ra: right ascension (in degrees)
        dec: declination (in degrees)
        t: datetime
        """
        observer = self.salt_observer(t)

        lst = observer.sidereal_time()
        hour_angle = lst - math.radians(
            ra
        )  # local sidereal time is stored as an angle in radians

        phi = math.radians(self.SUTH_LATITUDE)
        delta = math.radians(dec)
        h = hour_angle

        sin_altitude = math.sin(phi) * math.sin(delta) + math.cos(phi) * math.cos(
            delta
        ) * math.cos(h)

        return self.MIN_ALTITUDE_SINE < sin_altitude < self.MAX_ALTITUDE_SINE

    def visibility_windows_next_night(self, t):
        """
        Returns the visibility windows for the night following
        the given time.

        Parameters:
        t: time

        Returns:
        the visibility windows for the next night
        """

        # get the night data
        observer = self.salt_observer(t)
        sunset = observer.next_setting(ephem.Sun()).datetime()
        sunrise = observer.next_rising(ephem.Sun()).datetime()
        night_ephemerides = self.extract_ephemeris_data(sunset, sunrise)

        # get the visibility windows
        windows = []
        dt = datetime.timedelta(seconds=300)
        tp = sunset
        window_start = None
        window_end = None
        while tp <= sunrise:
            index = bisect.bisect(night_ephemerides[0], tp)
            ra = night_ephemerides[1][index]
            dec = night_ephemerides[2][index]
            visible = self.is_visible_with_salt(ra, dec, tp)
            if visible:
                if window_start is None:
                    window_start = tp
                window_end = tp
            else:
                if window_start is not None:
                    windows.append([window_start - dt, window_end + dt])
                    window_start = None
                    window_end = None
            tp += dt
        if window_start is not None:
            windows.append((window_start - dt, window_end + dt))

        return windows

    def visibility_windows(self, start_time, end_time):
        """
        Returns the visibility windows for all the nights in the
        specified time interval.

        Both the start and end time should be outside a night.

        Parameters:
        start_time: start time
        end_time: end time

        Returns:
        the visibility windows
        """
        if start_time >= end_time:
            raise Exception("The start time must be earlier than the end time")

        t = start_time
        dt = datetime.timedelta(days=1)
        windows = []
        while t < end_time:
            windows_for_night = self.visibility_windows_next_night(t)
            for window in windows_for_night:
                window_start = window[0]
                if window_start < end_time:
                    window_end = min(window[1], end_time)
                    windows.append((window_start, window_end))
            t += dt

        return windows

    def generate_non_sidereal_finder_charts(
        self,
        start_time,
        end_time,
        mode,
        pa,
        title,
        survey,
        output_format,
        output_dir,
        file_basename,
        slitwidth=None,
        basic_annotations=False,
    ):
        """
        Generates all the finder charts for the specified time interval.

        One chart is generated per visibility window, unless the target path
        won't fit into the field of view.

        Parameters:
        start_time: start time
        end_time: end time
        mode: mode
        pa: position angle
        title: finding chart title
        survey: survey
        output_format: output format
        out_dir: output directory
        file_basename: basename for the generated files
        basic_annotations: whether not to include ephemeris annotations
        """
        visibility_windows = self.visibility_windows(start_time, end_time)
        fov_windows = self.fov_fitting_intervals(visibility_windows)
        fc_paths = []
        for fov_window in fov_windows:
            window_start = fov_window[0]
            window_end = fov_window[1]
            times, ra, dec, ra_rate, dec_rate, mag = self.extract_ephemeris_data(
                window_start, window_end
            )

            # handle the transition from 360 to 0 degrees
            self.transform_ra(ra, True)

            # get the finding chart
            ephemerides = {
                "times": times,
                "ra": ra,
                "dec": dec,
                "ra_rate": ra_rate,
                "dec_rate": dec_rate,
                "mag": mag,
            }
            fc = FindingChart.non_sidereal_finding_chart(
                ephemerides=ephemerides,
                mode=mode,
                pa=pa,
                title=title,
                survey=survey,
                slitwidth=slitwidth,
                basic_annotations=basic_annotations,
            )
            fc_path = os.path.join(
                output_dir,
                file_basename
                + "_"
                + str(self.timestamp(window_start))
                + "_"
                + str(self.timestamp(window_end))
                + "."
                + output_format,
            )
            fout = open(fc_path, "w")
            fc.save(fout, output_format)
            fc_paths.append(fc_path)
            fc.discard()

        return fc_paths

    @staticmethod
    def timestamp(t):
        """
        Returns the Unix timestamp (i.e. the seconds since 1 January 1970)
        for the specified datetime.

        Params:
        t: datetime (in UT)

        Returns:
        the timestamp
        """
        return calendar.timegm(t.timetuple())


class FindingChart:
    STSCI_SURVEYS = [
        "poss2ukstu_red",
        "poss2ukstu_blue",
        "poss2ukstu_ir",
        "poss1_red",
        "poss1_blue",
    ]

    SKY_VIEW_SURVEYS = ["2mass-j", "2mass-h", "2mass-k"]

    @staticmethod
    def sidereal_finding_chart(
        mode,
        pa,
        title,
        ra=None,
        dec=None,
        survey=None,
        fits=None,
        slitwidth=None,
        mos_mask=None,
        basic_annotations=False,
    ):
        return FindingChart(
            mode=mode,
            ra=ra,
            dec=dec,
            pa=pa,
            title=title,
            survey=survey,
            fits=fits,
            slitwidth=slitwidth,
            mos_mask=mos_mask,
            basic_annotations=basic_annotations,
        )

    @staticmethod
    def non_sidereal_finding_chart(
        ephemerides, mode, pa, title, survey, slitwidth=None, basic_annotations=False
    ):
        center_ra, center_dec = FindingChart.center_position(ephemerides)
        if not pa:
            pa = FindingChart.position_angle(ephemerides)
        fc = FindingChart.sidereal_finding_chart(
            mode=mode,
            pa=pa,
            title=title,
            ra=center_ra,
            dec=center_dec,
            survey=survey,
            slitwidth=slitwidth,
            basic_annotations=basic_annotations,
        )
        fc.plot_ephem(ephemerides, basic_annotations)

        return fc

    @staticmethod
    def center_position(ephemerides):
        # find the RA, dec center
        center_ra = NonSiderealFindingCharts.find_center(ephemerides["ra"])
        center_dec = NonSiderealFindingCharts.find_center(ephemerides["dec"])

        # avoid right ascensions greater than 360 degrees
        if center_ra >= 360:
            center_ra -= 360

        return center_ra, center_dec

    @staticmethod
    def position_angle(ephemerides):
        ra_min, ra_max = NonSiderealFindingCharts.find_limits(ephemerides["ra"])
        dec_min, dec_max = NonSiderealFindingCharts.find_limits(ephemerides["dec"])

        dra = ra_max - ra_min
        ddec = dec_max - dec_min

        # The north direction is the direction of the x axis, the east direction is in
        # the direction of the y axis. So right ascensions are x values, and
        # declinations are y values. Right ascensions have to be corrected for
        # declination.
        dx = ddec
        dy = dra * math.cos(math.radians(0.5 * (dec_min + dec_max)))

        # get the position as radians between 0 and pi
        try:
            radians = math.atan(dy / dx)
        except ZeroDivisionError:
            radians = 0.5 * math.pi
        if radians < 0:
            radians += math.pi

        # return the position angle in degrees
        return math.degrees(radians)

    def __init__(
        self,
        target: Target,
        mode: Mode,
        pa: float = None,
        survey: Survey = None,
        fits=None,
        bandpass=None,
        mag_range=None,
        slitwidth=None,
        mos_mask=None,
        basic_annotations=False,
    ):
        # convert enum values to strings
        mode_value = mode.value
        survey_value = survey.value

        # get target details
        ra = target.ra
        dec = target.dec
        target_name = target.name

        # is the mode allowed?
        if mode_value != "imaging":
            raise ValueError("Currently only the IMAGING mode is supported.")
        if mode_value == "hrs":
            mode_value = "imaging"
        if mode_value not in ["imaging", "ls", "mos", "slot"]:
            raise ValueError("Unknown observation mode: " + mode_value)

        # an image server or a FITS file must be supplied
        if not survey_value and not fits:
            raise ValueError("An image server or a FITS file must be supplied")
        if survey_value and fits:
            raise ValueError("The survey and fits argument are mutually exclusive")

        # one source of ephemerides, please
        ephemeris_arguments = 0
        if ra is not None or dec is not None:
            ephemeris_arguments += 1
        if mos_mask:
            ephemeris_arguments += 1
        if ephemeris_arguments == 0:
            raise ValueError(
                "The ra and dec arguments or the mos_mask argument must be supplied"
            )
        if ephemeris_arguments > 1:
            raise ValueError(
                "The ra and dec arguments and the mos_mask argument are mutually exclusive."  # noqa: E501
            )

        # ra requires dec, and vice versa
        if ra is not None and dec is None:
            raise ValueError("The ra argument must be used with the dec argument")
        if dec is not None and ra is None:
            raise ValueError("The dec argument must be used with the ra argument")

        # the long slit mode and the slit width require each other
        if mode_value == "ls" and not slitwidth:
            slitwidth = 1.5
        if (mode_value == "ls" and not slitwidth) or (mode_value != "ls" and slitwidth):
            raise ValueError(
                "The slitwidth argument and the ls mode require each other"
            )

        # the mos mode and MOS mask require each other
        if (mode_value == "mos" and not mos_mask) or (mode_value != "mos" and mos_mask):
            raise ValueError(
                "The mos_mask argument and the mos mode require each other"
            )

        self.mode = mode_value
        self.title = target_name
        self.pa = pa
        self.survey = survey_value
        self.ra = ra
        self.dec = dec
        self.bandpass = bandpass
        self.mag_range = mag_range
        self.slitwidth = slitwidth
        self.mos_slits = None
        self.mos_refs = None
        self.plot = None
        self.basic_annotations = basic_annotations

        self.survey_figure = None

        if mode_value == "mos":
            doc = xml.dom.minidom.parseString(mos_mask)

            pars = doc.getElementsByTagName("parameter")
            self.mos_slits = doc.getElementsByTagName("slit")
            self.mos_refs = doc.getElementsByTagName("refstar")

            parameters = {}

            for par in pars:
                name = par.getAttribute("name")
                val = par.getAttribute("value")
                parameters[name] = val

            self.ra = float(parameters["CENTERRA"])
            self.dec = float(parameters["CENTERDEC"])
            self.pa = float(parameters["ROTANGLE"])
        else:
            self.ra = ra
            self.dec = dec
            self.pa = pa

        if fits:
            self.hdu = self.get_fits(fits)
        elif survey_value:
            self.hdu = self.get_dss(survey_value, self.ra, self.dec)
            self.survey_figure = aplpy.FITSFigure(self.hdu)
        else:
            raise Exception("A value must be supplied for the survey or fits argument.")

        self.init_plot()

    def discard(self):
        self.plot.close()
        if self.survey_figure:
            self.survey_figure.close()

    # save the plot
    def save(self):
        out = io.BytesIO()
        self.plot.save(out, format="png")
        img = Image.open(out)
        img.load()
        return img

    # grab MOS xml definition from WM given account and barcode
    @staticmethod
    def get_slitmask_xml(username, password, barcode):
        """
        Return the slit mask XML as a DOM document.
        """
        encoded_username = base64.encodestring(username).strip()
        encoded_password = base64.encodestring(password).strip()

        mask_url = "https://www.salt.ac.za/wm/downloads/SlitmaskXml.php"

        # We pass the parameters in a GET request.
        url = mask_url + "?username=%s&password=%s&Barcode=%s" % (
            encoded_username,
            encoded_password,
            barcode,
        )
        response = urllib.request.urlopen(url)
        dom = xml.dom.minidom.parse(response)

        # Handle the case that the request wasn't successful.
        if dom.documentElement.tagName == "Invalid":
            raise Exception("You are not allowed to view the slit mask XML.")

        return dom

    # grab 10' x 10' image from server and pull it into pyfits
    @staticmethod
    def get_dss(survey, ra, dec):
        url = None
        params = None
        if survey in FindingChart.STSCI_SURVEYS:
            url = "http://archive.stsci.edu/cgi-bin/dss_search"
            params = urllib.parse.urlencode(
                {
                    "v": survey,
                    "r": "%f" % ra,
                    "d": "%f" % dec,
                    "e": "J2000",
                    "h": 10.0,
                    "w": 10.0,
                    "f": "fits",
                    "c": "none",
                }
            ).encode("utf-8")
        elif survey in FindingChart.SKY_VIEW_SURVEYS:
            ra = Angle(ra, unit=units.degree)
            dec = Angle(dec, unit=units.degree)
            url = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
            params = urllib.parse.urlencode(
                {
                    "Position": "'%d %d %f, %d %d %f'"
                    % (
                        round(ra.hms[0]),
                        ra.hms[1],
                        ra.hms[2],
                        round(dec.dms[0]),
                        abs(dec.dms[1]),
                        abs(dec.dms[2]),
                    ),
                    "Survey": survey,
                    "Coordinates": "J2000",
                    "Return": "FITS",
                    "Pixels": 700,
                    "Size": 0.1667,
                }
            ).encode("utf-8")
        else:
            raise Exception("Unsupported survey: " + survey)
        fits_data = io.BytesIO()
        data = urllib.request.urlopen(url, params).read()
        fits_data.write(data)
        fits_data.seek(0)
        return pyfits.open(fits_data)

    # grab uploaded base64-encoded FITS
    @staticmethod
    def get_fits(b64str):
        fits_data = io.StringIO()
        fits_data.write(base64.b64decode(b64str))
        fits_data.seek(0)
        return pyfits.open(fits_data)

    # draw a line centered at ra,dec of a given length at a given angle
    # theta,ra,dec => deg; length => arcmin
    def draw_line(self, theta, length, ra, dec, color="b", linewidth=1.0, alpha=0.7):
        theta = theta * np.pi / 180.0
        length /= 2.0
        dx = np.sin(theta) * length / (np.cos(dec * np.pi / 180.0) * 60.0)
        dy = np.cos(theta) * length / 60.0
        coords = np.array([[ra + dx, ra - dx], [dec + dy, dec - dy]])
        self.plot.show_lines([coords], color=color, linewidth=linewidth, alpha=alpha)

    # draw a box centered at ra,dec of a given length and width at a given angle
    # theta,ra,dec => deg; width, height => arcmin
    def draw_box(self, theta, width, length, ra, dec, color, linewidth=1, alpha=0.7):
        theta = theta * np.pi / 180.0
        length /= 2.0
        width /= 2.0
        # position of line centers
        ra_l = ra + np.cos(theta) * width / (np.cos(dec * np.pi / 180.0) * 60.0)
        ra_r = ra - np.cos(theta) * width / (np.cos(dec * np.pi / 180.0) * 60.0)
        dec_l = dec - np.sin(theta) * width / 60.0
        dec_r = dec + np.sin(theta) * width / 60.0

        dx = np.sin(theta) * length / (np.cos(dec * np.pi / 180.0) * 60.0)
        dy = np.cos(theta) * length / 60.0
        coords = np.array(
            [
                [ra_l, ra_l + dx, ra_r + dx, ra_r - dx, ra_l - dx, ra_l],
                [dec_l, dec_l + dy, dec_r + dy, dec_r - dy, dec_l - dy, dec_l],
            ]
        )
        self.plot.show_lines([coords], color=color, linewidth=linewidth, alpha=alpha)

    # draw slits and reference boxes for MOS
    def draw_mos_mask(self, slits, refs, pa):
        # draw the slits
        for slit in slits:
            tilt = 0.0
            if "tilt" in slit.attributes.keys():
                tilt = float(slit.attributes["tilt"].value)
            self.draw_box(
                pa + tilt,
                float(slit.attributes["width"].value) / 60.0,
                float(slit.attributes["length"].value) / 60.0,
                float(slit.attributes["xce"].value),
                float(slit.attributes["yce"].value),
                color="r",
            )
        # make bigger boxes around the reference objects
        for ref in refs:
            self.draw_box(
                pa,
                5.0 / 60.0,
                5.0 / 60.0,
                float(ref.attributes["xce"].value),
                float(ref.attributes["yce"].value),
                color=(1, 1, 0),
                linewidth=2,
            )

    # plot the object positions as a function of time from the ephem file
    def plot_ephem(self, ephemerides, basic_annotations):
        center_ra, center_dec = FindingChart.center_position(ephemerides)

        # is the target moving much?
        ra_min, ra_max = NonSiderealFindingCharts.find_limits(ephemerides["ra"])
        dec_min, dec_max = NonSiderealFindingCharts.find_limits(ephemerides["dec"])
        mag_min, mag_max = NonSiderealFindingCharts.find_limits(ephemerides["mag"])
        ra_width = ra_max - ra_min
        dec_width = dec_max - dec_min
        if (
            ra_width > NonSiderealFindingCharts.MINIMUM_PATH_BOX_WIDTH
            or dec_width > NonSiderealFindingCharts.MINIMUM_PATH_BOX_WIDTH
        ):
            significant_movement = True
        else:
            significant_movement = False

        ra_pos = ephemerides["ra"]
        dec_pos = ephemerides["dec"]
        start_time = ephemerides["times"][0]
        end_time = ephemerides["times"][-1]

        dra_start_to_end = ra_pos[-1] - ra_pos[0]
        ddec_start_to_end = dec_pos[-1] - dec_pos[0]

        # plot the target's path
        if significant_movement and not basic_annotations:
            lv = np.vstack([ra_pos, dec_pos])
            self.plot.show_lines(
                [lv], layer="object_path_lines", edgecolor="blue", linestyle="solid"
            )

        # direction at the start and end
        ddec_start = dec_pos[1] - dec_pos[0]
        dra_end = ra_pos[-1] - ra_pos[-2]
        ddec_end = dec_pos[-1] - dec_pos[-2]

        if not basic_annotations:
            if significant_movement:
                # plot the arrow at the end time to show the direction
                self.show_arrow_head(ra_pos[-1], dec_pos[-1], dra_end, ddec_end)
            else:
                ra_correction = abs(math.cos(math.radians(center_dec)))
                v_x, v_y = dra_start_to_end * ra_correction, ddec_start_to_end
                length = math.sqrt(v_x ** 2 + v_y ** 2)
                v_x, v_y = v_x / length, v_y / length
                self.show_arrow_head(
                    center_ra + 0.0013 * v_x / ra_correction,
                    center_dec + 0.0013 * v_y,
                    dra_start_to_end,
                    ddec_start_to_end,
                )

            # the labels shouldn't overlap with the path
            abs_vertical_shift = 0.002
            if significant_movement:
                label_position_start = {
                    "horizontal_alignment": "center",
                    "horizontal_position": ra_pos[0],
                    "vertical_alignment": "top" if ddec_start > 0 else "bottom",
                    "vertical_position": dec_pos[0],
                    "vertical_shift": (-1 if ddec_start > 0 else 1)
                    * abs_vertical_shift,
                }
                label_position_end = {
                    "horizontal_alignment": "center",
                    "horizontal_position": ra_pos[-1],
                    "vertical_alignment": "bottom" if ddec_start > 0 else "top",
                    "vertical_position": dec_pos[-1],
                    "vertical_shift": (1 if ddec_end > 0 else -1) * abs_vertical_shift,
                }
            else:
                radius = 0.5 * NonSiderealFindingCharts.MINIMUM_PATH_BOX_WIDTH
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
            self.plot.add_label(
                label_position_start["horizontal_position"],
                label_position_start["vertical_position"]
                + label_position_start["vertical_shift"],
                start_time.strftime("%Y-%m-%d %H:%M UT"),
                relative=False,
                size="8",
                horizontalalignment=label_position_start["horizontal_alignment"],
                verticalalignment=label_position_start["vertical_alignment"],
                color=(0, 0, 1),
            )

            # add the end time label
            self.plot.add_label(
                label_position_end["horizontal_position"],
                label_position_end["vertical_position"]
                + label_position_end["vertical_shift"],
                end_time.strftime("%Y-%m-%d %H:%M UT"),
                relative=False,
                size="8",
                horizontalalignment=label_position_end["horizontal_alignment"],
                verticalalignment=label_position_end["vertical_alignment"],
                color=(0, 0, 1),
            )

            # add a "target circle" if the movement isn't significant
            if not significant_movement:
                self.plot.show_circles(
                    [center_ra],
                    [center_dec],
                    [NonSiderealFindingCharts.MINIMUM_PATH_BOX_WIDTH / 2.0],
                    edgecolor="b",
                )
        else:
            # output the time range
            self.plot.add_label(
                center_ra,
                center_dec - 4 / 60.0,
                start_time.strftime("%Y-%m-%d %H:%M UT")
                + " - "
                + end_time.strftime("%Y-%m-%d %H:%M UT"),
                relative=False,
                size="large",
                horizontalalignment="center",
                verticalalignment="bottom",
                color=(0, 0.5, 1),
            )

        # add the magnitude range
        if mag_min is not None and mag_max is not None:
            self.show_magnitudes(center_ra, center_dec, (mag_min, mag_max), "V")
        else:
            self.show_magnitudes(center_ra, center_dec, None, "V")

    def show_arrow_head(self, ra, dec, dra, ddec):
        h = 0.002
        w = 0.0013
        ra_correction = abs(math.cos(math.radians(dec)))
        v_x, v_y = dra * ra_correction, ddec
        length = math.sqrt(v_x ** 2 + v_y ** 2)
        v_x, v_y = (
            v_x / length,
            v_y / length,
        )  # v is normalised and points in the direction of the arrow
        u_x, u_y = -v_y, v_x  # u is normalised and orthogonal to v
        dx_1 = (-h * v_x + w * u_x) / ra_correction
        dy_1 = -h * v_y + w * u_y
        dx_2 = (-h * v_x - w * u_x) / ra_correction
        dy_2 = -h * v_y - w * u_y

        coords = np.array([[ra + dx_1, ra, ra + dx_2], [dec + dy_1, dec, dec + dy_2]])
        self.plot.show_lines([coords], color="b", linewidth=1, alpha=1)

    def show_magnitudes(self, ra, dec, mag_range, bandpass):
        if mag_range is not None:
            mag_min = mag_range[0]
            mag_max = mag_range[1]
            if mag_max - mag_min < 0.1:
                mag_text = bandpass + " = %.1f" % mag_max
            else:
                mag_text = bandpass + " = %.1f - %.1f" % (mag_min, mag_max)
        else:
            mag_text = "no magnitude available"
        self.plot.add_label(
            ra,
            dec - 4.8 / 60.0,
            mag_text,
            style="italic",
            weight="bold",
            size="large",
            relative=False,
            horizontalalignment="center",
            verticalalignment="bottom",
            color=(0, 0.5, 1),
        )

    # set up basic plot
    def init_plot(self):
        survey_names = {
            "poss2ukstu_red": "POSS2/UKSTU Red",
            "poss2ukstu_blue": "POSS2/UKSTU Blue",
            "poss2ukstu_ir": "POSS2/UKSTU IR",
            "poss1_blue": "POSS1 Blue",
            "poss1_red": "POSS1 Red",
            "2mass-j": "2MASS-J",
            "2mass-h": "2MASS-H",
            "2mass-k": "2MASS-K",
            "own_fits": "User-supplied",
        }

        out = sys.stdout
        sys.stdout = open("/dev/null", "w")
        self.plot = aplpy.FITSFigure(self.hdu)
        self.plot.show_grayscale()
        self.plot.set_theme("publication")
        sys.stdout = out

        self.plot.add_label(
            0.95,
            -0.05,
            "PA = %.1f" % self.pa,
            relative=True,
            style="italic",
            weight="bold",
        )

        self.plot.add_label(
            0.5,
            1.03,
            self.title,
            relative=True,
            style="italic",
            weight="bold",
            size="large",
        )
        self.plot.add_label(
            -0.05,
            -0.05,
            "%s" % survey_names[self.survey],
            relative=True,
            style="italic",
            weight="bold",
        )

        self.plot.add_grid()
        self.plot.grid.set_alpha(0.2)
        self.plot.grid.set_color("b")

        self.plot.show_circles(
            [self.ra, self.ra],
            [self.dec, self.dec],
            [4.0 / 60.0, 5.0 / 60.0],
            edgecolor="g",
        )

        self.plot.add_label(
            0.79,
            0.79,
            "RSS",
            relative=True,
            style="italic",
            weight="bold",
            size="large",
            horizontalalignment="left",
            color=(0, 0, 1),
        )
        self.plot.add_label(
            0.86,
            0.86,
            "SCAM",
            relative=True,
            style="italic",
            weight="bold",
            size="large",
            horizontalalignment="left",
            color=(0, 0, 1),
        )
        self.plot.add_label(
            self.ra,
            self.dec + 4.8 / 60.0,
            "N",
            style="italic",
            weight="bold",
            size="large",
            color=(0, 0.5, 1),
        )
        self.plot.add_label(
            self.ra + 4.8 / (np.abs(np.cos(self.dec * np.pi / 180.0)) * 60),
            self.dec,
            "E",
            style="italic",
            weight="bold",
            size="large",
            horizontalalignment="right",
            color=(0, 0.5, 1),
        )
        self.draw_line(0, 8, self.ra, self.dec, color="g", linewidth=0.5, alpha=1.0)
        self.draw_line(90, 8, self.ra, self.dec, color="g", linewidth=0.5, alpha=1.0)

        if self.bandpass and self.mag_range:
            self.show_magnitudes(self.ra, self.dec, self.mag_range, self.bandpass)

        if self.mode == "mos" and not self.basic_annotations:
            self.draw_mos_mask(self.mos_slits, self.mos_refs, self.pa)

        if self.mode == "ls" and not self.basic_annotations:
            self.draw_line(
                self.pa, 8.0, self.ra, self.dec, color="r", linewidth=3, alpha=0.5
            )

        if self.mode == "slot" and not self.basic_annotations:
            self.draw_box(
                self.pa + 90,
                2.0 / 6.0,
                10.0,
                self.ra,
                self.dec,
                color="r",
                linewidth=2,
                alpha=0.5,
            )

        if self.mode == "imaging":
            self.plot.show_circles([self.ra], [self.dec], [0.8 / 60.0], edgecolor="g")
            self.plot.add_label(
                0.57,
                0.57,
                "BVIT",
                relative=True,
                style="italic",
                weight="bold",
                size="large",
                horizontalalignment="left",
                color=(0, 0, 1),
            )


def finder_chart(
    mode: Mode, target: Target, pa: float, survey: Survey = Survey.POSS2UKSTU_RED
) -> Image:
    fc = FindingChart(
        mode=mode,
        target=target,
        pa=pa,
        survey=survey,
        bandpass=None,
        mag_range=None,
        slitwidth=None,
        mos_mask=None,
        basic_annotations=False,
    )
    return fc.save()
