import enum
from abc import ABC
from typing import Any, Optional, BinaryIO

import astropy.units as u
from astropy.units import Quantity

from salt_finder_charts.util import MOSMask


class Mode(enum.Enum):
    """
    Instrument mode.

    """

    HRS = "hrs"
    IMAGING = "imaging"
    LONGSLIT = "ls"
    MOS = "mos"
    SLOT = "slot"


class ModeDetails(ABC):
    """
    Observation mode details.

    Parameters
    ----------
    mode : Mode
        Observation mode.

    """

    def __init__(self, mode: Mode):
        self._mode = mode

    @property
    def mode(self) -> Mode:
        return self._mode

    def position_angle(self) -> Quantity:
        """
        The position angle.

        Returns
        -------
        Quantity
            The position angle.

        """

        raise NotImplementedError

    def annotate_finder_chart(self, finder_chart: Any) -> None:
        """
        Add the mode specific content to a finder chart.

        Parameters
        ----------
        finder_chart : FindingChart
            Finding chart.

        """

        raise NotImplementedError


class ImagingModeDetails(ModeDetails):
    """
    Details for the imaging mode.

    Parameters
    ----------
    pa : Optional[Quantity]
        Position angle.

    """

    def __init__(self, pa: Optional[Quantity]):
        super().__init__(Mode.IMAGING)
        self.pa = pa if pa is not None else 0 * u.deg

    def position_angle(self) -> Quantity:
        return self.pa

    def annotate_finder_chart(self, finder_chart: Any) -> None:
        # indicate field of view for BVIT
        finder_chart.draw_circle(
            finder_chart.ra, finder_chart.dec, 0.8 * u.arcmin, color="g"
        )
        finder_chart.plot.add_label(
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


class SlotModeDetails(ModeDetails):
    """
    Details for the slot mode.

    Parameters
    ----------
    pa : Optional[Quantity]
        Position angle.

    """

    def __init__(self, pa: Optional[Quantity]):
        super().__init__(Mode.SLOT)
        self.pa = pa if pa is not None else 0 * u.deg

    def position_angle(self) -> Quantity:
        return self.pa

    def annotate_finder_chart(self, finder_chart: Any) -> None:
        finder_chart.draw_centered_rectangle(
            self.pa + 90 * u.deg,
            2.0 * u.arcmin / 6.0,
            10.0 * u.arcmin,
            finder_chart.ra,
            finder_chart.dec,
            color="r",
            linewidth=2,
            alpha=0.5,
        )


class LongslitModeDetails(ModeDetails):
    """
    Details for the longslit mode.

    Parameters
    ----------
    slitwidth : Quantity
        Slitwidth, as an angle.
    pa : Optional[Quantity]
        Position angle.

    """

    def __init__(self, slitwidth: Quantity, pa: Optional[Quantity]):
        super().__init__(Mode.LONGSLIT)
        self.pa = pa if pa is not None else 0 * u.deg
        self.slitwidth = slitwidth

    def position_angle(self) -> Quantity:
        return self.pa

    def annotate_finder_chart(self, finder_chart: Any) -> None:
        # draw the slit
        finder_chart.draw_centered_rectangle(
            self.pa,
            self.slitwidth,
            8.0 * u.arcmin,
            finder_chart.ra,
            finder_chart.dec,
            color="r",
            linewidth=1,
            alpha=0.5,
        )


class MOSModeDetails(ModeDetails):
    """
    Details for the MOS mode.

    Parameters
    ----------
    mos_mask : MOSMask
        MOS mask.

    """

    def __init__(self, mos_mask: MOSMask):
        super().__init__(Mode.MOS)
        self.mos_mask = mos_mask

    def position_angle(self) -> Quantity:
        return self.mos_mask.position_angle

    def annotate_finder_chart(self, finder_chart: Any) -> None:
        # draw the slits
        pa = self.mos_mask.position_angle
        for slit in self.mos_mask.slits:
            finder_chart.draw_centered_rectangle(
                pa + slit.tilt, slit.width, slit.height, slit.ra, slit.dec, color="r"
            )

        # make bigger boxes around the reference objects
        for ref in self.mos_mask.reference_stars:
            finder_chart.draw_centered_rectangle(
                pa,
                5.0 * u.arcsec,
                5.0 * u.arcsec,
                ref.ra,
                ref.dec,
                color=(1, 1, 0),
                linewidth=2,
            )
