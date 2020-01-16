from typing import List

from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

from asteria.exceptions import AsteriaError
from asteria.instruments import (
    Star,
    BaseInstrument,
    GapInstrument,
)


def view_sky(instr: BaseInstrument, max_stars: int = -1) -> List[Star]:
    # Query the GAIA DR2 vizier catalogue (catalogue identifier: I/345/gaia2)
    # The '+' in "+Gmag" sorts the list by the brightest first, which is useful
    # for us later when it comes to picking the best candidate stars.
    # The '_r' column gives a handy radial distance of the star from the search
    # position ra_deg, dec_deg

    # The proper motion (pmRA and pmDE) filter is needed to remove some stars
    # that may be fast moving, e.g. block_id 78169 is an interesting test!
    # - There's a 13.5 mag star with high proper motion!
    # This is not so much a problem *now*, soon after GAIA catalogue was made,
    # but is more important for future proofing this code...
    # Limits of +-50mas/yr are selected in RA and Dec.
    vquery = Vizier(
        columns=["_r", 'Source', 'RA_ICRS', 'DE_ICRS', '+Gmag'],
        column_filters={
            "Gmag": (f"<={instr.faint_limit + 2:.2f}"),
            "pmRA": ("> -50 && < 50"),
            "pmDE": ("> -50 && < 50")
        },
        # column_filters={"Gmag":(">=13.0 && <=17.0")},
        row_limit=max_stars)
    results = vquery.query_region(instr.target,
                                  radius=instr.instr_fov,
                                  catalog="I/345/gaia2")[0]

    # if we find no stars in the field, raise an error
    if len(results) == 0:
        raise AsteriaError("No suitable guide stars in field")
    else:
        stars = []
        for row in results:
            stars.append(
                Star(
                    ra=row['RA_ICRS'] * u.deg,
                    dec=row['DE_ICRS'] * u.deg,
                    radius=row['_r'] * u.arcmin,
                    g_mag=row['Gmag'],
                ))

    return stars


def find_best_stars(instr: BaseInstrument, stars: List[Star]) -> List[Star]:
    '''
    Return the best stars as defined by the instrument's best_stars function
    '''
    stars = instr.filter(stars)

    return instr.best_stars(stars)


def produce_ds9_region(instr: BaseInstrument, stars: List[Star],
                       best_stars: List[Star]):
    # create a ds9 region to help visualise things - for testing purposes...
    # note this must be written to disk as a region file that must be loaded
    # from file by ds9 later with a system command

    # TODO: add hooks to BaseInstrument class so that each instrument can
    # define its own debug output

    def draw_circle(r, label, c, x, y, width=3, params=""):
        rline = (f"circle {x:.6f}d {y:.6f}d {r}\" "
                 f"# text={{{label}}} color={c} width={width} {params}\n")
        return rline

    # def draw_line(c, x1, y1, x2, y2, params=""):
    #     rline = (f"line {x1:.6f}d {y1:.6f}d {x2:.6f}d {y2:.6f}d "
    #              f"# color={c} {params}\n")
    #     return rline

    ra_deg = instr.target.ra.to_value(u.deg)
    dec_deg = instr.target.dec.to_value(u.deg)

    output = ""

    output += "global color=white\n"
    output += "fk5\n"
    # setup some markers for guidelines...
    # these are specific to SALT
    # 4.0 arcmin
    output += draw_circle(240,
                          "4 arcmin",
                          'green',
                          ra_deg,
                          dec_deg,
                          width=1,
                          params="select=0")
    # 5.0 arcmin
    output += draw_circle(300,
                          "5 arcmin",
                          'green',
                          ra_deg,
                          dec_deg,
                          width=1,
                          params="select=0")

    # draw inner exclusion distance
    output += draw_circle(
        instr.inner_excl_distance.to_value(u.arcsec),
        f"{instr.inner_excl_distance.to_value(u.arcmin)} arcmin",
        'red',
        ra_deg,
        dec_deg,
        width=2)

    # the actual target/centre of search; in pink so it stands out!
    output += draw_circle(10, "target", 'magenta', ra_deg, dec_deg)

    # a line near the target to make sure the slit gap width is OK
    if isinstance(instr, GapInstrument):
        gap_width = instr.slit_gap_radius.to_value(u.arcmin) * 2
        pa = instr.slit_gap_angle.to_value(u.rad)
        output += (f"box {ra_deg} {dec_deg} {gap_width}' 10' {pa}r #"
                   f" fill=0 color=red\n")

    # TODO: make these colours editable too
    star_merit_colours = {
        0: "red",  # out of range
        1: "orange",  # nearby stars
        2: "yellow",  # too dim
        3: "yellow",  # too bright
        4: "green",  # good candidate
    }

    for s in stars:
        # yellow = good candidates
        colour = "yellow"

        # red == excluded candidate guide star
        if (s.merit == 0):
            colour = "red"

        colour = star_merit_colours[s.merit]

        # write the star out to the region file, with the colour reflecting
        # its suitability
        output += "circle %.6fd %.6fd 2\" # text={%.2f} color=%s\n" % (
            s.ra.to_value(u.deg), s.dec.to_value(u.deg), s.g_mag, colour)

    # indicate the best stars
    for s in best_stars:
        if s is not None:
            output += draw_circle(
                10,
                "best star",
                'cyan',
                s.ra.to_value(u.deg),
                s.dec.to_value(u.deg),
                width=2,
            )

    return output
