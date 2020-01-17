from astropy import units as u

from asteria import asteria


class GuiderSingle(asteria.BaseInstrument):
    def best_stars(self, stars):
        return [s for s in stars if s.merit >= 4][:1]


fif = GuiderSingle(
    "FIF",  # instr_name
    5 * u.arcmin - 20 * u.arcsec,  # instr_fov, arcminutes radius
    1.5 * u.arcmin,  # inner_excl_distance, arcminutes
    10 * u.arcsec,  # probe_fov, arcseconds
    8,  # bright limit
    14,  # faint limit
)


class GuiderDual(asteria.GapInstrument):
    def best_stars(self, stars):
        # rank each star linearly by magnitude and distance from exclusion zone
        def score_star(star):
            radius = (star.radius - self.inner_excl_distance) / (
                self.instr_fov - self.inner_excl_distance - self.nearby_limit
            )
            mag = 1 - (star.g_mag - self.bright_limit) / (
                self.faint_limit - self.bright_limit
            )
            return radius * mag

        # seperate stars into left and right of gap
        stars_left = [
            s for s in stars if s.merit >= 4 and s.instr_coord.lon.to_value(u.deg) < 0
        ]
        stars_right = [
            s for s in stars if s.merit >= 4 and s.instr_coord.lon.to_value(u.deg) > 0
        ]

        # sort in place by star score
        stars_left.sort(key=lambda s: score_star(s), reverse=True)
        stars_right.sort(key=lambda s: score_star(s), reverse=True)

        # compare distances between the top 3 stars on each list
        # to get longest distance between a star pair
        max_seperation = (0, None, None)  # seperation, s, t
        for s in stars_left[:3]:
            for t in stars_right[:3]:
                sep = s.instr_coord.separation(t.instr_coord)
                if sep > max_seperation[0]:
                    max_seperation = sep, s, t

        # choose the stars with the greatest seperation
        return (max_seperation[1], max_seperation[2])


pfgs = GuiderDual(
    "PFGS",  # instr_name
    5 * u.arcmin - 22 * u.arcsec,  # instr_fov, arcminutes radius
    1.0 * u.arcmin,  # inner_excl_distance, arcminutes
    11 * u.arcsec,  # probe_fov, arcseconds
    9,  # bright limit
    15,  # faint limit
    0.65 * u.arcmin,  # arcmin
    0 * u.deg,  # degrees
)
pfgs.inner_excl_shape = "square"
