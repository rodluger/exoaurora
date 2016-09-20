# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 16:58:44 2016

@author: dflemin3 [David P. Fleming, University of Washington]

@email: dflemin3 (at) uw (dot) edu

Typical coronagraph noise routine parameters adapted from Robinson et al 2016
and Meadows et al 2016 for future coronagraphs.
"""

# Quantum efficiency
Q = 0.9

# Size of photometric apeture in lambda / D
# theta = X * lambda / D
X = 1.5

# Coronagraph throughput
T = 0.05

# Zodaical light surface brightness in mag/arcsec^2 in the V band
MZV = 23.0

# Number of exozodis
NEZ = 1.0

# Exozodaical light surface brightness in mag/arcsec^2 in the V band
MEZV = 22.0

# Coronagraph design constrast for space-based telescopes
C = 1.0e-10 # Expected for future corornagrpah (TMT,LUVOIR,from Meadows et al 2016)
#C = 1.0e-5 # Conservative coronagraph design contrast for ground-based telescopes

# Dark count rate [/s]
DE = 1.0e-4 # Value from Meadows et al 2016

# Number of pixels spectrum spread over in horizontal for IFS
DNHPIX = 3.0

# Read noise counts per pixel
RE = 0.1

# Maximum exposure time in hours
DTMAX = 1.0

# Telescope mirror temperature (Meadows et al 2016)
TSYS = 269.0

# Telescope emissivity (Meadows et al 2016)
EMIS = 0.9