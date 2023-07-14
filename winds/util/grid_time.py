"""
winds.util.grid_time
====================

"""

import numpy as np
from netCDF4 import num2date
from datetime import datetime

from pyart.config import get_field_name, get_metadata


def datetimes_from_grid(grid):
    """ Return datetimes of grid. """

    return num2date(grid.time['data'], grid.time['units'])
