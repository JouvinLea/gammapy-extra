import numpy as np
from numpy.testing import assert_allclose
from numpy.testing import assert_equal
from astropy.tests.helper import assert_quantity_allclose
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS
from astropy.units import Quantity
from astropy.coordinates import Angle, SkyCoord
from astropy.modeling.models import Gaussian1D
from gammapy.utils.testing import requires_dependency, requires_data
from gammapy.datasets import gammapy_extra
from gammapy.background import GaussianBand2D, CubeBackgroundModel, EnergyOffsetBackgroundModel
from gammapy.utils.energy import EnergyBounds
from gammapy.data import ObservationTable
from gammapy.data import DataStore, EventList
from gammapy.region import SkyCircleRegion
from gammapy.background.models import compute_pie_fraction, select_events_outside_pie
from gammapy.image import make_empty_image
from gammapy.image import coordinates, bin_events_in_image, make_empty_image


def make_excluded_sources():
    centers = SkyCoord([1, 0], [2, 1], unit='deg')
    radius = Angle('0.3 deg')
    sources = SkyCircleRegion(pos=centers, radius=radius)
    catalog = Table()
    catalog["RA"] = sources.pos.data.lon
    catalog["DEC"] = sources.pos.data.lat
    catalog["Radius"] = sources.radius
    return catalog


def test_select_events_outside_pie():
    """
    Create an empty image centered on the pointing position and all the radec position of the pixels will
    define one event. Thus we create a false EventList with these pixels. We apply the select_events_outside_pie()
    and we fill the image only with the events (pixels) outside the pie.
    """
    excluded_sources = make_excluded_sources()
    excluded_sources2 = make_excluded_sources()
    pointing_position = SkyCoord(0.5, 0.5, unit='deg')
    # Define an empty image centered on the pointing position
    image = make_empty_image(nxpix=1000, nypix=1000, binsz=0.01, xref=pointing_position.ra.deg,
                             yref=pointing_position.dec.deg,
                             proj='TAN', coordsys='CEL')
    ra, dec = coordinates(image)
    events = EventList()
    # Faked EventList with the radec of all the pixel in the empty image
    events["RA"] = ra.flat
    events["DEC"] = dec.flat

    # Test that if the sources are out of the fov, it gives the index for all the events since no event will be removed
    idx = select_events_outside_pie(excluded_sources, events, pointing_position, Angle(0.3, "deg"))
    assert_allclose(np.arange(len(events)), idx)

    # Test if after calling the select_events_outside_pie, the image is 0 inside the pie and 1 outside the pie
    idx = select_events_outside_pie(excluded_sources, events, pointing_position, Angle(5, "deg"))
    events_outside = events[idx]
    image = bin_events_in_image(events_outside, image)
    image.writeto("test_pie.fits", clobber=True)


if __name__ == '__main__':
    test_select_events_outside_pie()
