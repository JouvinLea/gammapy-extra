import numpy as np
from astropy.table import Table
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
from gammapy.data import DataStore
from gammapy.data import ObservationGroupAxis
from gammapy.data import ObservationGroups
from gammapy.datasets import gammapy_extra
import pylab as pt
pt.ion()

dir=str(gammapy_extra.dir)+'/datasets/hess-crab4'
data_store = DataStore.from_dir(dir)
Observation_Table=data_store.obs_table
alt = Angle([0, 30, 60, 90], 'degree')
az = Angle([-90, 90, 270], 'degree')
ntels = np.array([3, 4])
list_obs_group_axis = [ObservationGroupAxis('ALT', alt, 'bin_edges'),
                       ObservationGroupAxis('AZ', az, 'bin_edges'),
                       ObservationGroupAxis('N_TELS', ntels, 'bin_values')]
obs_groups = ObservationGroups(list_obs_group_axis)
obs_table_grouped = obs_groups.group_observation_table(Observation_Table)
print(obs_table_grouped)
