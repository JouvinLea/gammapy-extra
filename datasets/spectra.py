from astropy.coordinates import SkyCoord, Angle
from gammapy.datasets import gammapy_extra
from gammapy.image import ExclusionMask
from gammapy.obs import DataStore
from gammapy.region import SkyCircleRegion
from gammapy.spectrum import SpectrumAnalysis
from gammapy.utils.energy import EnergyBounds

center = SkyCoord(83.63, 22.01, unit='deg', frame='icrs')
radius = Angle('0.3 deg')
on_region = SkyCircleRegion(pos=center, radius=radius)

bkg_method = dict(type='reflected')

#J ai defini la variable GAMMAPY_EXTRA dans mon .bashrc qui pointe sur le repertoire gammapy-extra
#exclusion_file = gammapy_extra.filename("test_datasets/spectrum/"
#                                        "dummy_exclusion.fits")
exclusion_file ="/Users/jouvin/Desktop/these/test_Gammapy/gammapy-extra/datasets/exclusion_masks/tevcat_exclusion.fits"
excl = ExclusionMask.from_fits(exclusion_file)

bounds = EnergyBounds.equal_log_spacing(1, 10, 40, unit='TeV')

#store = gammapy_extra.filename("datasets/hess-crab4")
#ds = DataStore.from_dir(store)
obs = [23523, 23559, 23592, 23526]
ds = DataStore.from_dir('hess-crab4')
ana = SpectrumAnalysis(datastore=ds, obs=obs, on_region=on_region,
                       bkg_method=bkg_method, exclusion=excl, ebounds=bounds)

ana.write_ogip_data(outdir='ogip_data')


from gammapy.datasets import gammapy_extra
from gammapy.spectrum.spectrum_analysis import SpectralFit

pha23592 = gammapy_extra.filename("datasets/ogip_data/pha_run23592.pha")
pha23526 = gammapy_extra.filename("datasets/ogip_data/pha_run23526.pha")
pha23523 = gammapy_extra.filename("datasets/ogip_data/pha_run23523.pha")
pha23559 = gammapy_extra.filename("datasets/ogip_data/pha_run23559.pha")

pha_list = [pha23592, pha23526, pha23523, pha23559]
fit = SpectralFit(pha_list)
fit.model = 'PL'
fit.energy_threshold_low = '100 GeV'
fit.energy_threshold_high = '10 TeV'
fit.run(method='sherpa')

pha_group = ["group.pha"]
fit1 = SpectralFit(pha_group)
fit1.model = 'PL'
fit1.energy_threshold_low = '100 GeV'
fit1.energy_threshold_high = '10 TeV'
fit1.run(method='sherpa')


