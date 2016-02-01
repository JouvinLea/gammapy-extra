from astropy.coordinates import SkyCoord, Angle
from gammapy.datasets import gammapy_extra
from gammapy.image import ExclusionMask
from gammapy.obs import DataStore
from gammapy.region import SkyCircleRegion
from gammapy.spectrum import SpectrumAnalysis
from gammapy.utils.energy import EnergyBounds
from gammapy.spectrum.spectrum_analysis import SpectrumFit
from glob import glob

center = SkyCoord(83.63, 22.01, unit='deg', frame='icrs')
radius = Angle('0.3 deg')
on_region = SkyCircleRegion(pos=center, radius=radius)

bkg_method = dict(type='reflected', n_min=2)
exclusion_file = "/Users/jouvin/Desktop/these/test_Gammapy/gammapy-extra/datasets/exclusion_masks/tevcat_exclusion.fits"
excl = ExclusionMask.from_fits(exclusion_file)

bounds = EnergyBounds.equal_log_spacing(1, 10, 40, unit='TeV')

obs = [int(file[6:11]) for file in glob('run*.fits')]
ds = DataStore.from_dir('/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/CrabEventList/Crab')
ds1 = DataStore.from_dir('/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/CrabEventList/Crab')

ana = SpectrumAnalysis(datastore=ds, obs=obs, on_region=on_region,
                       bkg_method=bkg_method, exclusion=excl, ebounds=bounds)
anaband = SpectrumAnalysis(datastore=ds1, obs=obs, on_region=on_region,
                           bkg_method=bkg_method, exclusion=excl, ebounds=bounds)

#Band boudaries
Offmin = 0.
Offmax = 2.5
Offbin = 5.
Effmin = 0
Effmax = 100
# Effbin=15
Effbin = 5
Zenmin = 0
Zenmax = 70
# Zenbin=10
Zenbin = 5

anaband.define_spectral_groups([Offmin, Offmax], Offbin, [Effmin, Effmax], Effbin, [Zenmin, Zenmax], Zenbin)

ana.write_ogip_data(outdir='ogip_data')
anaband.write_ogip_data(outdir='groups')


print "ALL THE OBSERVATIONS"
pha_list = glob('ogip_data/pha*')
fit = SpectrumFit(pha_list)
fit.model = 'PL'
fit.energy_threshold_low = '100 GeV'
fit.energy_threshold_high = '10 TeV'
fit.run(method='sherpa')

print "GROUPING"
pha_band = glob('groups1/pha*')
fitband = SpectrumFit(pha_band)
fitband.model = 'PL'
fitband.energy_threshold_low = '100 GeV'
fitband.energy_threshold_high = '10 TeV'
fitband.run(method='sherpa')
