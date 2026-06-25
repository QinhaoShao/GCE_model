from astropy.io import ascii
from pylab import *
from imf import kimf
from astropy.table import Table
from accretion import *
from astropy.io import fits
import pandas as pd

'''
Main class for the Milky Way chemical evolution model. Workflow:
1. Read the gas accretion history, star formation efficiency history, and model parameters from run_MW.py.
2. Load the yield tables for AGB stars, CCSNe, SNe Ia, NSMs, and MRSNe.
3. Map the stellar mass formed at each time step to future ejecta and event rates using stellar lifetimes and DTDs.
4. Evolve gas mass, elemental masses, SFR, and event rates step by step in time.
5. Write one FITS file per radial annulus: the table extension stores log10(gasmass), log10(stellarmass), log10(sfr), [Fe/H], and event rates, while the image extension stores [X/H] for 84 elements.

variables:
- output_dir: Output directory.
- model_filenames: FITS output filename for each radial annulus.
- gas_accretion_history: Gas accretion-rate history accH and star formation efficiency history sfeH.
- initial_gas_mass: Initial gas mass.
- imf_type: Initial mass function type.
- outflow_mass_loading: Reserved gas outflow parameter; outflows are not included in the current model.
- time_step_gyr: Time step, in Gyr.
- ks_index: Power-law index of the Kennicutt-Schmidt relation.
- ccsn_yield_set: Core-collapse supernova yield table.
- ia_dtd_power: Power-law slope of the SN Ia DTD.
- ia_dtd_min_delay_gyr: Minimum delay time of the SN Ia DTD, in Gyr.
- n_radial_bins: Number of radial annuli.
- birth_radii_kpc: Annulus centers / birth radii, in kpc.
- migration_sigma_matrix: Radial migration strength matrix.
- burst_infall_abundance: Metallicity [M/H] of the gas accreted during the second phase.
- burst_start_gyr: Start time of the second phase, in Gyr.
- burst_duration_gyr: Duration of the second phase, in Gyr.
- fsn1a: Normalization factor for the SN Ia DTD.
- age_uni: Total evolutionary time, in Gyr.
- t: Time grid, in Gyr.
- mrange: Mass grid used for IMF normalization.
- mass_lft: Stellar mass grid matched to the lifetime and yield tables.
- area: Area conversion factor; each annulus is effectively treated per unit kpc^2.
- mgenhance: Placeholder parameter for metallicity enhancement; not used in the current version.
- zgrid: Metallicity grid [M/H] of the yield/lifetime tables, adopting Zsun = 0.014.
- masscal: Flag controlling the final-step reconstruction of the surviving stellar mass; 0 enables the reconstruction.
- dtd_1a_max: Upper integration limit for normalizing the SN Ia DTD, in Gyr.
- dtd_nsm_min: Minimum delay time of the neutron-star merger DTD, in Gyr.
- dtd_nsm_max: Upper integration limit for normalizing the neutron-star merger DTD, in Gyr.
- nsm_per_m: NSM event-number normalization per unit stellar mass formed.
- nsm_m_ej: Total ejecta mass per neutron-star merger.
- dtd_nsm_power: Power-law slope of the neutron-star merger DTD.
- mrsn_f: Fraction of the relevant massive stars that explode as magneto-rotational SNe; 0 disables this channel.
- solar: Solar abundance table.
'''

class galCE: 

 def __init__(self, output_dir, model_filenames, gas_accretion_history, initial_gas_mass, imf_type, outflow_mass_loading, time_step_gyr, ks_index, ccsn_yield_set, ia_dtd_power, ia_dtd_min_delay_gyr, n_radial_bins, birth_radii_kpc, migration_sigma_matrix, burst_infall_abundance, burst_start_gyr, burst_duration_gyr):
  self.gasAcc = gas_accretion_history
  self.outputFile = output_dir
  self.filename = model_filenames
  self.imf = imf_type
  self.fsn1a = 0.012
  self.twidth = time_step_gyr
  self.age_uni = 13.7
  self.N_bins = n_radial_bins
  self.Rb = birth_radii_kpc
  self.t = np.arange(0,self.age_uni,time_step_gyr)
  self.nks = ks_index
  self.mrange = np.arange(1200)*0.1+0.1
  self.mass_lft = np.append(np.append(np.arange(160)*0.025+1,np.arange(30)*0.5+5),np.arange(21)*5+20)
  self.area = 1.e6
  self.mgenhance = 0
  zgrid_yields = np.array([1.4e-5,2.e-5,5.e-5,1.e-4,3.e-4,1.e-3,2.e-3,3.e-3,6.e-3,8.e-3,1.e-2,1.4e-2,2.e-2])
  self.zgrid = log10(zgrid_yields/0.014)
  self.yields_ccsn = ccsn_yield_set
  self.mgas0 = initial_gas_mass
  self.masscal = 0
  self.dtd_1a_min = ia_dtd_min_delay_gyr
  self.dtd_1a_max = 21
  self.dtd_1a_power = ia_dtd_power
  self.dtd_nsm_min = 0.15
  self.dtd_nsm_max = 1.e6
  self.nsm_per_m = 2.e-5
  self.nsm_m_ej = 2.5e-2
  self.dtd_nsm_power = -1
  self.mrsn_f = 0
  self.solar = ascii.read('../solar_G07.txt',data_end=84)
  self.outflow = outflow_mass_loading
  self.sigma = migration_sigma_matrix
  self.acc_abun = burst_infall_abundance
  self.time_bur = burst_start_gyr
  self.dt = burst_duration_gyr


  temp_sampling = 10**np.arange(log10(self.dtd_1a_min),log10(self.dtd_1a_max),0.01)
  dt_sampling = temp_sampling[1:]-temp_sampling[:-1]
  self.norm_1a_dtd = 1./np.sum(temp_sampling[:-1]**self.dtd_1a_power*dt_sampling)
  
  temp_sampling = 10**np.arange(log10(self.dtd_nsm_min),log10(self.dtd_nsm_max),0.01)
  dt_sampling = temp_sampling[1:]-temp_sampling[:-1]
  self.norm_nsm_dtd = 1./np.sum(temp_sampling[:-1]**self.dtd_nsm_power*dt_sampling)

 def get_yields(self):
  #Load the yield tables for the different nucleosynthetic channels.
  self.agb = fits.getdata('../yields/AGB-Cristallo15-cube-galce.fits')

  if self.yields_ccsn=='l18':
   self.snii = fits.getdata('../yields/CCSN-LC18-cube-R-ave-galce.fits')
   self.snii[11, :, :] *= 10**(0.3) #Apply an empirical adjustment to the Mg yields from CCSNe.

  self.sn1a = ascii.read('../yields/yields_sn1a_i99_list.txt')

  self.nsm = ascii.read('../yields/yields_nsm_arnould07.txt')

  self.mrsn = ascii.read('../yields/yields_mrsn_nishimura15.txt')

 
 def streamline_yields(self,dmass,t,m_h):
  #Distribute the stellar mass formed in the current time step, dmass, into future ejecta and event rates.
  
  lft = fits.getdata('../lifetime-cube-230226.fits')/1.e9

  #For each radial annulus, select the nearest stellar-lifetime metallicity grid point according to the current gas [M/H].
  lft_reshape = lft.reshape(1,211,13)
  lft_di = np.tile(lft_reshape,(self.N_bins,1,1))
  zgrid_di = np.tile(self.zgrid,(self.N_bins,1))
  idz_grid_di = np.argmin(abs(zgrid_di-m_h[:, np.newaxis]), axis=1,)
  lft_z_di = lft_di[np.arange(lft_di.shape[0]), :, idz_grid_di]
  
  lft_z_di[lft_z_di > self.age_uni - t] = self.age_uni - t
  
  if self.imf=='kr':
   slope1 = 1.3
   slope2 = 2.3

   #Kroupa IMF: convert the total stellar mass formed into star counts on the adopted mass grid.
   nk = kimf(self.mrange,1,slope1,slope1,slope2)
   num_1 = nk[:10] * dmass[:, np.newaxis]
   num_120 = np.array([kimf(self.mass_lft, d, slope1, slope1, slope2) for d in dmass])

  numtot = num_1.sum(axis=1) + num_120.sum(axis=1)
  num = num_120

  id316 = (self.mrange>=3)&(self.mrange<=16)
  f316 = sum(nk[id316])


  #SNe Ia and NSMs follow delay-time distributions that map event delay times onto the model time grid.
  dtd_1a_id = np.arange(round((self.age_uni-t)/self.twidth))
  dtd_1a_sp = dtd_1a_id*self.twidth+self.twidth/2
  dtd_1a_valid = np.ones(len(dtd_1a_sp))
  dtd_1a_valid[dtd_1a_sp<=self.dtd_1a_min] = 0

  dtd_nsm_id = np.arange(round((self.age_uni-t)/self.twidth))
  dtd_nsm_sp = np.arange(0,13.7-t,0.001)
  dtd_nsm_valid = np.ones(len(dtd_nsm_sp))
  dtd_nsm_valid[dtd_nsm_sp<=self.dtd_nsm_min] = 0

  t_id = round(t/self.twidth)
  
  #SN Ia event-number weights in each future time step.
  weight_sn1a = numtot[:, np.newaxis]*self.fsn1a*f316*self.norm_1a_dtd*dtd_1a_sp**self.dtd_1a_power*self.twidth*dtd_1a_valid

  #MRSN event-number weights in each future time step.
  id_mrsn = (self.mass_lft>=13)&(self.mass_lft<=25)
  num_mrsn = num.copy()
  num_mrsn[:, ~id_mrsn] = 0
  weight_mrsn = num_mrsn*self.mrsn_f

  t_multi = int(self.twidth/0.001)

  #NSM ejecta-mass weights in each future time step.
  weight_nsm_dense = dmass[:, np.newaxis]*self.nsm_per_m*self.norm_nsm_dtd*dtd_nsm_sp**self.dtd_nsm_power*0.001*dtd_nsm_valid
  weight_nsm_dense[np.isnan(weight_nsm_dense)] = 0
  weight_nsm = np.zeros((self.N_bins,round(self.age_uni/self.twidth)-t_id))
  for m in range(round(self.age_uni/self.twidth)-t_id):
      weight_nsm[:,m] = np.sum(weight_nsm_dense[:,t_multi*m:t_multi*(m+1)], axis=1)*self.nsm_m_ej

  agbej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  sniiej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  sn1aej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  nsmej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  mrsnej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  
  #Mass-grid weights for the AGB event rate in each future time step.
  agb_mass = np.zeros((self.N_bins,len(self.mass_lft)))
  agb_mass[:,self.agb[0,:,0]>0] = 1
  agbrate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  weights_agbrate = agb_mass*num
  
  #Mass-grid weights for the CCSN/SN II event rate in each future time step.
  snii_mass = np.zeros((self.N_bins,len(self.mass_lft)))
  snii_mass[:,self.snii[0,:,0]>0] = 1
  sniirate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  weights_sniirate = snii_mass*num

  sn1arate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))

  nsmrate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))

  mrsnrate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  
  for di in range(self.N_bins):
   idz_grid = np.argmin(abs(self.zgrid-m_h[di]))
   lft_z = lft[:,idz_grid]
   self.lft_z = lft_z

   lft_tid = (lft_z[lft_z<self.age_uni-t]/self.twidth).astype(int)

   #For each element, accumulate the future gas return from AGB stars, CCSNe, SNe Ia, NSMs, and MRSNe.
   for el in range(84):
    weights_agb = self.agb[el,:,idz_grid]*num
   
    weights_snii = self.snii[el,:,idz_grid]*num*(1-self.mrsn_f)

    agbej[di,t_id:t_id+np.max(lft_tid)+1,el] = np.bincount(lft_tid, weights=weights_agb[di,lft_z<self.age_uni-t])
    sniiej[di,t_id:t_id+np.max(lft_tid)+1,el] = np.bincount(lft_tid, weights=weights_snii[di,lft_z<self.age_uni-t])
    sn1aej[di,t_id:,el] = np.bincount(dtd_1a_id, weights=self.sn1a['yields'][el]*weight_sn1a[di,:])
    nsmej[di,t_id:,el] = np.bincount(dtd_nsm_id, weights=self.nsm['yields'][el]*weight_nsm[di,:])
    mrsnej[di,t_id:t_id+np.max(lft_tid)+1,el] = np.bincount(lft_tid, weights=self.mrsn['yields'][el]*weight_mrsn[di,lft_z<self.age_uni-t])
    
    agbrate[di,t_id:t_id+np.max(lft_tid)+1] = np.bincount(lft_tid, weights=weights_agbrate[di,lft_z<self.age_uni-t])
    sniirate[di,t_id:t_id+np.max(lft_tid)+1] = np.bincount(lft_tid, weights=weights_sniirate[di,lft_z<self.age_uni-t])
    sn1arate[di,t_id:] = np.bincount(dtd_1a_id, weights=weight_sn1a[di,:])
    nsmrate[di,t_id:] = np.bincount(dtd_nsm_id, weights=weight_nsm[di,:]/self.nsm_m_ej)
    mrsnrate[di,t_id:t_id+np.max(lft_tid)+1] = np.bincount(lft_tid, weights=weight_mrsn[di,lft_z<self.age_uni-t])

  
  self.agbej = self.agbej + agbej
  self.sniiej = self.sniiej + sniiej
  self.sn1aej = self.sn1aej + sn1aej
  self.nsmej = self.nsmej + nsmej
  self.mrsnej = self.mrsnej + mrsnej
  
  self.agbrate = self.agbrate + agbrate
  self.sniirate = self.sniirate + sniirate
  self.sn1arate = self.sn1arate + sn1arate
  self.nsmrate = self.nsmrate + nsmrate
  self.mrsnrate = self.mrsnrate + mrsnrate

 def run(self):
  self.mgas = self.mgas0 #Current gas mass.
  self.fe_h = np.array(np.full((self.N_bins), -10)) #Current [Fe/H] in each radial annulus, initialized to a very low abundance.
  self.mass_el = 0 #Current total mass of each element in the gas.
  self.x_h = np.zeros(84)-10 #Current gas-phase number abundance X/H, initialized to a very low abundance.

  atomic_mass = ascii.read('../elements_mass.csv')
  atomic_mass_array = np.array(atomic_mass['AtomicMass'][:84])

  #Define the output fields: t is in Gyr; gasmass, stellarmass, and sfr are stored as log10 values; [Fe/H] and event rates are stored in linear units.
  output = np.zeros((self.N_bins,len(self.t)-1),dtype=[('t',float),('gasmass',float),('stellarmass',float),('sfr',float),('fe_h',float),('agb_rate',float),('snii_rate',float),('sn1a_rate',float),('nsm_rate',float),('mrsn_rate',float)]) 

  #output_abun stores [X/H] for all radial bins, time steps, and 84 elements.
  output_abun = np.zeros((self.N_bins,len(self.t)-1,84))

  
  self.x_h = np.zeros(84)-10
  self.get_yields()

  self.agbej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  self.sniiej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  self.sn1aej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  self.nsmej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  self.mrsnej = np.zeros((self.N_bins,round(self.age_uni/self.twidth),84))
  self.agbrate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  self.sniirate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  self.sn1arate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  self.nsmrate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  self.mrsnrate = np.zeros((self.N_bins,round(self.age_uni/self.twidth)))
  
  
  
  for i in range(len(self.t)-1):
   #Kennicutt-Schmidt law: compute the SFR from the gas surface density and SFE.
   sur_gas = self.mgas/self.area
   sfr = self.gasAcc.sfeH[:,i]*2.5*1.e-4*0.75**self.nks*sur_gas**self.nks*self.area*1.e-6 #KS law, Kennicutt et al. 1998
   if np.any(sfr*(self.t[i+1]-self.t[i])*1.e9>self.mgas):
    sfr = self.mgas/((self.t[i+1]-self.t[i])*1.e9)
   dmass_star = sfr*(self.t[i+1]-self.t[i])*1.e9

   #avgz estimates the stellar mass that remains alive; here it is reconstructed only near the final model step and written to stellarmass.
   avgz = np.zeros((self.N_bins,i+1,4))

   if (self.masscal==0)&(i>=len(self.t)-2):
    for j in range(i):
     idy = self.lft_z <= ((self.t[i]-self.t[j])*1.e9)
     numz = np.array([kimf(self.mass_lft,10**output['sfr'][x,j]*(self.t[j+1]-self.t[j])*1.e9,1.3,1.3,2.3) for x in range(self.N_bins)]) #Kroupa IMF.
     myoung = np.sum(numz[:, idy] * self.mass_lft[idy], axis=1)
     avgz[:,j,0] = 10**output['sfr'][:,j]*(self.t[j+1]-self.t[j])*1.e9-myoung
    avgz[:,0,1] = 0
    avgz[:,0,2] = 0

   self.streamline_yields(dmass_star, self.t[i], self.fe_h)

   output['t'][:,i] = self.t[i]
   output['gasmass'][:,i] = np.round(log10(self.mgas),4)
   idinf = ~np.isinf(avgz[:,:,0])
   output['stellarmass'][:,i] = np.round(log10(np.nansum(avgz[:,:,0] * idinf, axis=1)),4)
   output['sfr'][:,i] = np.round(np.log10(sfr),4)
   output['agb_rate'][:,i] = self.agbrate[:,i]
   output['snii_rate'][:,i] = self.sniirate[:,i]
   output['sn1a_rate'][:,i] = self.sn1arate[:,i]
   output['nsm_rate'][:,i] = self.nsmrate[:,i]
   output['mrsn_rate'][:,i] = self.mrsnrate[:,i]
   
   #Elemental mass returned to the interstellar medium by all channels in the current time step.
   rtn = self.agbej[:,i,:]+self.sniiej[:,i,:]+self.sn1aej[:,i,:]+self.nsmej[:,i,:]+self.mrsnej[:,i,:]
 
   if np.any(self.mgas>0):
    self.x_h = self.mass_el/self.mgas[:, np.newaxis]/0.75/atomic_mass_array
   
   #Gas-mass change = consumption by star formation + external accretion + stellar-evolution return.
   dgas = self.gasAcc.accH[:,i]*(self.t[i+1]-self.t[i])*1.e9
   dmass_gas = (-1)*sfr*(self.t[i+1]-self.t[i])*1.e9+dgas+np.nansum(rtn, axis=1)
   self.mgas = self.mgas+dmass_gas

   solar = np.tile(self.solar['col2'][:85], (self.N_bins, 1))
   
   #During the second accretion phase, assign a fixed [M/H] to the infalling gas; otherwise the infalling gas is assumed to be metal-free.
   if self.time_bur/self.twidth <= i <= (self.time_bur+self.dt)/self.twidth:
       dmgas_el = 10**(self.acc_abun)*(10**(solar-12))*dgas[:, np.newaxis]*0.75*atomic_mass_array
       dmgas_el = dmgas_el.data
   else:
       dmgas_el = 0
   
   self.mass_el = self.mass_el-self.mass_el/self.mgas[:, np.newaxis]*dmass_star[:, np.newaxis]+rtn+dmgas_el

   output_abun[:,i,:] = log10(self.x_h)-solar+12
   
   self.fe_h = output_abun[:,i,25]
   output['fe_h'][:,i] = self.fe_h
   

  for i in range(self.N_bins):
   hdu_1 = fits.BinTableHDU(data=output[i, :])
   hdu_2 = fits.ImageHDU(data=output_abun[i, :, :])
   hdu_list = fits.HDUList([fits.PrimaryHDU(), hdu_1, hdu_2])
   hdu_list.writeto(self.outputFile+self.filename[i], overwrite=True)
