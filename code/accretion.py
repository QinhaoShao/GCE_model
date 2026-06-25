from pylab import *

class gasAcc:
    def __init__(self, twidth, ageUni):
     self.twidth = twidth #Time step, in Gyr.
     self.ageUni = ageUni #Total evolutionary time, in Gyr.

    def mw_multi_burst(self, acc_initial, tau_initial, tau_secular, acc_second, tau_second, acc_postburst, sfe_initial, tau_sfe, sfe_secular, sfe_second, sfe_postburst, time_ini, time_bur, dt):
     #Construct the multi-phase gas accretion and star formation histories of the Milky Way.
     #acc_* are the accretion rates in each phase; tau_* are the corresponding exponential decay timescales; sfe_* are the star formation efficiency coefficients.

     t = np.arange(0, self.ageUni, self.twidth)[np.newaxis, :]
     
     dim = len(acc_initial)
     self.accH = np.zeros((dim, t.shape[1])) #Gas accretion-rate history.
     self.sfeH = np.zeros((dim, t.shape[1])) #Star formation efficiency history.

     #Initial phase: the first accretion episode, lasting until time_ini.
     id0 = (t <= time_ini)
     self.accH[:, id0[0]] = (acc_initial[:, np.newaxis] * exp(-1 * (t[:, id0[0]]) / tau_initial[:, np.newaxis]))
     self.sfeH[:, id0[0]] = sfe_initial[:, np.newaxis]

     #Secular phase: quiescent evolution in which the accretion rate continues to decline exponentially from the final initial-phase value.
     id1 = (t > time_ini) & (t <= time_bur)
     
     self.accH[:, id1[0]] = self.accH[:, id0[0]][:,-1][:, np.newaxis] * exp(-1 * (t[:,id1[0]]-time_ini) / tau_secular[:, np.newaxis])

     if np.any(sfe_secular < sfe_initial):
         sfeH_id1 = sfe_initial[:, np.newaxis] * exp(-1 * (t[:, id1[0]] - time_ini) / tau_sfe[:, np.newaxis])
         mask = sfeH_id1 < sfe_secular[:, np.newaxis]
         matrix = sfe_secular[:, np.newaxis] * np.ones_like(sfeH_id1)
         sfeH_id1[mask] = matrix[mask]
     else:
         sfeH_id1 = sfe_secular[:, np.newaxis] * np.ones((dim, np.sum(id1[0])))
     self.sfeH[:,id1[0]] = sfeH_id1

     #Second phase: the second accretion / starburst episode, lasting for dt.
     id2 = (t>time_bur)&(t<=time_bur+dt)
     self.accH[:, id2[0]] = acc_second[:, np.newaxis]*exp(-1 * (t[:, id2[0]]-time_bur) / tau_second[:, np.newaxis])
     self.sfeH[:, id2[0]] = sfe_second[:, np.newaxis]

     #Postburst phase: evolution after the second accretion / starburst episode.
     id3 = (t>time_bur+dt)
     self.accH[:, id3[0]] = acc_postburst[:, np.newaxis]
     self.sfeH[:, id3[0]] = sfe_postburst[:, np.newaxis]
