from pylab import *

def kimf(mrange,mtot,slope0,slope1,slope2):
    #Piecewise power-law Kroupa IMF. Given the total stellar mass formed, mtot (Msun), and the mass grid, mrange (Msun), return the number of stars in each mass bin.
    #slope0/slope1/slope2 apply to 0.01-0.08, 0.08-0.5, and 0.5-120 Msun, respectively; the normalization ensures that the integrated stellar mass equals mtot.
    a0 = slope0
    a1 = slope1
    a2 = slope2
    if a0==2:
     t0 = log(0.08)-log(0.01)
    if a0!=2:
     t0 = (0.08**(2-a0)-0.01**(2-a0))/(2-a0)
    if a1==2:
     t1 = log(0.5)-log(0.08)
    if a1!=2:
     t1 = (0.5**(2-a1)-0.08**(2-a1))/(2-a1)
    if a2==2:
     t2 = log(120)-log(0.5)
    if a2!=2:
     t2 = (120**(2-a2)-0.5**(2-a2))/(2-a2)
    coe = mtot/(t0+t1*0.08**(a1-a0)+t2*0.08**(a1-a0)*0.5**(a2-a1))
    dm = np.zeros(len(mrange))
    dm[0] = mrange[1]-mrange[0] 
    dm[1:] = mrange[1:]-mrange[:-1]
    num = np.zeros((len(mrange)))
    id0 = np.where((mrange<0.08)&(mrange>0.01))
    if len(id0[0])>0:
     num[id0] = coe*mrange[id0]**(-1*a0)*dm[id0] 
    id1 = np.where((mrange>=0.08)&(mrange<0.5))
    if len(id1[0])>0:
     num[id1] = coe*0.08**(a1-a0)*mrange[id1]**(-1*a1)*dm[id1]
    id2 = np.where(mrange>=0.5)
    if len(id2[0])>0:
     num[id2] = coe*0.08**(a1-a0)*0.5**(a2-a1)*mrange[id2]**(-1*a2)*dm[id2]
    return num   
