from astropy.io import ascii
from pylab import *
from imf import kimf
import sys
from astropy.io import fits

def mock(file_path, filename, num, delta_y, delta_o, age_idv):
    #Read the chemical evolution output for a single radial annulus: extension 1 stores the time-series table, and extension 2 stores the abundance matrix for 84 elements.
    file = fits.open(file_path + filename)
    model = file[1].data
    abun = file[2].data
    feh = abun[:,25]
    mgh = abun[:,11]
    sih = abun[:,13]
    mgfe = mgh-feh
    sife = sih-feh
    dt = model['t'][1:]-model['t'][:-1]

    mass_lft = np.append(np.append(np.arange(180)*0.025+0.5,np.arange(30)*0.5+5),np.arange(21)*5+20)

    #Mass range of stars to be sampled; this range affects the mass distribution in the mock stellar catalog.
    mvalid = mass_lft[(mass_lft>=0.5)&(mass_lft<40)] 
    lft = fits.getdata('../lifetime-cube-230313-lowm.fits')/1.e9
    zgrid_yields = np.array([1.4e-5,2.e-5,5.e-5,1.e-4,3.e-4,1.e-3,2.e-3,3.e-3,6.e-3,8.e-3,1.e-2,1.4e-2,2.e-2])
    zgrid = log10(zgrid_yields/0.014)

    pdf_ori = np.zeros((len(mvalid),len(model)))
    for i in range(len(model)-1):
        numz = kimf(mvalid,10**model['sfr'][i]*dt[i]*1.e9,1.3,1.3,2.3)
        idz_lft = np.argmin(abs(feh[i]-zgrid))
        mlife = lft[:,idz_lft][(mass_lft>=0.5)&(mass_lft<40)]
        idy = (mlife >= (13.7-model['t'][i]))
        if (len(mlife[idy])>0):
            pdf_ori[idy,i] = numz[idy]
    pdf_ori = pdf_ori/np.sum(pdf_ori)

    pdf_age = np.sum(pdf_ori,axis=0)

    pdf_age_norm = pdf_age/np.sum(pdf_age)


    mocksz = num #Number of stars in the mock stellar catalog.

    #Mock catalog: age is the age after adding observational uncertainty, age_0 is the intrinsic age, and the abundance columns include observational uncertainties.
    mock = np.zeros((mocksz,),dtype=[('age',float),('mass',float),('fe/h',float),('mg/fe',float),('si/fe',float),('age_0',float)])

    randices_age = np.random.choice(np.arange(len(model['t'])),mocksz,replace=True,p=pdf_age_norm)

    #Age and abundance observational uncertainties.
    sigma_age_y = delta_y #Age uncertainty for young stars.
    sigma_age_o = delta_o #Age uncertainty for old stars.
    sigma_feh = 0.07
    sigma_mgfe = 0.03
    sigma_sife = 0.03

    #Mock stellar ages: first sample from the age probability distribution, then add log-normal age uncertainties.
    ages = 13.7-model['t'][randices_age]
    mock['age_0'] = ages

    age_o = ages > age_idv
    age_y = ages <= age_idv
    mock['age'][age_y] = ages[age_y]*10**np.random.normal(0,sigma_age_y,len(mock[age_y]))
    mock['age'][age_o] = ages[age_o]*10**np.random.normal(0,sigma_age_o,len(mock[age_o]))


    idiv = np.where(mock['age']>13.7)
    for i in range(len(idiv[0])):
        #If the age uncertainty moves a star beyond the age of the Universe, resample it within the allowed age range.
        oldp = exp(-0.5*(log10(13.7-model['t'])-log10(ages[idiv[0]][i]))**2/sigma_age_o**2)
        oldp = oldp/sum(oldp)
        newage = model['t'][np.random.choice(np.arange(len(model['t'])),1,replace=True,p=oldp)]
        mock['age'][idiv[0][i]] = 13.7-newage[0]

    mock['fe/h'][:] = feh[randices_age]+np.random.normal(0,sigma_feh,len(mock))
    mock['mg/fe'][:] = mgfe[randices_age]+np.random.normal(0,sigma_mgfe,len(mock))
    mock['si/fe'][:] = sife[randices_age]+np.random.normal(0,sigma_sife,len(mock))

    #Mock stellar masses: sample from the surviving-mass probability distribution associated with the birth time.
    birth_star = 13.7-mock['age']
    for i in range(len(model)-1):
        idt = ((birth_star>=model['t'][i])&(birth_star<model['t'][i+1]))
        if len(birth_star[idt])>0:
            pdf_mass = pdf_ori[:,i]
            numz = kimf(mvalid,10**model['sfr'][i]*dt[i]*1.e9,1.3,1.3,2.3)
            idy = (mlife >= log10((13.7-model['t'][i])*1.e9))
            if (len(mlife[idy])>0):
                pdf_mass[idy] = numz[idy]
            if np.sum(pdf_mass)>0:
                pdf_mass = pdf_mass/np.sum(pdf_mass)
                mass = np.random.choice(mvalid,len(birth_star[idt]),replace=True,p=pdf_mass)
                mock['mass'][idt] = mass

    modelname = filename.split('fits')[0]
    ascii.write(mock, file_path + modelname +'txt',overwrite=1)
