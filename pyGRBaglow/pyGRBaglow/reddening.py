#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyGRBaglow import constants as cc
import numpy as np


def Pei92(wavelength, Av, z, Rv=-99.0, ext_law="smc", Xcut=False):
    """
    Extinction laws from Pei 1992 article

    Parameters
    ----------
    wavelength: `array` or `float`
        wavlength in angstroms

    Av: `float`
        amount of extinction in the V band

    z: `float`
        redshift

    Rv: `float`, optional, default: -99.
        selective attenuation Rv = Av / E(B-V)
        if 'd-99.' set values by default from article
        if a float is given, use this value instead

    ext_law: `str`
        type of extinction law to use.
        Choices: mw, lmc, smc

    Xcut: `boolean`, optional, default: False
         Whether to set attenuation to 0 for wavelength below 700 angstrom
         Useful when coupling with X-ray data

    Returns
    -------
    [Alambda_over_Av, Trans_dust]

    Alambda_over_Av : `array`
        atteanuation as a function of wavelength normalise by Av
        (attenuation in V band)

    Trans_dust: `array`
        transmission through dust as a function of wavelength

    """

    wvl = wavelength * 1e-4 / (1 + z)
    if ext_law.lower() == "smc":
        if Rv == -99.0:
            Rv = 2.93
        a = [185, 27, 0.005, 0.010, 0.012, 0.03]
        wvl_ = [0.042, 0.08, 0.22, 9.7, 18, 25]
        b = [90, 5.50, -1.95, -1.95, -1.80, 0.0]
        n = [2.0, 4.0, 2.0, 2.0, 2.0, 2.0]

    elif ext_law.lower() == "lmc":
        if Rv == -99.0:
            Rv = 3.16
        a = [175, 19, 0.023, 0.005, 0.006, 0.02]
        wvl_ = [0.046, 0.08, 0.22, 9.7, 18, 25]
        b = [90, 5.5, -1.95, -1.95, -1.8, 0.0]
        n = [2.0, 4.5, 2.0, 2.0, 2.0, 2.0]

    elif ext_law.lower() == "mw":
        if Rv == -99.0:
            Rv = 3.08
        a = [165, 14, 0.045, 0.002, 0.002, 0.012]
        wvl_ = [0.046, 0.08, 0.22, 9.7, 18, 25]
        b = [90, 4.0, -1.95, -1.95, -1.8, 0.0]
        n = [2.0, 6.5, 2.0, 2.0, 2.0, 2.0]

    sums = np.zeros(len(wvl))
    for i in range(len(a)):
        sums += a[i] / ((wvl / wvl_[i]) ** n[i] + (wvl_[i] / wvl) ** n[i] + b[i])

    # Need to check whether extrapolation is needed
    # outside the range defined in Pei92
    # convert Alambda_over_Ab to Alambda_over_Av
    Alambda_over_Av = (1.0 / Rv + 1.0) * sums

    # Applied a cut for wavelength below 700 angstrom
    # Useful when coupling with Xray data
    if Xcut:
        w = np.where(wvl < 0.07)
        Alambda_over_Av[w] = 0

    # Return optical depth due to dust reddening in funtion of wavelength
    Tau_dust = Av * Alambda_over_Av / 1.086

    Trans_dust = np.exp(-Tau_dust)

    Trans_dust[Trans_dust < 0] = 0
    Trans_dust[Trans_dust > 1] = 1

    return [Alambda_over_Av, Trans_dust]


def sne1(wavelength, Av, z, Xcut=False):
    """
    Extinction law for SNe

    Parameters
    ----------
    wavelength: `array` or `float`
        wavlength in angstroms

    Av: `float`
        amount of extinction in the V band

    z: `float`
        redshift

    Xcut: `boolean`, optional, default: False
         Whether to set attenuation to 0 for wavelength below 700 angstrom
         Useful when coupling with X-ray data

    Returns
    -------
    [Alambda_over_Av, Trans_dust]

    Alambda_over_Av : `array`
        atteanuation as a function of wavelength normalise by Av
        (attenuation in V band)

    Trans_dust: `array`
        transmission through dust as a function of wavelength

    """

    wvl = wavelength * 1e-4 / (1 + z)

    mask = wvl < 1000 * 1e-4
    Alambda_over_Av = np.zeros(len(wvl))
    # Fit available above 1000 Angstroms
    Alambda_over_Av[~mask] = (
        -2.2113e-05 / (wvl[~mask]) ** 8
        + 9.7507e-04 / (wvl[~mask]) ** 7
        - 1.7447e-02 / (wvl[~mask]) ** 6
        + 1.6186e-01 / (wvl[~mask]) ** 5
        - 8.2474e-01 / (wvl[~mask]) ** 4
        + 2.262 / (wvl[~mask]) ** 3
        - 3.13 / (wvl[~mask]) ** 2
        + 2.591 / (wvl[~mask])
        - 6.5916e-01
    )
    # Below 1000 Angstroms, use a rescaled SMC extinction law
    Alambda_over_Av[mask] = (
        Pei92(wavelength[mask], Av, z, ext_law="smc", Xcut=Xcut)[0] * 0.75794464
    )

    # Applied a cut for wavelength below 700 angstrom
    # Useful when coupling with Xray data
    if Xcut:
        w = np.where(wvl < 0.07)
        Alambda_over_Av[w] = 0

    # Return optical depth due to dust reddening in funtion of wavelength
    #Tau_dust = Av * Alambda_over_Av / 1.086
    #Trans_dust = np.exp(-Tau_dust)

    #Trans_dust[Trans_dust < 0] = 0
    #Trans_dust[Trans_dust > 1] = 1
    Trans_dust = 10**(-0.4 * Av * Alambda_over_Av)
    return [Alambda_over_Av, Trans_dust]
    
def sne(wavelength, Av, z):
    """
    Extinction law for SNe

    Parameters
    ----------
    wavelength: `array` or `float`
        wavlength in angstroms

    Av: `float`
        amount of extinction in the V band

    z: `float`
        redshift

    Returns
    -------
    [Alambda_over_Av, Trans_dust]

    Alambda_over_Av : `array`
        extinction as a function of wavelength normalise by Av
        (extinction in V band)
    Trans_dust: `array`
        transmission through dust as a function of wavelength
    """

    wvl = wavelength * 1e-4 / (1 + z)

    mask = wvl < 1000 * 1e-4
    Alambda_over_Av = np.zeros(len(wvl))
    # Fit available above 1000 Angstroms
    Alambda_over_Av[~mask] = (
        -2.2113e-05 / (wvl[~mask]) ** 8
        + 9.7507e-04 / (wvl[~mask]) ** 7
        - 1.7447e-02 / (wvl[~mask]) ** 6
        + 1.6186e-01 / (wvl[~mask]) ** 5
        - 8.2474e-01 / (wvl[~mask]) ** 4
        + 2.262 / (wvl[~mask]) ** 3
        - 3.13 / (wvl[~mask]) ** 2
        + 2.591 / (wvl[~mask])
        - 6.5916e-01
    )
    # Below 1000 Angstroms, physically incorrect but useful for the addition of X-ray data
    Alambda_over_Av[Alambda_over_Av < 0] = 0

    # Return optical depth due to dust reddening in funtion of wavelength
    Trans_dust = 10**(-0.4 * Av * Alambda_over_Av)

    return [Alambda_over_Av, Trans_dust]


def calzetti_1998 (wvl, z):
    """ Calzetti starbust reddening curve
    Wavelength: in angstrom
    z: redshift
    ebv 
    """
    wv = (wvl*1e-4)*(1+z)
    rv = 4.05
    mask1 = (wv >= 0.63) & (wv <= 1)
    mask2 = (wv >= 0.12) & (wv <= 0.63)
    k = np.zeros(len(wv))
    k [mask1]= ((1.86-0.48/wv[mask1])/wv[mask1] - 0.1)/wv[mask1] + 1.73 
    k [mask2] = 2.656*(-2.156+1.509/wv[mask2]-0.198/wv[mask2]**2 + 0.011/wv[mask2]**3)+4.88
    return k

def calzetti_2000 (wvl, z,Av,rv = 4.05):
    """ Calzetti starbust reddening curve
    Wavelength: in angstrom
    z: redshift
    es_bv = 0.44 * E(B-V) #not need to add this one
    Rv' = 4.05
    """
    wv = (wvl*1e-4)/(1+z)
    #esbv = ebv*0.44
    mask1 = (wv >= 0.63) & (wv <= 3) #& (wv <= 2.2)
    mask2 = (wv < 0.63)  &(wv >= 0) #&(wv >= 0.12)
    k = np.zeros(len(wv))
    k [mask1]= 2.659*(-1.857+1.040/wv[mask1]) + rv
    k [mask2] = 2.659*(-2.156+1.509/wv[mask2]-0.198/(wv[mask2]**2)+0.011/(wv[mask2]**3)) + rv
    k [(k == 0)] = 1
    Alambda_over_Av = k/rv
    
    #Tau_dust = Av * Alambda_over_Av / 1.086

    #Trans_dust = np.exp(-Tau_dust)

    #Trans_dust[Trans_dust < 0] = 0
    #Trans_dust[Trans_dust > 1] = 1
    
    trans = 10**(-0.4 * Av * Alambda_over_Av)
    trans[trans <= 0] = 1 
    return k/rv, trans #Trans_dust
    
def spl_law1(x,z,av,alpha):
    """ POWER LAW extinction law, Savaglio 2004
    Wavelength in angstrom
    z = redshift
    av = visual extinction
    """
    zx = x/(1+z)
    alambda_av = (5500/zx)**alpha
    #return (5500**alpha) * x**-alpha
    trans = 10**(-0.4*av*alambda_av)
    trans[trans <= 0] = 1    
    return alambda_av, trans
    


def spl_law(x,z,av,alpha):
    """ POWER LAW extinction law, Savaglio 2004
    Wavelength in angstrom
    z = redshift
    av = visual extinction
    Optical depth = 1 for wavelength < 100 nm
    """
    zx = x/(1+z)
    mask = x < 1000
    alambda_av = np.zeros(len(zx))
    alambda_av[~mask] = (5500/zx[~mask])**alpha
    #alambda_av = (5500/zx)**alpha
    #return (5500**alpha) * x**-alpha
    trans = 10**(-0.4*av*alambda_av)
    return alambda_av, trans

    


def spl_drude (wvl,z,av,alpha,eb = 1,delta_wvl = 350.,wvl_0 = 2200):
    """Simple power-law extinction model + Drude profile
    Wavelength in Angstrom

    E_b = Amplitude constant of the bump strength of the Drude profile (silicate feature)
    delta_wvl = 470 AA, width of the MW feature
    wvl_0 = 2175.8 AA, MW feature
    z = Redshift (and of the bump of the milky way)
    """
    wv = wvl/(1+z)
    mask = wvl < 1000
    alambda_av = np.zeros(len(wv))
    alambda_av[~mask] = ((5500/wv[~mask])**alpha ) + \
                        (eb*(wv[~mask]*delta_wvl)**2) / \
                        ((wv[~mask]**2-(wvl_0)**2)**2+(wv[~mask]*delta_wvl)**2)
    trans = 10**(-0.4*av*alambda_av)
    return alambda_av, trans 


def diff_ext (x,av,z,ext="0"):
    wv = (1/(x*1e-4)) * (1+z)
    f = np.poly1d(np.load('/home/nrakotondrainibe/Bureau/grb_git/codes____/xspec_v1/df%s.npy'%ext))
    g = np.poly1d(np.load('/home/nrakotondrainibe/Bureau/grb_git/codes____/xspec_v1/df%s_new.npy'%ext))
    alambda_av = np.zeros(len(wv))
    mask = (wv < 5.85)
    alambda_av [mask] = f(wv[mask])
    alambda_av [~mask] = g(wv[~mask])
    #alambda_av[alambda_av < 0] = 0    
    trans = 10**(-0.4*av*alambda_av)
    return alambda_av,trans



def gas_absorption(wavelength, z, NHx=0.2):
    """
    Compute the optical depth due to gas absorption

    Parameters
    ----------
    wavelength: `array` or `float`
        wavlength in angstroms

    z: `float`
        redshift

    NHx: `float`, optional, default: 0.2
         Metal column density from soft Xrays absortpion
         (in units 1e22 cm-2), expressed in units of equivalent
         hydrogen column density assuming solar abundances.
         In Milky Way:
         NHx/Av = 1.7 to 2.2 1e21 cm-2/mag (default set to 2)

    Returns
    ------
    Trans_gas : `array`
         Transmission coefficient due to the gas absorption either
         occuring in our galaxy or within the host.

    """
    nus = cc.c_light_m_s / (wavelength * 1e-10)
    Tau_gas = np.zeros(len(nus))

    for i in range(len(nus)):
        # photon frequency (Hz) in the rest frame
        nu = nus[i] * (1 + z)
        # photon energy (keV) in the rest frame
        E_kev = nu * cc.H_planck / (1e3 * cc.e_elec)
        E_kev2 = E_kev**2.0
        E_kev3 = E_kev**3.0

        # if E_kev < 13.6e-3:    #912 A (Lyman limit)
        #     c0=0; c1=0; c2=0
        # 41nm / 410A
        if E_kev < 0.030:
            coeffs = [0, 0, 0, "H"]
        # 12.4 nm
        elif E_kev < 0.100:
            coeffs = [17.3, 608.1, -2150, "He"]
        # 4.37 nm
        elif E_kev < 0.284:
            coeffs = [34.6, 267.9, -476.1, "C"]
        elif E_kev < 0.400:
            coeffs = [78.1, 18.8, 4.3, "N"]
        elif E_kev < 0.532:
            coeffs = [71.4, 66.8, -51.4, "O"]
        elif E_kev < 0.707:
            coeffs = [95.5, 145.8, -61.1, "Fe-L"]
        elif E_kev < 0.867:
            coeffs = [308.9, -380.6, 294.0, "Ne"]
        elif E_kev < 1.303:
            coeffs = [120.6, 169.3, -47.7, "Mg"]
        elif E_kev < 1.840:
            coeffs = [141.3, 146.8, -31.5, "Si"]
        elif E_kev < 2.471:
            coeffs = [202.7, 104.7, -17.0, "S"]
        elif E_kev < 3.210:
            coeffs = [342.7, 18.7, 0.0, "Ar"]
        elif E_kev < 4.038:
            coeffs = [352.2, 18.7, 0.0, "Ca"]
        elif E_kev < 7.111:
            coeffs = [433.9, -2.4, 0.75, "Fe"]
        elif E_kev < 8.331:
            coeffs = [629.0, 30.9, 0.0, "Ni"]
        # 124pm/1.24A
        elif E_kev < 10.0:
            coeffs = [701.2, 25.2, 0.0, "..."]
        else:
            coeffs = [0.0, 0.0, 0.0, "None"]

        # Figure of M&M
        sige3 = coeffs[0] + coeffs[1] * E_kev + coeffs[2] * E_kev2
        # cross section per hydrogen atom /cm2
        sig = sige3 / E_kev3  # * 1e-24

        # NHx is given in 1e22 cm-2, and sig should be multiplied by 1e-24
        Tau_gas[i] = sig * NHx * 1e-2

    Trans_gas = np.exp(-Tau_gas)

    Trans_gas[Trans_gas < 1e-5] = 0
    Trans_gas[Trans_gas > 1] = 1

    return Trans_gas
