import numpy as np
from Constants import *

betaT=-0.0065 #[K/m]
T0=288.15 #[K]
p0=101325 #[Pa]
rho0=1.225 #[kg/m3]
a0=(gamma*R*T0)**0.5 #[m/s]

z_trop=11000.0 #[ft]
T_trop=T0+z_trop*betaT #[K]
p_trop=p0*(T_trop/T0)**(-g0/(betaT*R)) #[Pa]
rho_trop=p_trop/(R*T_trop) #[kg/m3]

def Alt2Temp(Alt): #Alt[m]
    if Alt<z_trop:
        T=T0+Alt*betaT #[K]
    else:
        T=T0+z_trop*betaT #[K]

    return T #[K]

def Alt2Press(Alt): #Alt[m]
    if Alt<z_trop:
        T=T0+Alt*betaT #[K]
        p=p0*(T/T0)**(-g0/betaT/R) #[Pa]
    else:
        T=T0+z_trop*betaT #[K]
        p=p0*(T/T0)**(-g0/betaT/R)*np.exp(-g0/R/T*(Alt-z_trop)) #[Pa]
    
    return p #[Pa]
    
def Alt2Dens(Alt): #Alt[m]
    if Alt<z_trop:
        T=T0+Alt*betaT #[K]
        p=p0*(T/T0)**(-g0/betaT/R) #[Pa]
    else:
        T=T0+z_trop*betaT #[K]
        p=p0*(T/T0)**(-g0/betaT/R)*np.exp(-g0/R/T*(Alt-z_trop)) #[Pa]
    rho=p/R/T #[kg/m3]

    return rho #[kg/m3]

def tas2cas(tas,Alt): #tas[m/s],Alt[m]
    if Alt<z_trop:
        T=T0+Alt*betaT #[K]
        p=p0*(T/T0)**(-g0/betaT/R) #[Pa]
    else:
        T=T_trop #[K]
        p=p_trop*np.exp(-g0/R/T*(Alt-z_trop)) #[Pa]
    rho=p/R/T #[kg/m3]
    cas=(2*p0/mu/rho0*((1+p/p0*((1+mu/2*rho/p*tas**2)**(1/mu)-1))**mu-1))**(1/2) #[CAS_m/s]
    return cas #[m/s]

def cas2tas(cas,Alt):
    if Alt<z_trop:
        T=T0+Alt*betaT #[K]
        p=p0*(T/T0)**(-g0/betaT/R) #[Pa]
    else:
        T=T_trop #[K]
        p=p0*(T/T0)**(-g0/betaT/R)*np.exp(-g0/R/T*(Alt-z_trop)) #[Pa]
    rho=p/R/T #[kg/m3]
    tas=(2*p/mu/rho*((1+p0/p*((1+mu/2*rho0/p0*cas**2)**(1/mu)-1))**mu-1))**(1/2) #[TAS_m/s]
    return tas #[m/s]