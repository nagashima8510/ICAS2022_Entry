import numpy as np

#unit conversion factors
#length
m2ft=3.28084 #from meter to feet
ft2m=0.3048 #from feet to meter
ft2km=0.0003048 #from feet to kilo-meter

#speed
kt2mps=0.514444 #from knot to meter-per-sec
mps2kt=1.94384 #from meter-per-sec to knot

#angle
rad2deg=180/np.pi #from radian to degree
deg2rad=np.pi/180 #from degree to radian

#weight
lb2kg=0.453592 #from pond to kilo-gram
kg2lb=2.20462 #from kilo-gram to pond

#time
sec2min=1/60 #from sec to min
min2sec=60 #from min to sec

#atmospheric Constants
gamma=1.4 #specific heat ratio
mu=(gamma-1)/gamma
M=28.96442 #[kilo-gram/kilo-mole] molecular weight
R0=8.31446261815324 #[J/(K mole)] gas constant
R=287.05287 #[m2/(K kilo-gram2)] gas constant
g0=9.80665 #[m/s2] gravitational acceleration