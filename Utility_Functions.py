import numpy as np
import math

def CoordTrans1_angle(x1,y1,x2,y2,x3,y3,x4,y4): #Coordinate Transformation
    # phi is from 0-180deg
    phi=np.arccos(((x4-x3)*(x2-x1)+(y4-y3)*(y2-y1))/np.sqrt((x2-x1)**2+(y2-y1)**2)/np.sqrt((x4-x3)**2+(y4-y3)**2))
    
    # Classify by cross product
    if (x2-x1)*(y4-y3)-(y2-y1)*(x4-x3) >=0: # 0<=phi<=180deg
        return phi #[rad]
    else: # 180<phi<360deg
        return math.pi*2.0-phi #[rad]
    
#Point A: Latitude,Longitude= ja,da [deg]
#Point B: Latitude,Longitude=jb,db [deg]
#Target Point C: Latitude,Longitude=jc,dc [deg]
def gcr_lon2lat(ja,da,jb,db,dc):
    
    rja = np.radians(ja)
    rda = np.radians(da)
    rjb = np.radians(jb)
    rdb = np.radians(db)

    rdc = np.radians(dc)
    
    fjc = np.sin(rda)*np.sin(rjb-rdc)/np.cos(rda)/np.sin(rjb-rja)+np.sin(rdb)*np.sin(rja-rdc)/np.cos(rdb)/np.sin(rja-rjb)
    
    jc = np.degrees(np.arctan(fjc))
    
    return jc #[deg]

def gcr_distance(ja,da,jb,db):
    
    rja = np.radians(ja)
    rda = np.radians(da)
    rjb = np.radians(jb)
    rdb = np.radians(db)
    
    R = 6378.1 #[km] the earth radius
    d = R*np.arccos(np.sin(rja)*np.sin(rjb)+np.cos(rja)*np.cos(rjb)*np.cos(np.abs(rdb-rda)))
    
    return d #[km]