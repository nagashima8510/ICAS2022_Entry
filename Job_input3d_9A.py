import numpy as np
import math
from Constants import *
from ISA_Functions import *
from Utility_Functions import *

#Weight parameter
w1=0.0005
w2=0.9995

###風向風速設定
Wind_Speed=100.0 #[m/s] 風速
Wind_Direction=270.0 #[deg] 始点

###格子設定

#HND-TAURA
Dpt_Lat=35.31280556 #[deg]
Dpt_Long=139.7464722 #[deg]
Dpt_Alt=9000 #[ft]

#NHC-FOCUS
Arv_Lat=26.67555556 #[deg]
Arv_Long=127.7120556 #[deg]
Arv_Alt=4000 #[ft]

ja=Dpt_Long
da=Dpt_Lat
jb=Arv_Long
db=Arv_Lat

Distance=gcr_distance(ja,da,jb,db)

#X方向
imax=4
#Delta_x=200.0 #[km] 接点間隔
xmin=0.0 #[km] 始点
xmax=xmin+Distance #[km]　終点
imin=0
Delta_x=Distance/(imax-imin)
#imax=int(imin+(xmax-xmin)/Delta_x)

#Y方向
Delta_Theta=10 #[deg] y方向の刻み幅を設定する方位角
Theta_Max=30.0 #[deg] y方向の解析空間サイズを設定するための最大方位角
y_max=round((xmax-xmin)*0.5*np.tan(math.radians(Theta_Max)),1) #[km]
Delta_y=round(Delta_x*np.tan(math.radians(Delta_Theta)),1) #[km]
#Delta_y=0.1 #[km]
jmin=-int(y_max/Delta_y)
jmax=-jmin
ymin=jmin*Delta_y #[km]
ymax=jmax*Delta_y #[km]


#Z方向
Delta_z_=2500.0 #[ft]
zmin_=4000.0 #[ft]
zmax_=36500.0 #[ft]
kmin=0
kmax=int(kmin+(zmax_-zmin_)/Delta_z_)
Delta_z=Delta_z_*ft2m #[m]
zmin=zmin_*ft2m #[m]
zmax=zmax_*ft2m #[m]

#V方向
Delta_v=10.0 #[m/s_CAS]
vmin=150.0 #[m/s_CAS]
vmax=180.0 #[m/s_CAS]
lmin=0
lmax=int(lmin+(vmax-vmin)/Delta_v)

###始点終点設定

#始点
i_start=imin
j_start=0
k_start=2
l_start=1

#終点
i_end=imax
j_end=0
k_end=0
l_end=1