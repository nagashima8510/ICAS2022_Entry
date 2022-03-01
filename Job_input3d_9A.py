import numpy as np
import math
from Constants import *
from ISA_Functions import *
from Utility_Functions import *

#Weight parameter
w1=0.0005
w2=0.9995

###���������ݒ�
Wind_Speed=100.0 #[m/s] ����
Wind_Direction=270.0 #[deg] �n�_

###�i�q�ݒ�

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

#X����
imax=4
#Delta_x=200.0 #[km] �ړ_�Ԋu
xmin=0.0 #[km] �n�_
xmax=xmin+Distance #[km]�@�I�_
imin=0
Delta_x=Distance/(imax-imin)
#imax=int(imin+(xmax-xmin)/Delta_x)

#Y����
Delta_Theta=10 #[deg] y�����̍��ݕ���ݒ肷����ʊp
Theta_Max=30.0 #[deg] y�����̉�͋�ԃT�C�Y��ݒ肷�邽�߂̍ő���ʊp
y_max=round((xmax-xmin)*0.5*np.tan(math.radians(Theta_Max)),1) #[km]
Delta_y=round(Delta_x*np.tan(math.radians(Delta_Theta)),1) #[km]
#Delta_y=0.1 #[km]
jmin=-int(y_max/Delta_y)
jmax=-jmin
ymin=jmin*Delta_y #[km]
ymax=jmax*Delta_y #[km]


#Z����
Delta_z_=2500.0 #[ft]
zmin_=4000.0 #[ft]
zmax_=36500.0 #[ft]
kmin=0
kmax=int(kmin+(zmax_-zmin_)/Delta_z_)
Delta_z=Delta_z_*ft2m #[m]
zmin=zmin_*ft2m #[m]
zmax=zmax_*ft2m #[m]

#V����
Delta_v=10.0 #[m/s_CAS]
vmin=150.0 #[m/s_CAS]
vmax=180.0 #[m/s_CAS]
lmin=0
lmax=int(lmin+(vmax-vmin)/Delta_v)

###�n�_�I�_�ݒ�

#�n�_
i_start=imin
j_start=0
k_start=2
l_start=1

#�I�_
i_end=imax
j_end=0
k_end=0
l_end=1