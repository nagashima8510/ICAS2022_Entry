#3Dimensional Flight Path Optimization Program.
#風速分布は計算空間で一様
#空港までの距離を非常時の安全性評価指標に考慮

import concurrent.futures
import pandas as pd
import numpy as np
import math
import datetime
import requests

from A320 import *
from Constants import *
from ISA_Functions import *
from Job_input3d_9A import *
from Utility_Functions import *
from Runway import *

#Down Range定義
x=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[km]
y=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[km]

#Altitude定義
z=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[m]

#対気速度
vcas=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[CAS_m/s]

#真対気速度
v=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[TAS_m/s]

#風向風速分布
head_wind=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[m/s]
cross_wind=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[m/s]
wind_speed=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[m/s]

#速度制限がVMOかMMOかの仕分け
M_vcas_trans=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1), dtype=np.bool_) #Mach数評定の場合True

#有効接点識別タグ
Effective_Node=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1), dtype=np.bool_)

#親接点の格納配列にかかる相対位置
Previous_Node_i=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1), dtype=np.int32)
Previous_Node_j=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1), dtype=np.int32)
Previous_Node_k=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1), dtype=np.int32)
Previous_Node_l=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1), dtype=np.int32)

#重量
W=np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1), dtype=np.float64) #[kg]

#コスト配列
#Fuel_Burn = np.full((imax+1,jmax+1,kmax+1,lmax+1),99999.9)
Flight_Time = np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[sec] 飛行時間累積
Fuel_Burn = np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[kg] 燃料消費累積
#Flight_Quality = np.full((imax+1,jmax-jmin+1,kmax+1,lmax+1),999999.0) #[max_m/s] 機体動揺最大値
Flight_Quality = np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[max_m/s] 機体動揺最大値
Flight_Distance = np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[km] 飛行距離
Moving_Distance = np.zeros((imax+1,jmax-jmin+1,kmax+1,lmax+1)) #[km] 移動距離

#Gridデータ
print("Grid")
print("imin=",imin,"xmin=",xmin,"/ imax=",imax,"xmax=",xmax,"km")
print("jmin=",jmin,"ymin=",ymin,"/ jmax=",jmax,"ymax=",ymax,"km")
print("kmin=",kmin,"zmin=",zmin,"/ kmax=",kmax,"zmax=",zmax,"m")
print("lmin=",lmin,"vmin=",cas2tas(vmin,zmin),"/ lmax=",lmax,"vmax=",cas2tas(vmax,zmax),"m/s")

print("Initialize Grind Point Data")
for i in range(imin, imax+1):
    print("i=",i)
    for j in range(jmin, jmax+1):
        for k in range(kmin,kmax+1):
            for l in range(lmin,lmax+1):
                x[i][j][k][l]=xmin+Delta_x*(i-imin) #[km]
                y[i][j][k][l]=ymin+Delta_y*(j-jmin) #[km]
                z[i][j][k][l]=zmin+Delta_z*(k-kmin) #[m]
                vcas[i][j][k][l]=vmin+Delta_v*(l-lmin) #[CAS_m/s]
                v[i][j][k][l]=cas2tas(vcas[i][j][k][l],z[i][j][k][l]) #[TAS_m/s]
                wind_speed[i][j][k][l]=(((45000-z[i][j][k][l]*m2ft)/500*ft2m)/(10+((45000-z[i][j][k][l]*m2ft)/500*ft2m)**1.6)-0.15)/0.05*100/135.2274*Wind_Speed #[m/s]
                head_wind[i][j][k][l]=-wind_speed[i][j][k][l]*math.cos(math.radians(Wind_Direction))
                cross_wind[i][j][k][l]=-wind_speed[i][j][k][l]*math.sin(math.radians(Wind_Direction))
                
print("Initialize Air-Space Condition")
for i in range(imin, imax+1):
    print("i=",i)
    for j in range(jmin, jmax+1):
        for k in range(kmin,kmax+1):
            for l in range(lmin,lmax+1):
                x_cur=x[i][j][k][l] #ダウンレンジ
                z_cur=z[i][j][k][l] #圧力高度
                vcas_cur=vcas[i][j][k][l] #補正対気速度
                T_cur=Alt2Temp(z_cur) #[K]温度
                p_cur=Alt2Press(z_cur) #[Pa]気圧
                rho_cur=Alt2Dens(z_cur) #[kg/m3]密度
                a_cur=np.sqrt(gamma*R*T_cur) #[m/s]音速
                v_cur=cas2tas(vcas_cur,z_cur) #[m/s]真対気速度
                Mach_cur=v_cur/a_cur #Mach数
                
                p_trans=p0*((1+(gamma-1)/2*(vcas_cur/a0)**2)**(1/mu)-1)/((1+(gamma-1)/2*Mach_cur**2)**(1/mu)-1) #[Pa]
                
                if p_trans>=p_trop:
                    z_trans=T0/betaT*((p_trans/p0)**(-betaT*R/g0)-1) #[m]
                else:
                    z_trans=z_trop-R*T_trop/g0*math.log(p_trans/p_trop) #[m]
                
                if z_trans>=z_cur:
                    M_vcas_trans[i][j][k][l]=False #CAS評定
                else:
                    M_vcas_trans[i][j][k][l]=True #Mach数評定
                        
    #iminにおいて始点ノードを有効化
    m_start=mref #重量は初期値を代入
    Effective_Node[i_start][j_start][k_start][l_start]=True
    W[i_start][j_start][k_start][l_start]=m_start*1000 #[kg] 始点ノードにおける重量を指定
    
    #imaxにおいて終点ノードのみ有効化
    #Effective_Node[i_end][j_end][k_end][l_end]=True

#始点
print("Start Poit")
print("i_start=",i_start,"x=",x[i_start][j_start][k_start][l_start],"km")
print("j_start=",j_start,"y=",y[i_start][j_start][k_start][l_start],"km")
print("k_start=",k_start,"z=",z[i_start][j_start][k_start][l_start],"m")
print("l_start=",l_start,"vcas=",vcas[i_start][j_start][k_start][l_start],"CAS_m/s")
print("l_start=",l_start,"vtas=",v[i_start][j_start][k_start][l_start],"TCAS_m/s")
#終点
print("End Point")
print("i_end=",i_end,"x=",x[i_end][j_end][k_end][l_end],"km")
print("j_end=",j_end,"x=",y[i_end][j_end][k_end][l_end],"km")
print("k_end=",k_end,"z=",z[i_end][j_end][k_end][l_end],"m")
print("l_end=",l_end,"vcas=",vcas[i_end][j_end][k_end][l_end],"CAS_m/s")
print("l_end=",l_end,"vtas=",v[i_end][j_end][k_end][l_end],"TCAS_m/s")

#滑走路情報の読み込み
rw_table=read_runway()
print(rw_table)

start = time.time()

#最適経路探索ルーチン
#X方向にiを順にインクリメントしていく
#X方向に、ID=iのInitial接点とi+1のFinal接点の組み合わせに対し…
#接点間の高度/速度からコストを算出、コストミニマムとする経路を選定
for i in range(imin, imax):
    print("i=",i)
    for j in range(jmin, jmax+1):
        for k in range(kmin,kmax+1):
            for l in range(lmin,lmax+1):
                if Effective_Node[i][j][k][l]==True: #有効なInitial接点の場合、コストを算出する
                    x_init=x[i][j][k][l] #[km]
                    y_init=y[i][j][k][l] #[km]
                    z_init=z[i][j][k][l] #[m]
                    W_init=W[i][j][k][l] #[kg]
                    v_init=v[i][j][k][l] #[TAS_m/s]
                    T_init=Alt2Temp(z_init) #[K]
                    p_init=Alt2Press(z_init) #[Pa]
                    rho_init=Alt2Dens(z_init) #[kg/m3]
                    vcas_init=vcas[i][j][k][l] #[CAS_m/s]
                    a_init=np.sqrt(gamma*R*T_init) #[m/s]
                    Mach_init=v_init/a_init
                    q_init=0.5*rho_init*v_init**2 #[Pa] 動圧@出発点
                    CL_init=W_init*g0/q_init/Sref #揚力係数@出発点
                    CD_init=CD0_CR+CD2_CR*CL_init**2 #抵抗係数@出発点
                    Drag_init=CD_init*q_init*Sref
                    head_wind_init=head_wind[i][j][k][l] #[m/s]
                    cross_wind_init=cross_wind[i][j][k][l] #[m/s]
                    wind_speed_init=wind_speed[i][j][k][l] #[m/s]
                    #print("♪♪♪♪♪♪♪♪♪♪","j_init=",j,"k_init=",k,"l_init=",l,"W=",W_init,"TAS=",v_init,"CAS=",vcas_init,"Mach=",Mach_init,"CL1=",CL_init,"CD1=",CD_init,"Drag1=",Drag_init)
                    for jj in range(jmin, jmax+1):
                        for kk in range(kmin,kmax+1):
                            for ll in range(lmin,lmax+1):
                                
                                #if (i+1==i_end)and((jj!=j_end)or(kk!=k_end)or(ll!=l_end)):
                                #    break
                                
                                #elif (i+1==i_end)and((jj==j_end)and(kk==k_end)and(ll==l_end)):
                                    #print("Last Node!",i+1,jj,kk,ll)
                                        
                                #次のNodeとして現実的でない場合は、無駄に評価ないようBreakする。
                                #print("Final Node",i+1,jj,kk,ll)
                                x_fin=x[i+1][jj][kk][ll] #[km]
                                y_fin=y[i+1][jj][kk][ll] #[km]
                                
                                #横の変化は、最大Thetaを超えない
                                if (y_init+Delta_x*np.tan(math.radians(Theta_Max*1.05))<y_fin)or(y_init-Delta_x*np.tan(math.radians(Theta_Max*1.05))>y_fin):
                                    #print("Y:",y_init,y_fin)
                                    break
                                
                                z_fin=z[i+1][jj][kk][ll] #[m]
                                
                                #高度の変化は、（高度変化）/（水平距離）=1500(m)/20(km)を超えない
                                if (z_init+Delta_x*(1500/20)<z_fin)or(z_init-Delta_x*(1500/(20))>z_fin):
                                    #print("Z:",z_init,z_fin)
                                    break
                                
                                v_fin=v[i+1][jj][kk][ll] #[TAS_m/s]
                                
                                #TAS速度の変化は、（速度変化）/（水平距離）が、+20(m/s)/20(km)から-30(m/s)/20(km)の範囲
                                if (v_init+Delta_x*(20/20)<v_fin)or(v_init-Delta_x*(30/20)>v_fin):
                                    #print("V:",v_init,v_fin)
                                    break
                                #print("Init:",j,k,l,"Fin:",jj,kk,ll)
                                
                                #次のNodeの諸量を設定
                                T_fin=Alt2Temp(z_fin) #[K]
                                p_fin=Alt2Press(z_fin) #[Pa]
                                rho_fin=Alt2Dens(z_fin) #[kg/m3]
                                vcas_fin=vcas[i+1][jj][kk][ll] #[CAS_m/s]
                                a_fin=np.sqrt(gamma*R*T_fin) #[m/s]
                                Mach_fin=v_fin/a_fin
                                q_fin=0.5*rho_fin*v_fin**2 #[Pa] 動圧@到着点
                                M_vcas_fin=M_vcas_trans[i+1][jj][kk][ll] #最大速度制限の評定(Trueのときマッハ数)
                                head_wind_fin=head_wind[i+1][jj][kk][ll] #[m/s]
                                cross_wind_fin=cross_wind[i+1][jj][kk][ll] #[m/s]
                                wind_speed_fin=wind_speed[i+1][jj][kk][ll] #[m/s]
                                
                                #消費燃料を推算
                                
                                #接点間平均風速成分
                                cross_wind_temp = 0.5*(cross_wind_init + cross_wind_fin)
                                head_wind_temp = 0.5*(head_wind_init + head_wind_fin)
                                wind_speed_temp = 0.5*(wind_speed_init + wind_speed_fin)
                                
                                #座標変換
                                phi = CoordTrans1_angle(x_init,y_init,x_fin,y_fin,0.0,0.0,1.0,0.0)
                                head_wind_temp2 = -head_wind_temp*np.cos(phi)+cross_wind_temp*np.sin(phi)
                                cross_wind_temp2 = head_wind_temp*np.sin(phi)+cross_wind_temp*np.cos(phi)
                                #始点、終点の接点における諸量の算出完了
                                
                                #接点間距離
                                Ground_Distance = np.sqrt(\
                                                   (x_fin-x_init)**2+(y_fin-y_init)**2+((z_fin-z_init)/1000)**2\
                                                  )*1000 #[m] 距離
                                #print(Distance,"m")
                                
                                #接点間平均飛行速度
                                v_airspeed = 0.5*(v_init + v_fin) #[m/s]
                                #print("Air Speed",v_airspeed)
                                
                                #接点間対地速度
                                v_groundspeed = 0.5*(\
                                                    -2.0*head_wind_temp2+\
                                                    np.sqrt(4.0*head_wind_temp2**2-4.0*(head_wind_temp2**2+cross_wind_temp2**2-v_airspeed**2))\
                                                   ) #[m/s]
                                #print("Ground Speed",v_groundspeed)
                                
                                #接点間飛行距離（等加速度と仮定）
                                Distance = v_airspeed * Ground_Distance / v_groundspeed
                                #所要時間
                                Time_sec=Ground_Distance/v_groundspeed #[sec]
                                #Time_sec=2*Distance/(v_init-wind_speed_init+v_fin-wind_speed_fin) #[sec]
                                Time_min=Time_sec*sec2min #[min]
                                #print(Time_min,"min")
                                
                                #加速度
                                Accel=(v_fin**2-v_init**2)/(2*Distance) #[m/s2]
                                
                                #必要推力
                                ccc=Cf1/1000*(1+((v_init+v_fin)*0.5*mps2kt)/Cf2)*Time_min #[kg/N]
                                aa=CD2_CR*g0**2*ccc**2/q_fin/Sref
                                bb=-(2+2*CD2_CR*g0**2*W_init*ccc/q_fin/Sref+ccc*v_fin**2/Distance+2*ccc*g0*z_fin/Distance)
                                cc=CD0_CR*q_init*Sref+CD2_CR*W_init**2*g0**2/q_init/Sref+CD0_CR*q_fin*Sref+CD2_CR*g0**2*W_init**2/q_fin/Sref+W_init*(v_fin**2-v_init**2+2*g0*(z_fin-z_init))/Distance

                                Thr1=(-bb+np.sqrt(bb**2-4*aa*cc))/(2*aa) #[kN]
                                Thr2=(-bb-np.sqrt(bb**2-4*aa*cc))/(2*aa) #[kN]
                                
                                #最大推力（常時）
                                Thr_max_climb=(CTC1*(1-(z_fin*m2ft)/CTC2+CTC3*(z_fin*m2ft)**2))*CTcr
                                #print("Max Thrust=",Thr_max_climb)
                                
                                #降下時推力
                                if z_init<Hpdes*ft2m:
                                    Thr_des=Thr_max_climb*CTdes_low
                                else:
                                    Thr_des=Thr_max_climb*CTdes_high
                                #print("Thrust_nominal=",Thr2,"N"," Thrust_descent=",Thr_des,"N")
                                Thr=max(Thr2,Thr_des)
                                
                                #燃料消費
                                Fuel_Consumption_nom=Thr2*ccc #[kg] 必要推力を使用した場合の消費燃料
                                Fuel_Consumption_des=Cf3*(1-(z_init+z_fin)*0.5*m2ft/Cf4)*Time_min #[kg] 推力に余力がある場合(アイドル)の消費燃料
                                
                                Fuel_Consumption=max(Fuel_Consumption_nom,Fuel_Consumption_des) #[kg] 燃料消費は、大きい方の値
                                #print("Fuel Consumption=",Fuel_Consumption_nom,Fuel_Consumption_des)
                                
                                #print("♪♪♪♪♪♪♪♪♪♪　k_fin=",kk,"l_fin=",ll,"W_fin=",W[i+1][jj][kk][ll],"TAS_fin=",v_fin,"CAS_fin=",vcas_fin)
                                
                                W_fin=W_init-Fuel_Consumption
                                CL_fin=W_fin*g0/q_fin/Sref #揚力係数@到着点
                                CD_fin=CD0_CR+CD2_CR*CL_fin**2 #抵抗係数@到着点
                                Drag_fin=CD_fin*q_fin*Sref
                                
                                #print("Ed1=",0.5*W_init*g0*v_init**2,"Ed2=",0.5*W_fin*g0*v_fin**2,"Delta_Ed=",0.5*W_fin*g0*v_fin**2-0.5*W_init*g0*v_init**2)
                                #print("Ep1=",W_init*g0*z_init,"Ep2=",W_fin*g0*z_fin,"Delta_Ep=",W_fin*g0*z_fin-W_init*g0*z_init)
                                
                                #print("Drag=",0.5*(Drag_init+Drag_fin),"Drag1=",Drag_init,"Drag2=",Drag_fin,"N")
                                #print("Thrust=",Thr,"N")
                                #print((Thr-0.5*(Drag_init+Drag_fin))*Distance,"Nm")
                                
                                #次のNodeでの危険性を推定
                                emergency_=10000.0 #距離が閾値ぴったりの場合、1.0になる危険度指数
                                
                                #print(x_fin,y_fin,z_fin)
                                
                                for n in range(0,len(rw_table)):
                                    Name_rw=rw_table.iat[n,2]
                                    L_rw=rw_table.iat[n,5]
                                    z_rw=rw_table.iat[n,9]
                                    x_rw=rw_table.iat[n,10]
                                    y_rw=rw_table.iat[n,11]
                                    #print(n,Name_rw,x_rw,y_rw,z_rw,L_rw)
                                    threshold=(z_fin-z_rw)/1000*CL_fin/CD_fin
                                    distance_to_rw=np.sqrt((x_fin-x_rw)**2+(y_fin-y_rw)**2)
                                    emergency_temp=(np.exp(distance_to_rw/threshold)-1.0)/(np.exp(1.0)-1.0)
                                    if emergency_temp<=emergency_:
                                        emergency_=emergency_temp
                                
                                #print("Emergency=",emergency_)
                                
                                #最適経路を判定
                                if Thr_max_climb>=Thr: #搭載エンジンによりこの高度まで上昇可能
                                    #print("OK: Less than Max Thrust"," CAS=",vcas_fin,"V_stall=",Vstall_CR)
                                    if vcas_fin>=Vstall_CR: #失速速度以上
                                        #print("OK: more than Stall speed.")
                                        if ((M_vcas_fin==True)and(Mach_fin<MMO))or((M_vcas_fin==False)and(vcas_fin<VMO)): #速度制限以下
                                            #print("OK: Less than maximum speed constraint") #この場合、次の接点まで飛行化可能
                                            #print("Weight=",W[i+1][jj][kk][ll],W_init-Fuel_Consumption)
                                            #if W[i+1][jj][kk][ll]<W_init-Fuel_Consumption: #この場合が、最小燃料消費（機体重量最大）の可能性あり
                                            
                                            if Effective_Node[i+1][jj][kk][ll] == False:
                                                Effective_Node[i+1][jj][kk][ll] = True
                                                #print("ok")
                                                W[i+1][jj][kk][ll]=W_init-Fuel_Consumption #[kg]
                                                Flight_Time[i+1][jj][kk][ll]=Flight_Time[i][j][k][l]+Time_sec #[sec]
                                                Fuel_Burn[i+1][jj][kk][ll]=Fuel_Burn[i][j][k][l]+Fuel_Consumption #[kg]
                                                Flight_Quality[i+1][jj][kk][ll]=Flight_Quality[i][j][k][l]+emergency_
                                                Flight_Distance[i+1][jj][kk][ll] = Flight_Distance[i][j][k][l]+Distance #[km] 飛行距離
                                                Moving_Distance[i+1][jj][kk][ll] = Moving_Distance[i][j][k][l]+Ground_Distance #[km] 移動距離
                                                Previous_Node_i[i+1][jj][kk][ll]=i-(i+1)
                                                Previous_Node_j[i+1][jj][kk][ll]=j-jj
                                                Previous_Node_k[i+1][jj][kk][ll]=k-kk
                                                Previous_Node_l[i+1][jj][kk][ll]=l-ll
                                                #print(jj,kk,ll,"W",W[i+1][jj][kk][ll],"FT",Flight_Time[i+1][jj][kk][ll],"FB",Fuel_Burn[i+1][jj][kk][ll],"FQ",Flight_Quality[i+1][jj][kk][ll])
                                                #Object_Parameter_temp=w1*(Fuel_Burn[i+1][jj][kk][ll])+w2*(Flight_Quality[i+1][jj][kk][ll])
                                                #print("OK_0",jj,kk,ll,emergency_)
                                            else:
                                                #print("OKOK")
                                                if w1*(Fuel_Burn[i+1][jj][kk][ll])+w2*(Flight_Quality[i+1][jj][kk][ll])>w1*(Fuel_Burn[i][j][k][l]+Fuel_Consumption)+w2*(Flight_Quality[i][j][k][l]+emergency_):
                                                    #print("OKOK")
                                                    W[i+1][jj][kk][ll]=W_init-Fuel_Consumption #[kg]
                                                    Flight_Time[i+1][jj][kk][ll]=Flight_Time[i][j][k][l]+Time_sec #[sec]
                                                    Fuel_Burn[i+1][jj][kk][ll]=Fuel_Burn[i][j][k][l]+Fuel_Consumption #[kg]
                                                    Flight_Quality[i+1][jj][kk][ll]=Flight_Quality[i][j][k][l]+emergency_
                                                    Flight_Distance[i+1][jj][kk][ll] = Flight_Distance[i][j][k][l]+Distance #[km] 飛行距離
                                                    Moving_Distance[i+1][jj][kk][ll] = Moving_Distance[i][j][k][l]+Ground_Distance #[km] 移動距離
                                                    Previous_Node_i[i+1][jj][kk][ll]=i-(i+1)
                                                    Previous_Node_j[i+1][jj][kk][ll]=j-jj
                                                    Previous_Node_k[i+1][jj][kk][ll]=k-kk
                                                    Previous_Node_l[i+1][jj][kk][ll]=l-ll
                                                    #print(jj,kk,ll,"W",W[i+1][jj][kk][ll],"FT",Flight_Time[i+1][jj][kk][ll],"FB",Fuel_Burn[i+1][jj][kk][ll],"FQ",Flight_Quality[i+1][jj][kk][ll])
                                                    #Object_Parameter_temp=w1*(Fuel_Burn[i+1][jj][kk][ll])+w2*(Flight_Quality[i+1][jj][kk][ll])
                                                    #print("OK",jj,kk,ll,emergency_)
                                                    
print("Completed")
print("Final Weight=",W[i_end][j_end][k_end][l_end])
print("Initial Weight=",W[i_start][j_start][k_start][l_start])

end = time.time()
delta = end - start
print('処理時間:{}s'.format(round(delta,3)))

#最適経路を記録
import datetime
print("Record Solution")

path_x=np.zeros(imax+1, dtype=np.float64)
path_y=np.zeros(imax+1, dtype=np.float64)
path_z=np.zeros(imax+1, dtype=np.float64)
path_v=np.zeros(imax+1, dtype=np.float64)
path_vcas=np.zeros(imax+1, dtype=np.float64)
path_W=np.zeros(imax+1, dtype=np.float64)
path_FT=np.zeros(imax+1, dtype=np.float64)
path_FB=np.zeros(imax+1, dtype=np.float64)
path_FD=np.zeros(imax+1, dtype=np.float64)
path_MD=np.zeros(imax+1, dtype=np.float64)
path_FQ=np.zeros(imax+1, dtype=np.float64)

ipn=i_end
jpn=j_end
kpn=k_end
lpn=l_end

for i in range(imin,imax+1):
    path_x[i]=x[ipn][jpn][kpn][lpn]
    path_y[i]=y[ipn][jpn][kpn][lpn]
    path_z[i]=z[ipn][jpn][kpn][lpn]
    path_v[i]=v[ipn][jpn][kpn][lpn]
    path_vcas[i]=vcas[ipn][jpn][kpn][lpn]
    path_W[i]=W[ipn][jpn][kpn][lpn]
    path_FT[i]=Flight_Time[ipn][jpn][kpn][lpn]
    path_FB[i]=Fuel_Burn[ipn][jpn][kpn][lpn]
    path_FD[i]=Flight_Distance[ipn][jpn][kpn][lpn]
    path_MD[i]=Moving_Distance[ipn][jpn][kpn][lpn]
    path_FQ[i]=Flight_Quality[ipn][jpn][kpn][lpn]
    
    ipn_tmp=ipn+Previous_Node_i[ipn][jpn][kpn][lpn]
    jpn_tmp=jpn+Previous_Node_j[ipn][jpn][kpn][lpn]
    kpn_tmp=kpn+Previous_Node_k[ipn][jpn][kpn][lpn]
    lpn_tmp=lpn+Previous_Node_l[ipn][jpn][kpn][lpn]
    
    ipn=ipn_tmp
    jpn=jpn_tmp
    kpn=kpn_tmp
    lpn=lpn_tmp

Solution = pd.DataFrame(
    data={'X(km)':path_x, 
          'Y(km)':path_y,
          'Z(m)':path_z,
          'v(m/s_TAS)':path_v,
          'vcas(m/s_CAS)':path_vcas,
          'W(kg)':path_W,
          'Flight Time(sec)':path_FT,
          'Fuel Burn(kg)':path_FB,
          'Flight Distance(km)':path_FD,
          'Moving Distance(km)':path_MD,
          'Flight Quality(ND)':path_FQ})


dt_now = datetime.datetime.now()
dt_now_txt = dt_now.strftime('%Y%m%d_%H%M%S')
print(dt_now_txt)

file_name2="Result_"+dt_now_txt+".csv"
Solution.to_csv(file_name2)



#プロット作成
print("Output Graph")
import matplotlib.pyplot as plt

Solution = pd.read_csv(file_name2, dtype=str)

Solution['X(km)'] = Solution['X(km)'].astype(float)
Solution['Y(km)'] = Solution['Y(km)'].astype(float)
Solution['Z(m)'] = Solution['Z(m)'].astype(float)
Solution['v(m/s_TAS)'] = Solution['v(m/s_TAS)'].astype(float)
Solution['vcas(m/s_CAS)'] = Solution['vcas(m/s_CAS)'].astype(float)
Solution['W(kg)'] = Solution['W(kg)'].astype(float)
Solution['Flight Time(sec)'] = Solution['Flight Time(sec)'].astype(float)
Solution['Fuel Burn(kg)'] = Solution['Fuel Burn(kg)'].astype(float)
Solution['Flight Distance(km)'] = Solution['Flight Distance(km)'].astype(float)
Solution['Moving Distance(km)'] = Solution['Moving Distance(km)'].astype(float)
Solution['Flight Quality(ND)'] = Solution['Flight Quality(ND)'].astype(float)

path_Delta_W=np.zeros(imax, dtype=np.float64)
path_x2=np.zeros(imax, dtype=np.float64)

for i in range(imin,imax+1):
    path_x[i]=Solution.at[i,'X(km)']
    path_y[i]=Solution.at[i,'Y(km)']
    path_z[i]=Solution.at[i,'Z(m)']
    path_v[i]=Solution.at[i,'v(m/s_TAS)']
    path_vcas[i]=Solution.at[i,'vcas(m/s_CAS)']
    path_W[i]=Solution.at[i,'W(kg)']
    path_FT[i]=Solution.at[i,'Flight Time(sec)']
    path_FB[i]=Solution.at[i,'Fuel Burn(kg)']
    path_MD[i]=Solution.at[i,'Moving Distance(km)']
    path_FD[i]=Solution.at[i,'Flight Distance(km)']
    path_FQ[i]=Solution.at[i,'Flight Quality(ND)']

    if i>imin:
        path_Delta_W[i-1]=path_W[i]-path_W[i-1]
        path_x2[i-1]=(path_x[i]+path_x[i-1])*0.5

#平面図
fig2 = plt.figure(figsize=(24, 12))
plt.title("Air Space")
ax5= fig2.add_subplot(111)

ax5.grid(True)
ax5.set_xlabel("Down Range(km)")
ax5.set_ylabel("Lateral(km)")
ax5.set_ylim(-200,200)
ax5.plot(path_x,path_y,color='r',marker="x")

#側面図
fig1 = plt.figure(figsize=(24, 12))
plt.title("Air Space")

ax1 = fig1.add_subplot(211)
ax1.grid(True)
ax1.set_xlabel("Down Range(km)")
ax1.set_ylabel("Altitude(m)")
ax1.set_ylim(3000,12000)
ax2 = ax1.twinx()
ax2.grid(False)
ax2.set_ylabel("Velocity(m/s)")
ax2.set_ylim(180,300)

ax3 = fig1.add_subplot(212)
ax3.grid(True)
ax3.set_xlabel("Down Range(km)")
ax3.set_ylabel("Weight(kg)")
ax3.set_ylim(0,300000)
ax4 = ax3.twinx()
ax4.grid(False)
ax4.set_ylabel("Fuel Burn(kg)")
ax4.set_ylim(0,20000)

ax1.plot(path_x,path_z,color='r',marker="x",label="Altitude(m)")
ax1.tick_params(axis='y', colors="r")
for i in range(imin,imax+1):
    ax1.annotate(int(path_W[i]), (path_x[i],path_z[i]))

ax2.plot(path_x,path_v,color='b',marker="x",label="Air Speed(m/s)")
ax2.tick_params(axis='y', colors="b")
    
ax3.plot(path_x,path_W,color='g',marker="x",label="Weight(kg)")
ax3.tick_params(axis='y', colors="g")

ax4.bar(path_x2,path_Delta_W,color='b', alpha=0.4, width=20, label="Fuel Burn(kg)")
ax4.tick_params(axis='y', colors="b")

plt.show()