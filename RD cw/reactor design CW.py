# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 19:41:31 2023

@author: HUAWEI
"""
import math
import matplotlib.pyplot as plt
#x SO2 y N2 z O2
p0=1 #atomsphere pressure
Cp1x=50.87 #KJ/(kmol.k) as below
Cp1y=31.2
Cp1z=32.94
Cp1c=97.4
H=-98900   
F=0.654

 
S=3.8**2*3.14 

void=0.197 
Dp=0.0056 
gc=1


def f(T0,w,x0,y0,z0,c0):
    # 1st conversion
    m=F*x0*64+F*y0*28+F*z0*32
    c=0
    Density=x0*2.63+y0*1.3+z0*1.784 
    mu=12.82*10**(-12)*x0+14.72*10**(-12)*y0+32.66*10**(-6)*z0+c*1.3*10**-3
    row=650
    c1=c0 
    x1=x0 
    y1=y0 
    z1=z0
    X1=0.75
    
    Ft=F
    W=0
    h=10
    P=p0
    x00=x0
    while W <w:
        W=W+h
        # pressure drop
        mu=12.82*10**(-12)*x0+14.72*10**(-12)*y0+32.66*10**(-6)*z0+c1*1.3*10**-3
        T1=T0+(X1*F*x0*(-H))/((Cp1x*x0+Cp1y*y0+Cp1z*z0+Cp1c*c1)*Ft)
        # T1=T0+X1*(-H)/(Cp1x*x1+Cp1y*y1+Cp1z*z1+Cp1c*c1))
        def p(P):
            Beta=-(m/S)/Density/gc/Dp*((1-void)/void**3)*((150*(1-void)*mu/Dp+1.75*(m/S)))
            alpha=-Beta/(S*(1-void)*row*p0)*p0/(P/p0)
            #y=P/p0
            deta= -alpha/2*p0/(P/p0)*T1/T0*Ft/F
            return deta
        P=p0+p(P)*(W/(row*(1-void)*S))
    
    
        #Stoichiometry
        t=x0*(1-X1)+y0+z0+X1*x0/2
        Ft=F*t   
        x1=x0*(1-X1)/t
        y1=y0/t
        z1=(F*z0-F*x0*X1/2)/Ft
        c1=x0*X1/t
                                      
    
        
        
        K1=20.6*math.e**(-69.5*1000/T1/8.314)
        Kp=math.e**(-10.68+11300/T1)
        B=P*c1/(Kp*P*x1*(P*z1)**0.5)
        rso2=K1*P*x1*(1-B)
        
        
        X1=X1+rso2/F/x00*h
        
    return x1,y1,z1,c1,T1,X1,P
f(695,6000,0.022,0.852,0.06,0.065)
# plt.figure()
# x = list(range(1, 10000))
# y = [f(695,i,0.085,0.825,0.09,0)[5] for i in x]
# plt.plot(x,y)
# plt.xlabel('Mass of atalyst(kg)')
# plt.ylabel('Conversion')
# plt.title('Catalyst Profile')
# plt.show()
