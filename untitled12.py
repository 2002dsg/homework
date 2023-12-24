# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 15:11:55 2023

@author: HUAWEI
"""

import math
import numpy as np
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


def f(T0,w,x0,y0,z0,c0,X1):
    # 1st conversion
    m=F*x0*64+F*y0*28+F*z0*32
    Density=x0*2.63+y0*1.3+z0*1.784 
    mu=12.82*10**(-12)*x0+14.72*10**(-12)*y0+32.66*10**(-6)*z0+c0*1.3*10**-3
    row=650
    c1=c0 
    x1=x0 
    y1=y0 
    z1=z0
    X11=0
    #v0=m/Density
    Ft=F
    # X1=(x0*Cp1x*(T1-T0)+y0*Cp1y*(T1-T0)+z0*Cp1z*(T1-T0))/(-H)
    
    
    
    W=0
    h=10
    P=p0
    x00=x0
    while W <w:
        W=W+h
        # pressure drop
        mu=12.82*10**(-12)*x0+14.72*10**(-12)*y0+32.66*10**(-6)*z0+c1*1.3*10**-3
        T1=T0+(X11*F*x0*(-H))/((Cp1x*x0+Cp1y*y0+Cp1z*z0+Cp1c*c1)*Ft)
        # T1=T0+X1*(-H)/(Cp1x*x1+Cp1y*y1+Cp1z*z1+Cp1c*c1))
        # def p(P):
        #     Beta=-(m/S)/Density/gc/Dp*((1-void)/void**3)*((150*(1-void)*mu/Dp+1.75*(m/S)))
        #     alpha=-Beta/(S*(1-void)*row*p0)*p0/(P/p0)
        #     deta= -alpha/2*p0/(P/p0)*T1/T0*Ft/F
        #     return deta
        # P=p0+p(P)*(W/(row*(1-void)*S))
    
    
        #Stoichiometry
        t=x0*(1-X11)+y0+z0+X11*x0/2+c0
        Ft=F*t   
        
        x1=x0*(1-X11)/t
        y1=y0/t
        z1=(F*z0-F*x0*X11/2)/Ft
        c1=(c0+x0*X11)/t
                                      
    
        
        
        K1=20.6*math.e**(-69.5*1000/T1/8.314)
        Kp=math.e**(-10.68+11300/T1)
        B=P*c1/(Kp*P*x1*(P*z1)**0.5)
        rso2=K1*P*x1*(1-B)
        
        
        X11=X11+rso2/F/x00*h
        # x0=x1
        # y0=y1
        # z0=z1
        # c=c1
        # T0=T1
    X1=X11*x0
    return x1,y1,z1,c1,T1,X1,t
def ff(T0,x0,y0,z0,c0):
    x000=x0
    
    x1=f(T0,2500,x0,y0,z0,c0,0)[0]
    y1=f(T0,2500,x0,y0,z0,c0,0)[1]
    z1=f(T0,2500,x0,y0,z0,c0,0)[2]
    c1=f(T0,2500,x0,y0,z0,c0,0)[3]
    X1=f(T0,2500,x0,y0,z0,c0,0)[-1]
    # T1 = f(695,2500,x0,y0,z0,c0,0)[4]
    
    x2=f(695,6000,x1,y1,z1,c1,X1)[0]
    y2=f(695,6000,x1,y1,z1,c1,X1)[1]
    z2=f(695,6000,x1,y1,z1,c1,X1)[2]
    c2=f(695,6000,x1,y1,z1,c1,X1)[3]
    X2=f(695,6000,x1,y1,z1,c1,X1)[-1]
    # T1 = f(695,6000,x1,y1,z1,c1,X1)[4]
    
    x3=f(695,10000,x2,y2,z2,c2,X2)[0]
    y3=f(695,10000,x2,y2,z2,c2,X2)[1]
    z3=f(695,10000,x2,y2,z2,c2,X2)[2]
    c3=f(695,10000,x2,y2,z2,c2,X2)[3]
    X3=f(695,10000,x2,y2,z2,c2,X2)[-1]
    
    x4=f(695,15000,x3,y3,z3,c3,X3)[0]
    y4=f(695,15000,x3,y3,z3,c3,X3)[1]
    z4=f(695,15000,x3,y3,z3,c3,X3)[2]
    c4=f(695,15000,x3,y3,z3,c3,X3)[3]
    T1 = f(695,15000,x3,y3,z3,c3,X3)[4]
    X4=1-x4/x000

    return T1,x4,y4,z4,c4,X4
    
ff(695,0.085,0.825,0.09,0) 
    
# f(700,20000,0.085,0.825,0.09,0,0)
plt.figure()
x = [i for i in range(650, 750)]

y = [ff(i,0.085,0.825,0.09,0)[1] for i in x]
plt.plot(x,y)
plt.xlabel('Feed temperature(K)')
plt.ylabel('Output sulfur dioxide composition (100%)')
plt.title('Effect of Temperature change')
plt.show()
# plt.plot(x,y)
# plt.xlabel('Feed temperature(K)')
# plt.ylabel('Output temperature(K)')
# plt.title('Effect of Temperature change')
# plt.show()
# f(695,2500,0.085,0.825,0.09,0)
# f(695,6000,0.022,0.852,0.06,0.065,0)
# f(695,10000,0.0085,0.858,0.0535,0.08)
# f(695,15000,0.0031,0.86,0.051,0.085)

# f(695,2500,0.085,0.825,0.09,0,0)
# f(695,30000,0.022,0.852,0.06,0.065,0)
# f(695,10000,0.0085,0.858,0.0535,0.08)
# f(695,15000,0.0031,0.86,0.051,0.085)