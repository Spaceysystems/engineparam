# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:51:24 2017

@author: Alison
"""
import math
import numpy

#Now to create a text file
title = 'Engine Parameter Calculator\n\nNomenclature'
#Variables 
V1 = '\nm = Mass Flow Rate (kg/s)'
V2 = '\nP = Pressure (Pa)'
V3 = '\np = density (kg/m^3)'
V4 = '\nv = velcity (m/s)'
V5 = '\nA = Area (m^3)'
V6 = '\nR = Gas Constant ()'
V7 = '\nW ='
V8 = '\nT = Temerpature (K)'
V9 = '\nM = mass'
V10 = '\nF = Thrust (N)'
V11 = '\nR = Gas Constant (J/kmol-K)'
V12 = '\nG = Acceleration due to gravity (m/s^2)'
End = '\n\nSubscripts denote chamber conditions (_c), throat conditions (_t) or exit conditions (_e)'
Intro = title+V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+End

#creates and opens a file
saveFile = open('ReadmeFile.txt','w')
#write the text in
saveFile.write(Intro)
#closes file
saveFile.close

#Assigning Global Variables for standard terrestrial conditions, to be overriden for non-standard conditions
g = 9.8 
P_a = 101,325 
k = 1.4 
R = 8,314.4621
p_a = 1.223

class engine:
###Overall Engine Paramters and Characteristics

#method one for Isp
    def Isp(g,It,mp):
        Isp = It/(g*mp)
        print(Isp,'s')
#method two for Isp
    def Isp2(EffectiveExVelC,g):
        Isp2 = EffectiveExVelC/g
        print(Isp2,'s')
#method three for Isp 
    def Isp3(A_e,g,m,P_e,P_a,v_e):
        Isp3 = (1/g)*((v_e + ((P_e-P_a)*A_e))/m)
        print(Isp3,'s')
#method four for Isp
    def Isp4(Cf,CharVel,g):
        Isp4=(CharVel*Cf)/g
        print(Isp4,'s')
#method 5 for Isp
    def Isp5(g,k,P_e,P_c,R,T_c,W):
        Isp5=(1/g)*(((2*k*R*T_c)/(W*(k-1)))*((1-(P_e/P_c)**(k - 1/k)))**(0.5))
        print(Isp5,'s')

#characteristic exhaust velocity (engine performance measure independant of nozzle performance)
    def cstar(A_t,m,P_c):
        cstar= (P_c*A_t)/m
        print(cstar,'m/s')
#characteristic exhaust velocity Method 2
    def cstar2(k,R,T_c,W):
        cstar2 = ((k*R*T_c/W)**.5)/k*(((2/(k+1))**((k+1)/(k-1)))**.5)
        print(cstar2,'m/s')
           
#deltaV
    def deltaV(EffectiveExVelC,Mf,Mo):
        deltaV = EffectiveExVelC*numpy.log(Mo/Mf)
        print(deltaV,'m/s')
        
#launchvehiclemassfraction
    def massfraction(Mo,Mf):
        m = 1 - (Mf/Mo)
        print(m)

#propulsive efficiency
    def efficiency(EffectiveExVel,VehicleVel):
        efficiency = ((2*VehicleVel)/EffectiveExVel)/(1+(VehicleVel/EffectiveExVel)**2)
        print(efficiency,'%')

#Thrust approximation
    def thrust(A_e,m,P_a,P_e,v_e):
        thrust = m*v_e + ((P_e-P_a)*A_e)
        print(thrust,'N')
#More Accurate Approximation
    def thrust2(A_e,A_t,P_a,P_c,k,P_e):
        thrust2 = A_t*P_c*(((2*k**2/(k-1))*((2/(k+1))**((k+1)/(k-1)))*(1-((P_e/P_c)**((k-1)/k))))**.5) + A_e*(P_e-P_a)
        print(thrust2,'N')
#Less Accurrate Approximation using Average Impulse and Thrust Time
    def thrust1(I_a,t_t):
        thrust1 = I_a/t_t
        print(thrust1,'N')

#massflow calculation 
    def massflow(burnArea,burnRate,p):
        massflow = burnArea*burnRate*p
        print(massflow,'kg/s')
#Massflow Calculation method 2
    def massflow2(A_t,P_c,k,R,T_c,W):
        massflow2 = (A_t*P_c*k/((k*R*T_c/W)**0.5))*(((2/(k+1))**((k+1)/(k-1)))**.5)
        print(massflow2,'kg/s')
#MassFlow Method three
    def massflow3(m_b,m_i,t_b):
        massflow3 = (m_i-m_b)/t_b
        print(massflow3,'kg/s')
        
#Burnrate
    def burnrate(BurnRateCoeff,BurnRateExponent,P_c):
        burnrate = BurnRateCoeff*(P_c**BurnRateExponent)
        print(burnrate)
            
# Fuel consumption
    def fueluse(burnrate,burntime):
        fueluse = burnrate*burntime
        print(fueluse,'kg/s')
        
#CF the Engine Thrust Coefficient 
    def Cf(A_e,A_t,k,P_a,P_c,P_e):
        Cf = ((((2*k**2)/(k-1))*(((2/(k+1)**(k+1)/(k-1))))*(1 - ((P_e/P_c)**((k-1)/1))))**.5) + (A_e/A_t)*((P_e-P_a)/P_c)
        print(Cf)
        
#Throat-Exit Area Ratio
    def AreaRatio(k,P_c,P_e):
        AreaRatio = (((k+1)/2)**(1/(k-1)))*((P_e/P_c)**(1/k))*(((k+1)/(k-1))*(1-((P_e/P_c)**((k-1)/k)))**.5)
        print(AreaRatio)

###Engine Chamber Characteristics###
#Chamber pressure
    def PressureC(p_c,R,T_c,W):
        PressureC = p_c*R*T_c/W
        print(PressureC,'Pa')      
#Throat Temperature
    def TempT(k,T_c):
        TempT = (2*T_c)/(k + 1)
        print(TempT,'K')
#flow rate throat
    def massflowthroat(A_t,p_t,v_t):
        massflowthroat = p_t*A_t*v_t
        print(massflowthroat,'kg/s')
#velcotiy throat
    def VelT(k,R,T_c,W):
        VelT = ((2*k*T_c*R/(W*(k +1)))**.5) 
        print(VelT,'m/s')
#throat vel method 2
    def VelT2(k,R,T_t,W):
        VelT2 = ((k*T_t*R/W)**.5)
        print(VelT2,'m/s')
#throat density
    def densityt(k,p_c):
        densityt = p_c*(2/(k+1)**(1/(k+1)))
        print(densityt,'kg/m^3')
#throat density approximation 
    def densityt2(k,p_c):
        densityt2 = p_c*((2/(k+1))**(1/(k+1)))
        print(densityt2,'kg/m^3') 
#throat Area
    def areaT(A_e,p_e,p_t,v_e,v_t):
        areaT = A_e*p_e*v_e/p_t*v_t
        print(areaT,'m^2')
#throat Area Method when exit area is unknown
    def areaT2(A_e,k,P_c,P_e):
        areaT2 = A_e*(((k+1)/2)**(1/(k-1)))*((P_e/P_c)**(1/k))*(((k+1)/(k-1))*(1-((P_e/P_c)**((k-1)/k)))**.5)
        print(areaT2,'m^2')
        
#massflow at throat
    def massflowt2(A_t,k,P_c,R,T_c,W):
        massflowt2 = (A_t*P_c*k/((k*R*T_c/W)**0.5))*((2/(k+1))**((k+1)/(k-1))**0.5)
        print(massflowt2,'kg/s')
        
#velocity at throat
    def VelT3(k,R,T_c,W):
        VelT3= (2*k*T_c*R/(W*(k+1)))**.5
        print(VelT3,'m/s')

###Condtions at Engine Nozzle Exit
#Exit Area 
    def exitarea(A_t,p_e,p_t,v_e,v_t):
        exitarea = A_t*p_t*v_t/p_e*v_e
        print(exitarea, 'm^2')
#Exit Area Method Two
    def exitarea2(A_t,k,P_c,P_e):
        exitarea2 = A_t/(((k+1)/2)**(1/(k-1)))*((P_e/P_c)**(1/k))*(((k+1)/(k-1))*(1-((P_e/P_c)**((k-1)/k)))**.5)
        print(exitarea2,'m^2')
#vel exit
    def ExhaustVel(k,P_c,P_e,R,T_c,W):
        ExhaustVel = ((R*2*k*T_c/W*(k-1))*(1-(P_e/P_c)**((k-1)/k)))**.5  
        print(ExhaustVel,'m/s')
#exit Area 
    def AreaE(A_t,p_e,p_t,v_e,v_t):
        AreaE = A_t*p_t*v_t/p_e*v_e
        print(AreaE,'m^2')
#effective exhaust velocity
    def Veff1(g,Isp):
        Veff1 = g*Isp
        print(Veff1,'m/s')
#method 3 effective exgaust vel
    def Veff(F,m):
        Veff = F/m
        print(Veff,'m/s')
#method 2 effective exhaust vel
    def Veff2(A_e,m,P_a,P_e,v_e):
        Veff2 = (v_e+((P_e-P_a)*A_e))/m
        print(Veff2,'m/s')
#atmosphericdrag
    def drag(A_s,Cd,p_a):
        drag = (p_a*Cd*A_s)/2
        print(drag,'kg/m')

#Creating an index to search for equations
searchindex = []
searchindex.append(engine.Isp.__code__.co_varnames)
searchindex.append(engine.Isp2.__code__.co_varnames)
searchindex.append(engine.Isp3.__code__.co_varnames)
searchindex.append(engine.Isp4.__code__.co_varnames)
searchindex.append(engine.Isp5.__code__.co_varnames)
searchindex.append(engine.cstar.__code__.co_varnames)
searchindex.append(engine.cstar2.__code__.co_varnames)
searchindex.append(engine.deltaV.__code__.co_varnames)
searchindex.append(engine.efficiency.__code__.co_varnames)
searchindex.append(engine.thrust.__code__.co_varnames)
searchindex.append(engine.thrust2.__code__.co_varnames)
searchindex.append(engine.massflow.__code__.co_varnames)
searchindex.append(engine.massflow2.__code__.co_varnames)
searchindex.append(engine.massflow3.__code__.co_varnames)
searchindex.append(engine.burnrate.__code__.co_varnames)
searchindex.append(engine.fueluse.__code__.co_varnames)
searchindex.append(engine.Cf.__code__.co_varnames)
searchindex.append(engine.AreaRatio.__code__.co_varnames)
searchindex.append(engine.PressureC.__code__.co_varnames)
searchindex.append(engine.TempT.__code__.co_varnames)
searchindex.append(engine.massflowthroat.__code__.co_varnames)
searchindex.append(engine.VelT.__code__.co_varnames)
searchindex.append(engine.densityt.__code__.co_varnames)
searchindex.append(engine.densityt2.__code__.co_varnames)
searchindex.append(engine.VelT2.__code__.co_varnames)
searchindex.append(engine.areaT.__code__.co_varnames)
searchindex.append(engine.areaT2.__code__.co_varnames)
searchindex.append(engine.massflowt2.__code__.co_varnames)
searchindex.append(engine.VelT3.__code__.co_varnames)
searchindex.append(engine.exitarea.__code__.co_varnames)
searchindex.append(engine.exitarea2.__code__.co_varnames)
searchindex.append(engine.ExhaustVel.__code__.co_varnames)
searchindex.append(engine.Veff.__code__.co_varnames)
searchindex.append(engine.Veff1.__code__.co_varnames)
searchindex.append(engine.Veff2.__code__.co_varnames)
searchindex.append(engine.AreaE.__code__.co_varnames)
searchindex.append(engine.drag.__code__.co_varnames)
print('Search for right function using search(Var1,Var1,etc) with variables inputted in alphabetical order \nAbbreviations are in the ReadmeFile')

##Creating a function to find relavent equations within the class
def search(*args):
    #x = 1
    for sublist in searchindex:
        #if args in sublist[:-1]:
            if sublist[:-1] == args:
                print('Here is a potentially useful equation -', sublist[-1])


##Unused Eq's might come in handy
#T0 = Tt + (Vt**.5/2*Cp) = Tt + ((M*(k*R*Tt/W)**.5)**2)/2*Cp
#h0 = he + (ve**2/2)
##CpTo = CpTe + (Ve**2/2)
#Ta/Tb = (Pa/Pb)**(k-1)/k) = (pa/pb)**(k-1)
#Tt/Tc = (pt/pc)**(k-1) = 2/(K+1)
#pt = pc*(2/(k+1))**(1/(k+1))
#(dA/A) = (M**2 - 1)*(dU/U)
#dP/p = -UdU
#gamma k = Cp/Cv
# Cp = (R/W)*(k/(k-1))
#ve**2 = 2*Cp*T0*(1-(Te/T0)) = 2*Cp*T0*(1-((Pe/P0)**((k-1)/k)))
#temperature sensitivity (pressure) of burn reaction def TempSens(): = (1/r)*(dR/dTb)_pc = (dlnR/dTb)_Pc
#K (not gamma) = Ab/At 
#temperature sensitivity of pressure Pi_k = (dlnPc/dTb)_k = (1/Pc)*(dPc/dTb)_k
#burn rate reaction = a*G_ox**n = a*(massflow/N*Pi*Rp**2)**n
