# Joel Surfleet's Isentropic Toolbox

import math as m

pi = m.pi

gamma = 1.4

class isentropic:
    def __init__(fluid,k = 1.4):
        fluid.k = k

    def assumptions():
        print("Isentropic Assumptions")

    def M(fluid,M):
        # T/T0 etc...
        fluid.Tratio = 1 / (1 + (((fluid.k - 1) / 2) * (M ** 2)))
        fluid.T(fluid.Tratio)
        
    def P(fluid,x):
        fluid.Pratio = x
        fluid.Dratio = x ** (1/fluid.k)
        fluid.Tratio = x ** ((fluid.k-1)/fluid.k)
        
    def D(fluid,x):
        fluid.Pratio = x ** (fluid.k)
        fluid.Dratio = x
        fluid.Tratio = x ** (fluid.k-1)
        
    def T(fluid,x):
        fluid.Pratio = x ** (fluid.k/(fluid.k-1))
        fluid.Dratio = x ** (1/(fluid.k-1))
        fluid.Tratio = x

    def getStag(fluid,P,D,T):
        fluid.Ps = P
        fluid.Ds = D
        fluid.Ts = T

        fluid.P0 = P / fluid.Pratio
        fluid.D0 = D / fluid.Dratio
        fluid.T0 = T / fluid.Tratio

        fluid.a = (fluid.k * P / D) ** (1/2)

    def getStatic(fluid,P0,D0,T0):
        fluid.P0 = P0
        fluid.D0 = D0
        fluid.T0 = T0

        fluid.Ps = P0 * fluid.Pratio
        fluid.Ds = D0 * fluid.Dratio
        fluid.Ts = T0 * fluid.Tratio

class quasi1D():
    def __init__(flow,k=1.4):
        flow.k = k

    def M(flow,M):
        flow.Aratio = 1/M * (2/(flow.k + 1) * (1 + ((flow.k - 1) / 2) * M**2)) ** ((flow.k + 1) / (2 * (flow.k - 1)))

    def getM(flow,x):
        def Af(M):
            return (1/M * (2/(flow.k + 1) * (1 + ((flow.k - 1) / 2) * M**2)) ** ((flow.k + 1) / (2 * (flow.k - 1)))) - x
        
        flow.Msup = bisection(1,25,Af)
        flow.Msub = bisection(0.001,1,Af)
    
    def getMwithNS(flow,Pratio,Aratio):
        flow.MNS = (-(1/(flow.k-1)) + (1/((flow.k-1)**2) + ((2/(flow.k-1))*(2/(flow.k+1))**((flow.k+1)/(flow.k-1)) * (1/(Pratio * Aratio))**2)) ** (1/2)) ** (1/2)

class rayleigh:
    def __init__(flow,k=1.4):
        flow.k = k

    def M(flow,M):
        M = M ** 2

        flow.Pratio = (1 + flow.k) / (1 + flow.k * M)
        flow.Tratio = M * ((1 + flow.k) / (1 + flow.k * M)) ** 2
        flow.Dratio = 1 / M * ((1 + flow.k * M) / (1 + flow.k))
        flow.T0ratio = (((flow.k + 1) * M) / ((1 + flow.k * M) ** 2)) * (2 + (flow.k - 1) * M)
        flow.P0ratio = ((1 + flow.k) / (1 + flow.k * M)) * (((2 + (flow.k - 1) * M) / (flow.k + 1)) ** (flow.k / (flow.k - 1)))

    def getStar(flow,P,D,T,P0,T0):
        flow.Pstar = P / flow.Pratio
        flow.Tstar = T / flow.Tratio
        flow.Dstar = D / flow.Dratio
        flow.P0star = P0 / flow.P0ratio
        flow.T0star = T0 / flow.T0ratio

        flow.P = P
        flow.D = D
        flow.T = T
        flow.P0 = P0
        flow.T0 = T0

    def getReal(flow,Pstar,Dstar,Tstar,P0star,T0star):
        flow.P = Pstar * flow.Pratio
        flow.D = Dstar * flow.Dratio
        flow.T = Tstar * flow.Tratio
        flow.P0 = P0star * flow.P0ratio
        flow.T0 = T0star * flow.T0ratio

        flow.Pstar = Pstar
        flow.Tstar = Tstar
        flow.Dstar = Dstar
        flow.P0star = P0star
        flow.T0star = T0star

class fanno:
    def __init__(flow,k = 1.4):
        flow.k = k

    def M(flow,M):
        M = M ** 2
        k = flow.k
        flow.fL4D = (1 - M) / (k * M) + (k + 1) /(2 * k) * m.log(((k + 1) * M) / (2 + (k - 1)*M))

        flow.Pratio = (1 / M ** (1/2)) * ((k + 1) / (2 + (k - 1) * M)) ** (1/2)
        flow.Dratio = (1 / M ** (1/2)) * ((2 + (k - 1) * M) / (k + 1)) ** (1/2)
        flow.Tratio = (k + 1) / (2 + (k - 1) * M)
        flow.P0ratio = (1 / M ** (1/2)) * ((2 + (k - 1) * M) / (k + 1)) ** ((k + 1) / (2 * (k-1)))

    def getStar(flow,P,D,T,P0):
        flow.Pstar = P / flow.Pratio
        flow.Tstar = T / flow.Tratio
        flow.Dstar = D / flow.Dratio
        flow.P0star = P0 / flow.P0ratio

        flow.P = P
        flow.D = D
        flow.T = T
        flow.P0 = P0

    def getReal(flow,Pstar,Dstar,Tstar,P0star):
        flow.P = Pstar * flow.Pratio
        flow.D = Dstar * flow.Dratio
        flow.T = Tstar * flow.Tratio
        flow.P0 = P0star * flow.P0ratio

        flow.Pstar = Pstar
        flow.Tstar = Tstar
        flow.Dstar = Dstar
        flow.P0star = P0star

class normal:
    def __init__(shock,k = 1.4):
        shock.k = k

    def M(shock,M):
        shock.Dratio = ((shock.k+1) * M ** 2) / (2 + (shock.k - 1) * M ** 2)
        shock.Pratio = 1 + (((2 * shock.k) / (shock.k + 1)) * ((M ** 2) - 1))
        shock.Tratio = shock.Pratio * shock.Dratio
        shock.M2 = ((1 + ((shock.k - 1) / 2) * M ** 2) / (shock.k * M ** 2 - (shock.k - 1) / 2)) ** (1/2)

    def get1(shock,P2,D2,T2):
        shock.P1 = P2 / shock.Pratio
        shock.T1 = T2 / shock.Tratio
        shock.D1 = D2 / shock.Dratio

        shock.P2 = P2
        shock.D2 = D2
        shock.T2 = T2

    def get2(shock,P1,D1,T1):
        shock.P2 = P1 * shock.Pratio
        shock.D2 = D1 * shock.Dratio
        shock.T2 = T1 * shock.Tratio

        shock.P1 = P1
        shock.T1 = T1
        shock.D1 = D1

class oblique():
    def __init__(shock,k = 1.4):
        shock.k = k

        shock.Normal = normal(shock.k)

    def getBeta(shock,theta,M):
        if theta < 0 or theta > 45 or M < 1 or M > 21:
            print("My chart is NOT good enough for that")
            pass

        f = open(r"C:\Users\joels\Documents\Python\ThetaMachBeta.txt",'r')

        dataString = f.read()
        dataString = dataString.splitlines()        
        f.close()

        for i in range(360):
            dataString[i] = dataString[i].split()[0:3]
            for j in range(3):
                dataString[i][j] = float(dataString[i][j])

        for j in range(0,360,40):
            if dataString[j][0] == theta:
                for k in range(0,40):
                    if dataString[j+k][1] > M:
                        if dataString[j+k][2] == '-1':
                            print('detached wave')
                            pass
                        shock.b = linInterp(M,dataString[j+k-1][1],dataString[j+k][1],dataString[j+k-1][2],dataString[j+k][2])
                        return shock.b
                    
            if dataString[j][0] > theta:
                lowerTheta = dataString[j-40][0]
                upperTheta = dataString[j][0]
                for k in range(0,40):
                    if dataString[j+k][1] > M:
                        if dataString[j+k][2] == '-1':
                            print('detached wave')
                            pass
                        lowerb = linInterp(M,dataString[j+k-1-40][1],dataString[j+k-40][1],dataString[j+k-1-40][2],dataString[j+k-40][2])
                        upperb = linInterp(M,dataString[j+k-1][1],dataString[j+k][1],dataString[j+k-1][2],dataString[j+k][2])
                        shock.b = linInterp(theta,lowerTheta,upperTheta,lowerb,upperb)
                        return shock.b
                    
    def getTheta(shock,M,beta):
        if beta < 0 or beta > 67 or M < 1 or M > 21:
            print("My chart is NOT good enough for that")
            pass

        f = open(r"C:\Users\joels\Documents\Python\ThetaMachBeta.txt",'r')

        dataString = f.read()
        dataString = dataString.splitlines()        
        f.close()

        for i in range(360):
            dataString[i] = dataString[i].split()[0:3]
            for j in range(3):
                dataString[i][j] = float(dataString[i][j])

        for j in range(0,40):
            if dataString[j][1] > M:
                lowerM = dataString[j-1][1]
                upperM = dataString[j][1]
                for k in range(0,360,40):
                    if dataString[k+j][2] >= beta:
                        uppertheta = linInterp(beta,dataString[j+k][2],dataString[j+k-40][2],dataString[j+k][0],dataString[j+k-40][0])
                        lowertheta = linInterp(beta,dataString[j+k-1][2],dataString[j+k-41][2],dataString[j+k-1][0],dataString[j+k-41][0])
                        # print(lowertheta)
                        # print(uppertheta)
                        shock.theta = linInterp(M,lowerM,upperM,lowertheta,uppertheta)
                        return shock.theta

    def getMach(shock,theta,beta):
        if theta < 0 or theta > 45 or beta < 0 or beta > 67:
            print("My chart is NOT good enough for that")
            pass

        f = open(r"C:\Users\joels\Documents\Python\ThetaMachBeta.txt",'r')

        dataString = f.read()
        dataString = dataString.splitlines()        
        f.close()

        for i in range(360):
            dataString[i] = dataString[i].split()[0:3]
            for j in range(3):
                dataString[i][j] = float(dataString[i][j])

        for j in range(0,360,40):
            if dataString[j][0] == theta:
                for k in range(0,40):
                    if dataString[j+k][2] < beta and dataString[j+k][2] != -1:
                        if dataString[j+k-1][2] < 0 or dataString[j+k-40][2] < 0 or dataString[j+k-41][2] < 0:
                            return 0
                        shock.Mach = linInterp(beta,dataString[j+k-1][2],dataString[j+k][2],dataString[j+k-1][1],dataString[j+k][1])
                        return shock.Mach
                    
            if dataString[j][0] > theta:
                lowerTheta = dataString[j-40][0]
                upperTheta = dataString[j][0]
                for k in range(0,40):
                    if dataString[j+k][2] < beta and dataString[j+k][2] != -1:
                        if dataString[j+k-1][2] < 0 or dataString[j+k-40][2] < 0 or dataString[j+k-41][2] < 0:
                            return 0
                        # print(dataString[j+k])
                        # print(dataString[j+k-1])
                        # print(dataString[j+k-40])
                        # print(dataString[j+k-41])
                        lowerM = linInterp(beta,dataString[j+k-40][2],dataString[j+k-1-40][2],dataString[j+k-1-40][1],dataString[j+k-40][1])
                        upperM = linInterp(beta,dataString[j+k][2],dataString[j+k-1][2],dataString[j+k-1][1],dataString[j+k][1])
                        # print(lowerM)
                        # print(upperM)
                        shock.Mach = linInterp(theta,lowerTheta,upperTheta,lowerM,upperM)
                        return shock.Mach
            
    def M(shock,M):
        shock.Mn = M * m.sin(m.radians(shock.b))
        shock.Normal.M(shock.Mn)

class expansion():
    def __init__(flow,k = 1.4):
        flow.k = k

    def M(flow,M):
        flow.v = (((flow.k+1)/(flow.k-1))**(1/2)) * m.degrees(m.atan(((flow.k-1)/(flow.k+1) * (M**2-1))**(1/2))) - m.degrees(m.atan((M**2-1)**(1/2)))
        flow.mu = m.degrees(m.asin(1/M))

    def getM(flow,x):
        def vf(M):
            return ((((flow.k+1)/(flow.k-1))**(1/2)) * m.degrees(m.atan(((flow.k-1)/(flow.k+1) * (M**2-1))**(1/2))) - m.degrees(m.atan((M**2-1)**(1/2)))) - x
        
        flow.M2 = bisection(1,25,vf)

def linInterp(x1,y1,y2,z1,z2):
    # x1: your data point
    # y1, y2: the values bracketing your data point
    # z1, z2: the values bracketing your answer
    return (x1 - y1) / (y2 - y1) * (z2 - z1) + z1

def bisection(a,b,f,TOL = 1e-8):

    if f(a) * f(b) > 0 or a > b:
        print("Bisection Method not possible")
        return 0

    for i in range(1000):

        c = (a + b) / 2

        if f(c) == 0 or (b - a) / 2 < TOL:
            return (a+b) / 2
        elif f(a) * f(c) > 0:
            a = c
        else:
            b = c

    print("No convergence")
    pass

def newtons(f, fprime, x0=0,TOL=1e-8,maxIterations=10000):

    for i in range(maxIterations):
        
        x1 = x0 - f(x0)/fprime(x0)
        
        if abs(x1 - x0) < TOL:
            return x1
        
        x0 = x1

    return None

def dataRead(path):

    # Input: A full file path
    # *** MAKE SURE FILE PATH IS A RAW STRING
    # File can only be either comma seperated or space seperated

    # Output: A dictionary with each data column in file stored as lists of float values.

    f = open(path)
    DATA = f.read()
    f.close()
    
    data = {
        'master' : DATA.splitlines()
    }

    sep = ' '
    if ',' in data['master'][1]:
        sep = ','
    elif '\t' in data['master'][1]:
        sep = '\t'
        
    for j in range(0,len(data['master'][1].split(sep))):
        data['{}'.format(j)] = [float(data['master'][i].split(sep)[j]) for i in range(1,len(data['master']))]

    data.pop('master')

    return data

def feasable(Isp,g0,Finert):
    dV = [Isp*g0*m.log(1/Finert[i]) for i in range(96)]
    return dV

def exitPressure(p1,Aratio,k):
    flow = quasi1D(k)
    flow.getM(Aratio)
    exhaust = isentropic(k)
    exhaust.M(flow.Msup)
    exhaust.getStatic(p1,1,1)
    return exhaust.Ps

def thrustCoeff(p1,p3,Aratio,k):
    # p1 = chamber
    # p2 = exhaust
    # p3 = ambient
    # Aratio = A2/At (exit to throat)
    # k = gamma
    p2 = exitPressure(p1,Aratio,k)
    Cf = ( (2*k**2)/(k-1) * (2/(k+1)) ** ((k+1)/(k-1)) * (1 - ((p2/p1) ** ((k-1)/k))) ) ** (1/2) + (p2-p3)/p1 * Aratio
    return Cf

def getMaxIndex(list1):
    curMax = 0
    for i in range(len(list1)):
        if list1[i] > curMax:
            curMax = list1[i]
            iMax = i

    return iMax

def characteristVelocity(k,R,T):
    cStar = (k*R*T) ** (1/2) / (k * ((2/(k+1)) ** ((k+1)/(k-1))) ** (1/2))
    return cStar