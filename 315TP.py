import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import timeit
import os
'''
Q- discharge
H - pump head
A - amps
V - voltage
RPM - speed
D - diameter
IP - input power (W)
FP - fluid power (W)
LP - lower pressure (Pa)
HP - higher pressure (Pa)
'''

waterDens = 1000
gravity = 9.81
SW = waterDens * gravity
visc = 1.0016 * 10**(-3)
LP = 310264.2
HP = 413685.6
InToCm = 0.0254
pipeL = 3200
pipeR = 0.0000015
pumpTestD = 4.5
eCost = 0.000077

def calcSingleParams(df):
    df['IP'] = df['V'] * df['A']
    df['FP'] = SW * df['Q'] * df['H']
    df['CQ'] = df['Q'] / (df['RPM'] * (df['D']*InToCm)**3)
    df['CH'] = df['H'] / ((df['RPM'])**2 * (df['D']*InToCm)**2)
    df['Eff'] = df['FP'] / df['IP']

    return df


def pumpScaling(scaledDF, testpump, diameter, motors):

    data = []
    diameter1 = str(diameter)
    diameter = diameter * InToCm
    for idx, row in testpump.iterrows():
        newQ = row['CQ']*row['RPM']*diameter**3
        newH = row['CH']*row['RPM']**2 * diameter**2
        newTEff = 1-(1-row['Eff']) * (row['D'] / diameter)**0.2
        if newQ != 0.0:
            newIEff = 0.95 - (row['Q'] / newQ)**0.32 * (0.95 - row['Eff'])
        else:
            newIEff = 0.0
        data.append([diameter, row['RPM'], newQ, newH, newTEff, newIEff])

    tmp_df = pd.DataFrame(data, columns=['pumpD', 'RPM', 'Q', 'H', 'TEff', 'IEff'])
    scaledDF = pd.concat([scaledDF, tmp_df], ignore_index=True)

    '''fig, ax = plt.subplots()
    RPMOps = [3450, 3202, 2898, 2501, 2001, 1598]
    for RPM in RPMOps:
        sysSeg = scaledDF[(scaledDF['RPM'] == RPM) & (scaledDF['pumpD'] == diameter)]
        ax.plot(sysSeg['Q'], sysSeg['H'], label=RPM)
    ax.set_title('Pump Diameter: ' + diameter1 + 'in')
    ax.legend(title='Pump Speed (RPM)')
    ax.set_xlabel('Q (m^3/s)')
    ax.set_ylabel('H (m)')
    plt.savefig('PumpCurve\\' + diameter1 + '.png')'''

    return scaledDF


def systemCurve(sysCurves, pipeD):
    tmp_df = pd.DataFrame(columns=['pipeD', 'a', 'c'])
    data = []
    pipeD = pipeD * InToCm
    for q in np.arange(0.005, 0.1, 0.005):
        Re = q / (np.pi * (pipeD / 2)**2) * waterDens * pipeD / visc
        fric = 0.0055 * (1 + (2 * 10**4 * pipeR / pipeD + 10**6 / Re)**(1/3))
        noleftoverH = 8 * q**2 / (gravity * np.pi**2 * pipeD**4) * (fric * pipeL / pipeD + 1)
        LH = noleftoverH + LP / SW
        UH = noleftoverH + HP / SW
        data.append([pipeD, q, LH, UH])

    tmp_df = pd.DataFrame(data, columns=['pipeD', 'Q', 'LH', 'UH'])
    sysCurves = pd.concat([sysCurves, tmp_df])
    return sysCurves


def findOPV(pipeD, pumpD, sysCurves, pumpCurves, speed, noPumps, graph):
    sysSub = sysCurves[sysCurves['pipeD'] == pipeD*InToCm]
    pumpSub = pumpCurves[(pumpCurves['pumpD'] == pumpD*InToCm)]

    sysFit = np.polyfit(sysSub['Q'], sysSub['LH'], 2)
    pumpFit = np.polyfit(pumpSub['Q'], pumpSub['H'], 2)

    def effCurve(x, a, b):
        return a * x**2 + b * x

    effFit, effcov = scipy.optimize.curve_fit(effCurve, pumpSub['Q'], pumpSub['TEff'])

    #effFit = np.polyfit(pumpSub['Q'], pumpSub['TEff'], 2)
    if graph:
        fig, ax1 = plt.subplots()

        xs = np.linspace(0, 0.3, 500)
        polysys = np.poly1d(sysFit)
        polypump = np.poly1d(pumpFit)
        #polyeff = np.poly1d(effFit)

        ln1 = ax1.plot(xs, polysys(xs), color='r', label='System Curve')
        #ax1.plot(sysSub['Q'], sysSub['LH'], 'o', color='r')
        ln2 = ax1.plot(xs, polypump(xs), color='b', label='Pump Curve')
        #ax1.plot(pumpSub['Q'], pumpSub['H'], 'o', color='b')
        ax1.annotate('Operating Point\nQ=0.185 Efficiency=0.54', xy=(0.185, 73.6), xytext=(0.15, 120), arrowprops={"arrowstyle":"->"})
        ax1.set_xlabel('Q (m3/s)')
        ax1.set_ylabel('H (m)')
        plt.xlim(0, 0.3)
        plt.ylim(0,160)

        ax2 = ax1.twinx()
        ln3 = ax2.plot(xs, effCurve(xs, *effFit) * 100, color='g', label='Efficiency Curve')
        ax2.set_ylabel('Efficiency (%)')
        plt.ylim(0,100)

        lns = ln1+ln2+ln3
        labs = [l.get_label() for l in lns]
        plt.legend(lns, labs, loc='upper left')

        plt.show()
        #plt.savefig('Graphs\\Pipe' + str(pipeD) + 'pump' + str(pumpD) + 'speed' + str(speed) + '.jpg')

    a = sysFit[0] - pumpFit[0]
    b = sysFit[1] - pumpFit[1]
    c = sysFit[2] - pumpFit[2]

    discriminant = b**2 - 4 * a * c
    if discriminant >= 0:
        OPQ1 = (-b + np.sqrt(discriminant)) / (2 * a)
        OPQ2 = (-b - np.sqrt(discriminant)) / (2 * a)
        if OPQ1 > 0:
            OPQ = OPQ1
        elif OPQ2 > 0:
            OPQ = OPQ2
        if OPQ1 > 0 or OPQ2 > 0:
            OPH = sysFit[0] * OPQ**2 + sysFit[1] * OPQ + sysFit[2]
            OPEff = effFit[0] * OPQ**2 + effFit[1] * OPQ
            if OPEff <= 0:
                tmp_df = pd.DataFrame([[speed, pipeD, pumpD, np.nan, np.nan, np.nan, noPumps]], columns=['RPM', 'Pipe', 'Pump', 'Q', 'H', 'TEff', 'NoPumps'])

            else:
                tmp_df = pd.DataFrame([[speed, pipeD, pumpD, OPQ, OPH, OPEff, noPumps]], columns=['RPM', 'Pipe', 'Pump', 'Q', 'H', 'TEff', 'NoPumps'])
            return tmp_df

        else:
            tmp_df = pd.DataFrame([[speed, pipeD, pumpD, np.nan, np.nan, np.nan, noPumps]], columns=['RPM', 'Pipe', 'Pump', 'Q', 'H', 'TEff', 'NoPumps'])
            #operatingPoints = pd.concat([operatingPoints, tmp_df])
            return tmp_df

    else:
        tmp_df = pd.DataFrame([[speed, pipeD, pumpD, np.nan, np.nan, np.nan, noPumps]], columns=['RPM', 'Pipe', 'Pump', 'Q', 'H', 'TEff', 'NoPumps'])
        #operatingPoints = pd.concat([operatingPoints, tmp_df])
        return tmp_df


def costing(pumpCurvesVar, sysCurves, pipeOps, pumpDiameters, waterDemandH, numPumps):

    electricityPumps = pd.DataFrame(columns=['pipeD', 'pumpD', 'Power24', 'CostYear', 'NoPumps'])

    # Cycle through pipe diameters system curves
    for pipeD in pipeOps:

        # Cycle through pump diameters
        for idx, pumpD in enumerate(pumpDiameters):
            print('Pipe' + str(pipeD) + 'PumpD' + str(pumpD))
            opPSub = pumpCurvesVar[(pumpCurvesVar['pumpD'] == pumpD*InToCm)].copy(deep=True)
            opPSub.sort_values(by='RPM', inplace=True)
            opPSub.reset_index(inplace=True, drop=True)
            tmp_df = pd.DataFrame(columns=['RPM', 'Pipe', 'Pump', 'Q', 'H', 'TEff', 'NoPumps'])
            for PNum in range(1, numPumps + 1):
                for n in range(1000,3451,10):
                    pumpCurves = opPSub.copy(deep=True)
                    pumpCurves['Q'] = PNum * pumpCurves['Q'] * (n / pumpCurves['RPM'])
                    pumpCurves['H'] = pumpCurves['H'] * (n / pumpCurves['RPM'])**2

                    tmp_df = pd.concat([tmp_df, findOPV(pipeD, pumpD, sysCurves, pumpCurves, n, PNum, False)])

                #tmp_df.to_csv('OPsvarN.csv')

                # Find power for each hour
                powerPerHour = []
                head = []
                hour_df = pd.DataFrame(columns=['PumpD', 'PipeD', 'Hour', 'TargetQ', 'Q', 'H', 'Eff', 'NPumps', 'RPM', 'Power'])
                invalid = False
                for idx, hourFlow in enumerate(waterDemandH):
                    tmp_df['dif'] = np.abs(tmp_df['Q'] - hourFlow)
                    tmp_df.sort_values(by='dif', inplace=True)
                    OP = tmp_df.iloc[0]
                    powerHour = SW * OP['Q'] * OP['H'] / OP['TEff']
                    if powerHour < 0.0:
                        powerHour = 0.0
                    powerPerHour.append(powerHour)
                    tmp_df3 = pd.DataFrame([[pumpD, pipeD, idx, hourFlow, OP['Q'], OP['H'], OP['TEff'], OP['NoPumps'], OP['RPM'], powerHour]],
                                           columns=['PumpD', 'PipeD', 'Hour', 'TargetQ', 'Q', 'H', 'Eff', 'NPumps', 'RPM', 'Power'])
                    hour_df = pd.concat([hour_df, tmp_df3])

                    # Graphing
                    if PNum == 2:
                        pumpCurves = opPSub.copy(deep=True)
                        pumpCurves['Q'] = OP['NoPumps'] * pumpCurves['Q'] * (OP['RPM'] / pumpCurves['RPM'])
                        pumpCurves['H'] = pumpCurves['H'] * (OP['RPM'] / pumpCurves['RPM'])**2
                        findOPV(pipeD, pumpD, sysCurves, pumpCurves, OP['RPM'], OP['NoPumps'], True)

                hour_df.to_csv('Hourly' + str(pipeD) + '-' + str(pumpD) + '.csv')
                print(np.sum(powerPerHour) * eCost * 365)
                if invalid:
                    totalPower24 = np.nan
                    cost = np.nan
                else:
                    totalPower24 = np.sum(powerPerHour)
                    cost = totalPower24 * eCost * 365

                tmp_df2 = pd.DataFrame([[pipeD, pumpD, totalPower24, cost, PNum]], columns=['pipeD', 'pumpD', 'Power24', 'CostYear', 'NoPumps'])
                electricityPumps = pd.concat([electricityPumps, tmp_df2])
                if totalPower24 != np.nan:
                    print('pipeD: ' + str(pipeD) + ' pumpD: ' + str(pumpD) + ' NoPumps: ' + str(PNum))
    return electricityPumps


def pumpOptimization():
    fixed = pd.read_csv('Pump2-Fixed.csv')
    variable = pd.read_csv('Pump1-Variable.csv')

    fixed = calcSingleParams(fixed)
    variable = calcSingleParams(variable)

    pumpDiameters = [4.5, 5, 5.5, 8.5, 9.25, 10, 10.25, 11, 12, 13.25, 14, 15.5, 16.75]
    pipeOps = [4, 6, 8,12,15,18,24,30,36,42,48]
    waterDemandH = [0.0083,0.0367,0.1850,0.1000,0.1050,0.0917,0.0900,0.0933,0.1100,0.1333,0.0900,
                    0.0600,0.0600,0.0667,0.0750,0.0350,0.0167,0.0150,0.0117,0.0083]
    maxNoMotors = 2

    pipeOps = [12,15,18,24,30,36]
    pumpDiameters = [8.5, 9.25, 10, 10.25, 11, 12, 13.25, 14, 15.5, 16.75]
    pumpDiameters = [8.5, 9.25, 10, 10.25, 11, 12, 13.25, 14, 15.5, 16.75]

    waterDemandH = [i * 1.15 for i in waterDemandH]
    pipeOps = [12]
    pumpDiameters = [15.5]

    #Pump Curves
    pumpCurvesFixed = pd.DataFrame(columns=['pumpD', 'RPM', 'Q', 'H', 'TEff', 'IEff'])
    pumpCurvesVar = pd.DataFrame(columns=['pumpD', 'RPM', 'Q', 'H', 'TEff', 'IEff'])
    for pumpD in pumpDiameters:
        pumpCurvesFixed = pumpScaling(pumpCurvesFixed, fixed, pumpD, maxNoMotors)
    fig, axs = plt.subplots(7, 2)
    plt.subplots_adjust(hspace=0.5)
    plt.suptitle("Variable Speed Pump Curves", fontsize=18, y=0.95)
    for idx, pumpD in enumerate(pumpDiameters):
        pumpCurvesVar = pumpScaling(pumpCurvesVar, variable, pumpD, maxNoMotors)

    #System Curves
    sysCurves = pd.DataFrame(columns=['pipeD', 'Q', 'LH', 'UH'])
    for pipeD in pipeOps:
        sysCurves = systemCurve(sysCurves, pipeD)

    pumpsOPs = costing(pumpCurvesVar, sysCurves, pipeOps, pumpDiameters, waterDemandH, maxNoMotors)

    #pumpsOPs.to_csv('OPV.csv')

if __name__ == '__main__':
    pumpOptimization()

