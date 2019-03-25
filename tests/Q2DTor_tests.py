#!/usr/bin/env python

"""
"""

import os

SOFTWARES = {}
SOFTWARES["gaussian"] = ("hf 3-21g","1.000") 
SOFTWARES["orca"    ] = ("hf 3-21g","1.000") 
#SOFTWARES["mpwb1k"  ] = ("mpwb95/6-31+G(d,p) IOp(3/76=0560004400) int=ultrafine","0.964") 


#-----------------------#
# Some useful functions #
#-----------------------#
def check_answer(answer,maxint,minint=1):

    try   :
       answer = int(answer)
       if (answer<minint) or (answer>maxint): return -1
    except: return -1
    return answer
#-----------------------#
def get_diff(list1,list2):
    return [ a-b for a,b in zip(list1,list2) ]
#-----------------------#
def get_norm(alist):
    norm = 0.0
    for comp in alist:
        norm += comp**2
    return norm**0.5
#-----------------------#


#----------------------------------------#
# Read/Write files                       #
#----------------------------------------#
def read_file(the_file):
    data  = open(the_file,'r')
    lines = data.readlines()
    data.close()
    return lines
#----------------------------------------#
def select_lines(lines,start,end):
    keep     = False
    selected = []
    for line in lines:
        if line.startswith(end)  :
           keep = False
           break
        if keep:
           selected.append(line)
        if line.startswith(start):
           keep = True
    return selected
#----------------------------------------#

#----------------------------------------#
# Get data from SXX.out files            #
#----------------------------------------#
def lines_between_args(lines,arg):
    target1 = "Executing Q2DTor with option: %s"
    target2 = "End of execution of Q2DTor with %s"

    # Get data in section: --pes
    idx1 = lines.rfind(target1%arg)
    idx2 = lines.rfind(target2%arg)
    if (idx1>idx2) or (idx1==-1) or (idx2==-1):
        return None
    return lines[idx1:idx2].split("\n")
#----------------------------------------#
def data_from_PES(lines):
    pes_minE = None
    pes_V00  = None
    DATA = lines_between_args(lines,"--pes")
    if DATA is None: return pes_minE, relE_00
    for line in DATA:
        if "Minimum energy (EMIN) in calculated points is" in line:
           pes_minE = float(line.split()[-2])
        if " 0.00 |   0.00 " in line:
           pes_V00 = float(line.split()[-1])
    return pes_minE, pes_V00
#----------------------------------------#
def data_from_FOURIER(lines):
    fourier_minE = None
    fourier_aveE = None
    fourier_maxE = None
    fourier_ctnt = None
    DATA = lines_between_args(lines,"--fourier")
    if DATA is None: return fourier_minE, fourier_aveE, fourier_maxE, fourier_ctnt
    for line in DATA:
        if "* Min. energy in potential is:" in line or ("* Minimum energy in potential is:" in line):
            fourier_minE = float(line.split()[-2])
        if "* Avg. energy in potential is:" in line or ("* Average energy in potential is:" in line):
            fourier_aveE = float(line.split()[-2])
        if "* Max. energy in potential is:" in line or ("* Maximum energy in potential is:" in line):
            fourier_maxE = float(line.split()[-2])
        if "V(phi1,phi2) = " in line:
            fourier_ctnt = float(line.split()[-2])
    return fourier_minE, fourier_aveE, fourier_maxE, fourier_ctnt
#----------------------------------------#
def data_from_OPTSP(lines,arg="--optxtr"):
    dict_SPs = {}
    DATA = lines_between_args(lines,arg)
    if DATA is None: return dict_SPs
    keep = False
    for idx in range(len( DATA )):
        line = DATA[idx]
        if keep and "Storing updated" in line:
           keep = False
           break
        if (":" in line) and ("|" in line) and keep:
           name, colon, phi1, xx, phi2, xx, V = line.split()
           dict_SPs[name] = [float(phi1),float(phi2),float(V)]
        if "Summary table (relative energy in cm":
           keep = True
    DATA = "\n".join(DATA)
    for name in dict_SPs.keys():
        idx = DATA.rfind("==> %s"%name)
        DATA2 = DATA[idx:].split("\n")
        for idx2 in range(len(DATA2)):
            line = DATA2[idx2]
            if "* Vibrational frequencies (cm^-1)" in line:
                freq_first = DATA2[idx2+1].split()[ 0]
                if "i" in freq_first: freq_first = -float(freq_first[:-1])
                else                : freq_first = +float(freq_first[:-1])
                dict_SPs[name] += [freq_first]
            if "* ZPE ="             in line:
                freq_last  = DATA2[idx2-2].split()[-1]
                if "i" in freq_last : freq_last  = -float(freq_last[:-1])
                else                : freq_last  = +float(freq_last[:-1])
                dict_SPs[name] += [freq_last]
                dict_SPs[name] += [float(line.split()[-5])]
            if "* reduced I1 ="      in line: dict_SPs[name] += [float(line.split()[-1])]
            if "* reduced I2 ="      in line: dict_SPs[name] += [float(line.split()[-1])]
            if "* lambda_12  ="      in line: dict_SPs[name] += [float(line.split()[-1])]
            if "* non-torsional ZPE" in line: dict_SPs[name] += [float(line.split()[-5])]; break
    return dict_SPs
#----------------------------------------#
def data_from_TOR2DNS(lines):
    tor2dns_evals = None
    DATA = lines_between_args(lines,"--tor2dns")
    if DATA is None: return tor2dns_evals
    for idx in range(len(DATA)):
        if "2DNS eigenvalues (cm^-1):" in DATA[idx]: break
    tor2dns_evals = [float(evalue) for evalue in  DATA[idx+1].split()]
    return tor2dns_evals
#----------------------------------------#
def data_from_ROVIBPF(lines):
    zpe_1who = None
    zpe_msho = None
    zpe_2dns = None
    zpe_ehr  = None
    zpe_e2dt = None
    rovibpf_100 = []
    rovibpf_700 = []
    rovibpf_2500 = []
    DATA = lines_between_args(lines,"--rovibpf")
    if DATA is None: return zpe_1who, zpe_msho, zpe_2dns, zpe_ehr, zpe_e2dt, rovibpf_100, rovibpf_700, rovibpf_2500

    for idx in range(len(DATA)):
        line = DATA[idx]
        if "* 1WHO =>" in line: zpe_1who = float(line.split()[3])
        if "* MSHO =>" in line: zpe_msho = float(line.split()[3])
        if "* 2DNS =>" in line: zpe_2dns = float(line.split()[3])
        if "* EHR  =>" in line: zpe_ehr  = float(line.split()[3])
        if "* E2DT =>" in line: zpe_e2dt = float(line.split()[3])
        if "  100.00 | " in line:
            rovibpf_100  += [float(data) for data in line.split("|")[1:]]
        if "  700.00 | " in line:
            rovibpf_700  += [float(data) for data in line.split("|")[1:]]
        if " 2500.00 | " in line:
            rovibpf_2500 += [float(data) for data in line.split("|")[1:]]

    # Ignore ratios
    if rovibpf_100  != []: rovibpf_100  = rovibpf_100 [ 0: 3]+rovibpf_100 [ 4:10]+rovibpf_100 [11:14]+rovibpf_100 [15:]
    if rovibpf_700  != []: rovibpf_700  = rovibpf_700 [ 0: 3]+rovibpf_700 [ 4:10]+rovibpf_700 [11:14]+rovibpf_700 [15:]
    if rovibpf_2500 != []: rovibpf_2500 = rovibpf_2500[ 0: 3]+rovibpf_2500[ 4:10]+rovibpf_2500[11:14]+rovibpf_2500[15:]

    return zpe_1who, zpe_msho, zpe_2dns, zpe_ehr, zpe_e2dt, rovibpf_100, rovibpf_700, rovibpf_2500
#----------------------------------------#


#----------------------------------------#
# String to Data ; Data to String        #
#----------------------------------------#
def data_to_string(name, pes_minE, pes_V00, fourier_minE, fourier_aveE, fourier_maxE,\
               fourier_ctnt, dict_SPs,tor2dns_evals, zpe_1who, zpe_msho, zpe_2dns,\
               zpe_ehr, zpe_e2dt,rovibpf_100, rovibpf_700, rovibpf_2500):
    string  = "start_%s\n"%name
    string += "  minE   %.6f\n"%pes_minE
    string += "  Vmin   %7.2f\n"%fourier_minE
    string += "  Vave   %7.2f\n"%fourier_aveE
    string += "  Vmax   %7.2f\n"%fourier_maxE
    string += "  V(0,0) %7.2f\n"%pes_V00
    string += "  const  %7.2f\n"%fourier_ctnt
    SPs  = sorted([SP for SP in dict_SPs.keys() if SP.startswith("MIN")])
    SPs += sorted([SP for SP in dict_SPs.keys() if SP.startswith("TS" )])
    SPs += sorted([SP for SP in dict_SPs.keys() if SP.startswith("MAX")])
    for SP in SPs:
        phi1, phi2, relE, freq_first, freq_last, ZPE, I1, I2, L12, ZPEnt = dict_SPs[SP]
        string += "  SP  %-5s  %6.2f  %6.2f  %8.2f  %+7.2f  %+7.2f  %8.2f  %5.3f  %5.3f  %+6.3f  %8.2f\n"%(SP,phi1, phi2, relE, freq_first, freq_last, ZPE, I1, I2, L12, ZPEnt)
    string += "  TOR2DNS  "+" ".join(["%7.2f"%evalue for evalue in tor2dns_evals])+"\n"
    string += "  zpe1WHO  %6.3f\n"%zpe_1who
    string += "  zpeMSHO  %6.3f\n"%zpe_msho
    string += "  zpe2DNS  %6.3f\n"%zpe_2dns
    string += "  zpeEHR   %6.3f\n"%zpe_ehr
    string += "  zpeE2DT  %6.3f\n"%zpe_e2dt
    for T,rovibpf in [("100",rovibpf_100),("700",rovibpf_700),("2500",rovibpf_2500)]:
       #string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[ 0: 3]+rovibpf[ 4: 7]]) + "\n"
       #string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[ 7:10]+rovibpf[11:14]]) + "\n"
       #string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[15:17]                    ]) + "\n"
       #string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[17:20]                    ]) + "\n"
       #string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[20:25]                    ]) + "\n"
       #string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[25:30]                    ]) + "\n"
       #string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[30:35]                    ]) + "\n"
        string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[ 0: 6]]) + "\n"
        string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[ 6:12]]) + "\n"
        string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[12:13]]) + "\n"
        string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[14:17]]) + "\n"
        string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[17:22]]) + "\n"
        string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[22:27]]) + "\n"
        string += "  PF%-4s   "%T+" ".join(["%9.3E"%pf for pf in rovibpf[27:32]]) + "\n"
    string += "end_%s\n"%name

    return string

def string_to_data(string):
    pes_minE = None
    pes_V00 = None
    fourier_minE = None
    fourier_aveE = None
    fourier_maxE = None
    fourier_ctnt = None
    dict_SPs = {}
    tor2dns_evals = None
    zpe_1who = None
    zpe_msho = None
    zpe_2dns = None
    zpe_ehr = None
    zpe_e2dt = None
    rovibpf_100 = []
    rovibpf_700 = []
    rovibpf_2500 = []
    for line in string.split("\n"):
        if   " minE "    in line: pes_minE = float(line.split()[1])
        elif " Vmin "    in line: fourier_minE = float(line.split()[1])
        elif " Vave "    in line: fourier_aveE = float(line.split()[1])
        elif " Vmax "    in line: fourier_maxE = float(line.split()[1])
        elif " const "   in line: fourier_ctnt = float(line.split()[1])
        elif " V(0,0) "  in line: pes_V00 = float(line.split()[1])
        elif " SP "      in line:
           sptype = line.split()[1]
           phi1, phi2, relE, freq_first, freq_last, ZPE, I1, I2, L12, ZPEnt = [float(xx) for xx in line.split()[2:]]
           dict_SPs[sptype] = phi1, phi2, relE, freq_first, freq_last, ZPE, I1, I2, L12, ZPEnt
        elif " TOR2DNS " in line: tor2dns_evals = [float(xx) for xx in line.split()[1:]]
        elif " zpe1WHO " in line: zpe_1who = float(line.split()[1])
        elif " zpeMSHO " in line: zpe_msho = float(line.split()[1])
        elif " zpeEHR "  in line: zpe_ehr  = float(line.split()[1])
        elif " zpe2DNS " in line: zpe_2dns = float(line.split()[1])
        elif " zpeE2DT " in line: zpe_e2dt = float(line.split()[1])
        elif " PF100"    in line: rovibpf_100  += [float(xx) for xx in line.split()[1:]]
        elif " PF700"    in line: rovibpf_700  += [float(xx) for xx in line.split()[1:]]
        elif " PF2500"   in line: rovibpf_2500 += [float(xx) for xx in line.split()[1:]]

    return pes_minE, pes_V00, fourier_minE, fourier_aveE, fourier_maxE,\
           fourier_ctnt, dict_SPs,tor2dns_evals, zpe_1who, zpe_msho, zpe_2dns,\
           zpe_ehr, zpe_e2dt,rovibpf_100, rovibpf_700, rovibpf_2500
#----------------------------------------#


#----------------------------------------#
# Get data for each system               #
#----------------------------------------#
def data_for_SXX(the_file,SXX,arg="--optxtr"):
    lines = "".join(read_file(the_file))
    pes_minE, pes_V00 = data_from_PES(lines)
    fourier_minE, fourier_aveE, fourier_maxE, fourier_ctnt = data_from_FOURIER(lines)
    dict_SPs = data_from_OPTSP(lines,arg)
    tor2dns_evals = data_from_TOR2DNS(lines)
    zpe_1who, zpe_msho, zpe_2dns, zpe_ehr, zpe_e2dt, rovibpf_100, rovibpf_700, rovibpf_2500 = data_from_ROVIBPF(lines)

    return pes_minE, pes_V00, fourier_minE, fourier_aveE, fourier_maxE,\
                            fourier_ctnt, dict_SPs,tor2dns_evals,zpe_1who, zpe_msho, zpe_2dns, zpe_ehr, \
                            zpe_e2dt,rovibpf_100,rovibpf_700,rovibpf_2500
#----------------------------------------#


#-------------------------------------#
# Generate diccionary with geometries #
#-------------------------------------#
def get_geometries():
    dict_SXX = {}

    string = ""
    string = string + "7 \n"
    string = string + "reference geometry for S01\n"
    string = string + "C    -0.74269   -0.75074   +0.00000\n"
    string = string + "O    -0.13161   -1.78126   +0.00000\n"
    string = string + "H    -1.83532   -0.68751   +0.00000\n"
    string = string + "C    +0.00000   +0.57565   +0.00000\n"
    string = string + "O    -0.58302   +1.61618   +0.00000\n"
    string = string + "O    +1.31054   +0.44447   +0.00000\n"
    string = string + "H    +1.52423   -0.49711   +0.00000\n"
    dict_SXX["S01"] = (string, "5-4-1-2", "7-6-4-1", 1, 1, "b", "min",0,1)

    string = ""
    string = string + "8 \n"
    string = string + "reference geometry for S02\n"
    string = string + "C    +0.69378   +0.56203   +0.01417\n"
    string = string + "C    +1.40341   -0.55173   +0.00906\n"
    string = string + "H    +1.11748   +1.55437   +0.03981\n"
    string = string + "H    +2.47643   -0.47848   +0.03766\n"
    string = string + "H    +0.94214   -1.52184   -0.03578\n"
    string = string + "O    -0.65351   +0.69125   -0.02124\n"
    string = string + "O    -1.26541   -0.57100   -0.09208\n"
    string = string + "H    -1.76783   -0.57788   +0.72547\n"
    dict_SXX["S02"] = (string, "7-6-1-2", "8-7-6-1", 1, 1, "b", "min",0,1)
    
    string = ""
    string = string + "8 \n"
    string = string + "reference geometry for S03\n"
    string = string + "C    +0.81848   +0.48847   +0.00000\n"
    string = string + "O    +1.33518   -0.59628   +0.00000\n"
    string = string + "H    +1.41883   +1.40896   -0.00000\n"
    string = string + "C    -0.66716   +0.64244   +0.00000\n"
    string = string + "O    -1.31127   -0.57958   -0.00000\n"
    string = string + "H    -0.94337   +1.23619   +0.87712\n"
    string = string + "H    -0.94338   +1.23619   -0.87712\n"
    string = string + "H    -0.63136   -1.25985   +0.00000\n"
    dict_SXX["S03"] = (string, "5-4-1-2", "8-5-4-1", 1, 1, "b", "min",0,1)
    
    string = ""
    string = string + "9 \n"
    string = string + "reference geometry for S04\n"
    string = string + "C    +1.20990   -0.21967   -0.00000\n"
    string = string + "C    -0.08551   +0.54320   -0.00000\n"
    string = string + "O    -1.13699   -0.39168   +0.00000\n"
    string = string + "H    +1.27507   -0.85339   -0.88034\n"
    string = string + "H    +2.05647   +0.46265   -0.00010\n"
    string = string + "H    +1.27518   -0.85322   +0.88046\n"
    string = string + "H    -0.14160   +1.18547   +0.88227\n"
    string = string + "H    -0.14160   +1.18547   -0.88228\n"
    string = string + "H    -1.97389   +0.06530   +0.00001\n"
    dict_SXX["S04"] = (string,"9-3-2-1","4-1-2-3",1,3,"bc","min",0,1)
    
    string = ""
    string = string + "9 \n"
    string = string + "reference geometry for S05\n"
    string = string + "C    +0.00000   +0.48283   +0.00000\n"
    string = string + "O    +1.15009   +0.83192   +0.00000\n"
    string = string + "O    -1.01376   +1.35186   +0.00000\n"
    string = string + "C    -0.47453   -0.91089   +0.00000\n"
    string = string + "C    +0.39897   -1.90612   +0.00000\n"
    string = string + "H    -1.54133   -1.06868   +0.00000\n"
    string = string + "H    +0.08080   -2.93627   +0.00000\n"
    string = string + "H    +1.45791   -1.69471   +0.00000\n"
    string = string + "H    -0.63457   +2.23447   +0.00000\n"
    dict_SXX["S05"] = (string,"5-4-1-2","9-3-1-2",1,1,"b", "min",0,1)
    
    string = ""
    string = string + "9 \n"
    string = string + "reference geometry for S06\n"
    string = string + "C    +1.22844   +0.35901   +0.00000\n"
    string = string + "O    +1.28140   -0.86759   +0.00000\n"
    string = string + "H    +2.16345   +0.93079   +0.00000\n"
    string = string + "C    +0.00000   +1.09210   +0.00000\n"
    string = string + "C    -1.17537   +0.41901   +0.00000\n"
    string = string + "H    +0.00533   +2.16782   +0.00000\n"
    string = string + "H    -2.12986   +0.92647   +0.00000\n"
    string = string + "O    -1.28040   -0.88226   +0.00000\n"
    string = string + "H    -0.36537   -1.24700   +0.00000\n"
    dict_SXX["S06"] = (string,"5-4-1-2","9-8-5-4",1,1,"b", "min",0,1)
    
    string = ""
    string = string + "10 \n"
    string = string + "reference geometry for S07\n"
    string = string + "C    +0.59880   +1.38351   +0.00001\n"
    string = string + "C    +0.51502   +0.06302   +0.00000\n"
    string = string + "H    -0.30050   +1.97275   +0.00001\n"
    string = string + "C    -0.75821   -0.71073   +0.00002\n"
    string = string + "O    +1.55120   -0.80441   -0.00001\n"
    string = string + "H    +1.55090   +1.89133   -0.00000\n"
    string = string + "F    -1.83911   +0.12990   -0.00002\n"
    string = string + "H    -0.81081   -1.34194   -0.88505\n"
    string = string + "H    -0.81082   -1.34189   +0.88512\n"
    string = string + "H    +2.37987   -0.32889   -0.00000\n"
    dict_SXX["S07"] = (string,"10-5-2-1","7-4-2-1",1,1,"b", "min",0,1)
    
    string = ""
    string = string + "10 \n"
    string = string + "reference geometry for S08\n"
    string = string + "C    +0.73137   -0.56534   +0.00001\n"
    string = string + "C    +1.42689   +0.56583   -0.00001\n"
    string = string + "H    +1.21087   -1.53193   +0.00002\n"
    string = string + "H    +2.50170   +0.51211   +0.00003\n"
    string = string + "H    +0.97415   +1.54288   -0.00003\n"
    string = string + "O    -0.60093   -0.70747   -0.00002\n"
    string = string + "C    -1.35663   +0.47115   +0.00001\n"
    string = string + "H    -1.14577   +1.06727   -0.88801\n"
    string = string + "H    -2.39755   +0.17229   +0.00001\n"
    string = string + "H    -1.14577   +1.06723   +0.88805\n"
    dict_SXX["S08"] = (string,"7-6-1-2","8-7-6-1",1,3,"bc", "min",0,1)
    
    string = ""
    string = string + "10 \n"
    string = string + "reference geometry for S09\n"
    string = string + "C    +1.49085   -0.51073   -0.03344\n"
    string = string + "C    +0.75861   +0.58964   +0.03955\n"
    string = string + "H    +1.03468   -1.47833   -0.18135\n"
    string = string + "C    -0.73176   +0.62497   -0.03268\n"
    string = string + "H    +1.23939   +1.55203   +0.16457\n"
    string = string + "H    +2.56659   -0.47405   +0.03244\n"
    string = string + "O    -1.33550   -0.63466   -0.05576\n"
    string = string + "H    -1.11594   +1.23085   +0.79380\n"
    string = string + "H    -1.04226   +1.12461   -0.94924\n"
    string = string + "H    -1.10467   -1.10113   +0.74528\n"
    dict_SXX["S09"] = (string,"7-4-2-1","10-7-4-2",1,1,"b", "min",0,1)
    
    string = ""
    string = string + "10 \n"
    string = string + "reference geometry for S10\n"
    string = string + "C    -0.98540   -1.02200   +0.00000\n"
    string = string + "C    -0.08254   -0.04719   -0.00000\n"
    string = string + "H    -0.67379   -2.05159   -0.00000\n"
    string = string + "C    +1.38890   -0.23851   +0.00000\n"
    string = string + "O    -0.39480   +1.26991   -0.00000\n"
    string = string + "H    -2.04578   -0.81814   +0.00000\n"
    string = string + "H    +1.64414   -1.29178   -0.00005\n"
    string = string + "H    +1.82612   +0.23427   -0.87669\n"
    string = string + "H    +1.82610   +0.23418   +0.87674\n"
    string = string + "H    -1.34414   +1.38001   +0.00001\n"
    dict_SXX["S10"] = (string,"10-5-2-1","7-4-2-1",1,3,"bc", "min",0,1)
    
    string = ""
    string = string + "10 \n"
    string = string + "reference geometry for S11\n"
    string = string + "C    -0.77079   +0.57053   -0.00000\n"
    string = string + "C    +0.55152   +0.68508   +0.00000\n"
    string = string + "H    -1.41622   +1.43379   -0.00000\n"
    string = string + "H    +0.94210   +1.69108   +0.00001\n"
    string = string + "C    +1.53777   -0.43195   +0.00000\n"
    string = string + "O    -1.49514   -0.56814   -0.00000\n"
    string = string + "H    +1.07424   -1.41738   -0.00020\n"
    string = string + "H    +2.18314   -0.38810   -0.87589\n"
    string = string + "H    +2.18289   -0.38834   +0.87608\n"
    string = string + "H    -0.91597   -1.32794   +0.00002\n"
    dict_SXX["S11"] = (string,"10-6-1-2","7-5-2-1",1,3,"bc", "min",0,1)
    
    string = ""
    string = string + "11 \n"
    string = string + "reference geometry for S12\n"
    string = string + "C    +1.49776   -0.93176   +0.00000\n"
    string = string + "C    +0.48303   -0.07088   -0.00000\n"
    string = string + "H    +2.52681   -0.60635   -0.00001\n"
    string = string + "C    +0.60493   +1.41142   +0.00000\n"
    string = string + "C    -0.87779   -0.63284   +0.00000\n"
    string = string + "H    +1.31909   -1.99783   +0.00000\n"
    string = string + "H    +1.64461   +1.72215   -0.00008\n"
    string = string + "H    +0.10750   +1.83199   -0.87089\n"
    string = string + "H    +0.10764   +1.83198   +0.87098\n"
    string = string + "O    -1.87733   +0.03706   -0.00000\n"
    string = string + "H    -0.93460   -1.73394   +0.00001\n"
    dict_SXX["S12"] = (string,"10-5-2-1","7-4-2-1",1,3,"bc", "min",0,1)
    
    string = ""
    string = string + "11 \n"
    string = string + "reference geometry for S13\n"
    string = string + "C    -1.30267   +0.54201   +0.00001\n"
    string = string + "C    +0.02113   +0.71518   -0.00001\n"
    string = string + "H    -1.98257   +1.37817   +0.00002\n"
    string = string + "H    +0.37510   +1.73459   -0.00001\n"
    string = string + "C    +1.01115   -0.33762   +0.00000\n"
    string = string + "O    -1.97468   -0.62075   -0.00001\n"
    string = string + "C    +2.32789   -0.13857   +0.00000\n"
    string = string + "H    +0.66514   -1.36563   +0.00001\n"
    string = string + "H    +3.02211   -0.96306   +0.00003\n"
    string = string + "H    +2.74449   +0.85807   -0.00002\n"
    string = string + "H    -1.37183   -1.36216   -0.00001\n"
    dict_SXX["S13"] = (string,"7-5-2-1","11-6-1-2",1,1,"b", "min",0,1)
    
    string = ""
    string = string + "11 \n"
    string = string + "reference geometry for S14\n"
    string = string + "C    -1.64872   -0.73920   -0.00000\n"
    string = string + "C    -0.51013   -0.04517   +0.00000\n"
    string = string + "H    -2.61485   -0.25839   -0.00001\n"
    string = string + "O    -0.44957   +1.30456   -0.00001\n"
    string = string + "C    +0.80035   -0.67805   +0.00000\n"
    string = string + "H    -1.62450   -1.81486   +0.00001\n"
    string = string + "C    +1.94605   -0.00775   +0.00000\n"
    string = string + "H    +0.79157   -1.75790   -0.00001\n"
    string = string + "H    +2.89211   -0.52481   -0.00000\n"
    string = string + "H    +1.95883   +1.07066   +0.00000\n"
    string = string + "H    -1.33180   +1.66990   +0.00005\n"
    dict_SXX["S14"] = (string,"7-5-2-1","11-4-2-1",1,1,"b", "min",0,1)
    
    string = ""
    string = string + "11 \n"
    string = string + "reference geometry for S15\n"
    string = string + "C    -0.16354   +0.75126   -0.00000\n"
    string = string + "C    +1.16625   +0.63542   -0.00000\n"
    string = string + "H    -0.62164   +1.72965   +0.00000\n"
    string = string + "H    +1.73680   +1.55606   +0.00000\n"
    string = string + "C    +1.98798   -0.60304   +0.00000\n"
    string = string + "C    -1.10604   -0.36330   -0.00000\n"
    string = string + "H    +1.40605   -1.51672   -0.00005\n"
    string = string + "H    +2.63936   -0.61446   -0.87253\n"
    string = string + "H    +2.63927   -0.61451   +0.87260\n"
    string = string + "O    -2.30267   -0.21036   +0.00000\n"
    string = string + "H    -0.68633   -1.37927   -0.00001\n"
    dict_SXX["S15"] = (string,"10-6-1-2","7-5-2-1",1,3,"bc", "min",0,1)
    
    string = ""
    string = string + "12 \n"
    string = string + "reference geometry for S16\n"
    string = string + "C    -0.71130   -0.30091   +0.32927\n"
    string = string + "C    -1.84687   +0.02149   -0.27194\n"
    string = string + "H    -0.65661   -1.23417   +0.87927\n"
    string = string + "H    -2.71531   -0.61713   -0.22255\n"
    string = string + "H    -1.94464   +0.94019   -0.83357\n"
    string = string + "C    +0.53537   +0.51557   +0.30561\n"
    string = string + "C    +1.71186   -0.24042   -0.29061\n"
    string = string + "H    +0.35709   +1.43338   -0.25310\n"
    string = string + "H    +0.78554   +0.81359   +1.32542\n"
    string = string + "H    +1.51263   -0.51593   -1.32349\n"
    string = string + "H    +2.61698   +0.36188   -0.26927\n"
    string = string + "H    +1.90990   -1.15616   +0.26334\n"
    dict_SXX["S16"] = (string,"7-6-1-2","10-7-6-1",1,3,"bc", "min",0,1)
    
    string = ""
    string = string + "12 \n"
    string = string + "reference geometry for S17\n"
    string = string + "C    +0.66455   +0.65999   -0.00000\n"
    string = string + "C    -0.66455   +0.65999   +0.00000\n"
    string = string + "H    +1.15839   +1.62378   +0.00000\n"
    string = string + "H    -1.15839   +1.62378   +0.00000\n"
    string = string + "C    -1.57474   -0.51894   +0.00000\n"
    string = string + "C    +1.57474   -0.51894   +0.00000\n"
    string = string + "H    -1.04150   -1.46393   +0.00011\n"
    string = string + "H    -2.22496   -0.50301   +0.87388\n"
    string = string + "H    -2.22481   -0.50313   -0.87399\n"
    string = string + "H    +1.04150   -1.46393   -0.00011\n"
    string = string + "H    +2.22496   -0.50301   -0.87388\n"
    string = string + "H    +2.22481   -0.50313   +0.87399\n"
    dict_SXX["S17"] = (string,"7-5-2-1","10-6-1-2",3,3,"abc", "min",0,1)
    
    string = ""
    string = string + "12 \n"
    string = string + "reference geometry for S18\n"
    string = string + "C    +1.82681   -0.42387   +0.00000\n"
    string = string + "C    +2.98307   +0.23549   -0.00000\n"
    string = string + "H    +1.82189   -1.50549   +0.00000\n"
    string = string + "H    +3.92854   -0.28262   -0.00000\n"
    string = string + "H    +3.01369   +1.31549   -0.00000\n"
    string = string + "C    +0.54298   +0.23700   +0.00000\n"
    string = string + "C    -0.63834   -0.39022   +0.00000\n"
    string = string + "H    +0.55469   +1.32171   +0.00000\n"
    string = string + "C    -1.88129   +0.36873   +0.00000\n"
    string = string + "H    -0.71785   -1.46805   +0.00000\n"
    string = string + "O    -2.97997   -0.12621   -0.00000\n"
    string = string + "H    -1.76050   +1.46585   +0.00000\n"
    dict_SXX["S18"] = (string,"7-6-1-2","11-9-7-6",1,1,"b", "min",0,1)
    
    string = ""
    string = string + "16 \n"
    string = string + "reference geometry for S19\n"
    string = string + "C    -0.44444   +0.24701   +0.13754\n"
    string = string + "C    +0.46516   +1.28342   -0.01406\n"
    string = string + "C    +1.82115   +1.02336   -0.12141\n"
    string = string + "C    +2.28065   -0.28193   -0.08887\n"
    string = string + "C    +1.37751   -1.32311   +0.05344\n"
    string = string + "C    +0.02381   -1.06026   +0.16805\n"
    string = string + "H    +0.10950   +2.30369   -0.05152\n"
    string = string + "H    +2.51721   +1.83975   -0.23979\n"
    string = string + "H    +3.33624   -0.48770   -0.17833\n"
    string = string + "H    +1.72959   -2.34325   +0.07756\n"
    string = string + "H    -0.68458   -1.86738   +0.28042\n"
    string = string + "C    -1.91162   +0.52997   +0.29023\n"
    string = string + "O    -2.73031   -0.43978   -0.30627\n"
    string = string + "H    -2.18446   +0.52613   +1.34367\n"
    string = string + "H    -2.13344   +1.52798   -0.09450\n"
    string = string + "H    -2.52089   -0.49170   -1.23682\n"
    dict_SXX["S19"] = (string,"13-12-1-2","16-13-12-1",2,1,"bc","min",0,1)
    
    string = ""
    string = string + "16 \n"
    string = string + "reference geometry for S20\n"
    string = string + "C    +0.58045   -0.63849   -0.00036\n"
    string = string + "C    +0.41164   +0.74418   -0.00030\n"
    string = string + "C    -0.85087   +1.31081   +0.00166\n"
    string = string + "C    -1.97258   +0.49845   +0.00078\n"
    string = string + "C    -1.82949   -0.87613   -0.00146\n"
    string = string + "C    -0.55741   -1.42728   -0.00094\n"
    string = string + "O    +1.54120   +1.49218   -0.00096\n"
    string = string + "H    -0.95630   +2.38711   +0.00299\n"
    string = string + "H    -2.95417   +0.94658   +0.00181\n"
    string = string + "H    -2.69765   -1.51606   -0.00070\n"
    string = string + "H    -0.43961   -2.50142   -0.00009\n"
    string = string + "C    +1.95778   -1.21274   +0.00126\n"
    string = string + "H    +2.51704   -0.88597   +0.87559\n"
    string = string + "H    +1.92271   -2.29804   +0.00010\n"
    string = string + "H    +2.52010   -0.88355   -0.87016\n"
    string = string + "H    +1.32121   +2.42104   -0.00559\n"
    dict_SXX["S20"] = (string,"16-7-2-1","13-12-1-2",1,3,"bc","min",0,1)

   #string = ""
   #string = string + "10 \n"
   #string = string + "reference geometry for S21\n"
   #string = string + "C   -1.237523   -0.277621   -0.075859\n"
   #string = string + "C    0.076034    0.360923    0.302580\n"
   #string = string + "H   -2.072613    0.367993    0.215400\n"
   #string = string + "H   -1.280452   -0.440853   -1.157158\n"
   #string = string + "H   -1.355331   -1.248498    0.420537\n"
   #string = string + "H    0.119522    1.489713   -0.272361\n"
   #string = string + "H    0.158940    0.617952    1.368390\n"
   #string = string + "O    1.154393   -0.450795   -0.161076\n"
   #string = string + "H    2.018133   -0.091247    0.122917\n"
   #string = string + "H    0.145591    2.411489   -0.769446\n"
   #dict_SXX["S21"] = (string,"9-8-2-1","3-1-2-8",1,3,"c", "ts",0,2)
     
   #string = ""
   #string = string + "8 \n"
   #string = string + "reference geometry for S22\n"
   #string = string + "C   -0.567877    0.612865    0.00000\n"
   #string = string + "H    0.568256    0.265135    0.00000\n"
   #string = string + "H   -0.654475    1.214970   -0.90104\n"
   #string = string + "H   -0.654471    1.214958    0.90104\n"
   #string = string + "O   -1.469907   -0.436822    0.00000\n"
   #string = string + "H   -1.038608   -1.285242   -0.00000\n"
   #string = string + "H    1.826068   -1.031844   -0.00000\n"
   #string = string + "O    1.889969   -0.070074   -0.00000\n"
   #dict_SXX["S22"] = (string,"6-5-1-2","5-1-8-7",1,1,"b", "ts",0,2)
    return dict_SXX
#-------------------------------------#
def get_defaults():
    '''
    Default values as would be read from input file
    This means that angles are in degrees
    '''
    #---------------#
    torsion1        = None
    torsion2        = None
    tsigma1         = None
    tsigma2         = None
    #---------------#
    level           = "hf sto-3g"
    charge          = "0"
    multiplicity    = "1"
    ttype           = "min"
    #---------------#
    t1step          = "10.0"
    t2step          = "10.0"
    symmetry        = "none"
    #---------------#
    weight          = "0.9"
    ignore          = "0.0"
    cos1            = "1-6"
    cos2            = "1-6"
    sin1            = "none"
    sin2            = "none"
    cos1cos2        = "1-6 , 1-6"
    sin1sin2        = "1-6 , 1-6"
    cos1sin2        = "none , none"
    sin1cos2        = "none , none"
    #---------------#
    tolerance       = "1.0"
    freqscal        = "1.000"
    #---------------#
    dijvar          = "yes"
    kmax            = "100"
    maxeigen        = "1e4"
    #---------------#
    interpolation   = "fourier"
    integrationstep = "1.0"
    #---------------#

    return torsion1     , torsion2       , tsigma1     , tsigma2,  \
           level        , charge         , multiplicity, ttype  ,  \
           t1step       , t2step         , symmetry    ,           \
           weight       , ignore         ,                         \
           cos1         , cos2           , cos1cos2    , sin1sin2, \
           sin1         , sin2           , cos1sin2    , sin1cos2, \
           tolerance    , freqscal       ,                         \
           dijvar       , kmax           , maxeigen    ,           \
           interpolation, integrationstep
#----------------------------------------------------------#
def string_inp(torsion1     , torsion2       , tsigma1     , tsigma2,  \
               level        , charge         , multiplicity, ttype  ,  \
               t1step       , t2step         , symmetry    ,           \
               weight       , ignore         ,                         \
               cos1         , cos2           , cos1cos2    , sin1sin2, \
               sin1         , sin2           , cos1sin2    , sin1cos2, \
               tolerance    , freqscal       ,                         \
               dijvar       , kmax           , maxeigen    ,           \
               interpolation, integrationstep,                         \
               Tlist                                                   ):

    string  = "#----------------------------------#\n"
    string += "# Torsional information            #\n"
    string += "#----------------------------------#\n"
    string += "start_torsions\n"
    string += "  torsion1         %-11s     # atoms involved in torsion 1\n"%torsion1
    string += "  torsion2         %-11s     # atoms involved in torsion 2\n"%torsion2
    string += "  tsigma1          %-3s             # torsional symmetry number of hindered rotor 1\n"%tsigma1
    string += "  tsigma2          %-3s             # torsional symmetry number of hindered rotor 2\n"%tsigma2
    string += "end_torsions\n"
    string += "#----------------------------------#\n"
    string += "# Calculations                     #\n"
    string += "#----------------------------------#\n"
    string += "start_calcs\n"
   #string += "  type             %-3s             # type of structure to analyze torsional anharmonicity (min,ts)\n"%ttype
    string += "  level            %-11s     # the calculation level\n"%level
    string += "  charge           %-2s              # charge\n"%charge
    string += "  multiplicity     %-2s              # spin multiplicity\n"%multiplicity
    string += "end_calcs\n"
    string += "#----------------------------------#\n"
    string += "# Torsional PES                    #\n"
    string += "#----------------------------------#\n"
    string += "start_pes\n"
    string += "  t1step           %-4s            # step in phi1 for scan calculation [degrees]\n"%t1step
    string += "  t2step           %-4s            # step in phi2 for scan calculation [degrees]\n"%t2step
    string += "  symmetry         %-4s            # Symmetry condition for PES: [a,b,c,ab,ac,bc,abc] or none\n"%symmetry
    string += "end_pes\n"
    string += "#----------------------------------#\n"
    string += "# Fitting details                  #\n"
    string += "#----------------------------------#\n"
    string += "start_fourier\n"
    string += "  weight           %-5s           #\n"%weight
    string += "  ignore           %-5s           # Set to zero coefficients with smaller absolute value (in cm^-1)\n"%ignore
    string += "  # Fourier Terms (Even)           #\n"
    string += "  cos1             %-4s            # i values in cos(i*Phi_1)\n"%cos1
    string += "  cos2             %-4s            # j values in cos(j*Phi_2)\n"%cos2
    string += "  cos1cos2         %-11s     # i,j values in cos(i*Phi_1) * cos(j*Phi_2)\n"%cos1cos2
    string += "  sin1sin2         %-11s     # i,j values in sin(i*Phi_1) * sin(j*Phi_2)\n"%sin1sin2
    string += "  # Fourier Terms (Odd)            #\n"
    string += "  sin1             %-4s            # i values in sin(i*Phi_1)\n"%sin1
    string += "  sin2             %-4s            # j values in sin(j*Phi_2)\n"%sin2
    string += "  cos1sin2         %-11s     # i,j values in cos(i*Phi_1) * sin(j*Phi_2)\n"%cos1sin2
    string += "  sin1cos2         %-11s     # i,j values in sin(i*Phi_1) * cos(j*Phi_2)\n"%sin1cos2
    string += "end_fourier\n"
    string += "#----------------------------------#\n"
    string += "# Search and Opt stationary points #\n"
    string += "#----------------------------------#\n"
    string += "start_statpoint\n"
    string += "  tolerance        %-5s           # step (in degrees) to explore torsional PES when looking for CPs\n"%tolerance
    string += "  freqscal         %-5s           # scaling factor for frequencies\n"%freqscal
    string += "end_statpoint\n"
    string += "#----------------------------------#\n"
    string += "# 2D-NS Hamiltonian                #\n"
    string += "#----------------------------------#\n"
    string += "start_tor2dns\n"
    string += "  dijvar           %-3s             # yes (dij not constant) or no (dij constant)\n"%dijvar
    string += "  kmax             %-3s             # check 2013-JChemPhys_138_134112, eq (14)\n"%kmax
    string += "  maxeigen         %-3s             # threshold for H eigenvalues (in cm^-1)\n"%maxeigen
    string += "end_tor2dns\n"
    string += "#----------------------------------#\n"
    string += "# Partition functions              #\n"
    string += "#----------------------------------#\n"
    string += "start_rovibpf\n"
    string += "  interpolation    %-7s         # fourier or spline order (1,3,5)\n"%interpolation
    string += "  integrationstep  %-7s         # integration dphi\n"%integrationstep
    string += "end_rovibpf\n"
    string += "#----------------------------------#\n"
    string += "# Working temperatures             #\n"
    string += "#----------------------------------#\n"
    string += "start_temperatures                 #\n"
    for idx in range(0,len(Tlist),3):
        line = "  ".join(["%6.1f"%T for T in Tlist[idx:idx+3]])
        line = "  " + line + " "*(33-len(line))+"#\n"
        string += line
    string += "end_temperatures                   #\n"
    string += "#----------------------------------#\n"
    string += "\n"
    return string
#-------------------------------------#



def check_rovibpf(data,ref,descrip):
       eps = 1.0
       if len(data) != len(ref):
          print "%sWARNING"%(descrip)
          print "Unable to find all data... Did you run Q2DTor with '--rovipf thermo'?"
       else:
          reldiff = [100.0*abs(a-b)/b for a,b in zip(data,ref)]
          if max(reldiff) > eps:
             print "%sWARNING"%(descrip)
             idx =  0
             if reldiff[idx] > eps: print "           rv(1WHO) : %.3E  vs  %.3E"%(data[idx],ref[idx])
             idx =  1
             if reldiff[idx] > eps: print "           rv(MSHO) : %.3E  vs  %.3E"%(data[idx],ref[idx])
             idx =  2
             if reldiff[idx] > eps: print "           rv(E2DT) : %.3E  vs  %.3E"%(data[idx],ref[idx])
             idx =  3
             if reldiff[idx] > eps: print "           Q_2DNS   : %.3E  vs  %.3E"%(data[idx],ref[idx])
             idx =  4
             if reldiff[idx] > eps: print "           Qtorclas : %.3E  vs  %.3E"%(data[idx],ref[idx])
             idx =  5
             if reldiff[idx] > eps: print "           Q_EHR    : %.3E  vs  %.3E"%(data[idx],ref[idx])

             idx = 12
             if reldiff[idx] > eps: print "           Q_trans  : %.3E  vs  %.3E"%(data[idx],ref[idx])
             idx = 13
             if reldiff[idx] > eps: print "           Q_ele    : %.3E  vs  %.3E"%(data[idx],ref[idx])

             idx = 17
             if reldiff[idx] > eps: print "           U^0(1WHO): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 18
             if reldiff[idx] > eps: print "           H^0(1WHO): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 19
             if reldiff[idx] > eps: print "           S^0(1WHO): %6.2f  vs  %6.2f  (cal/mol/K)"%(data[idx],ref[idx])
             idx = 20
             if reldiff[idx] > eps: print "           G^0(1WHO): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 21
             if reldiff[idx] > eps: print "           Cp (1WHO): %6.2f  vs  %6.2f  (cal/mol/K)"%(data[idx],ref[idx])

             idx = 22
             if reldiff[idx] > eps: print "           U^0(MSHO): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 23
             if reldiff[idx] > eps: print "           H^0(MSHO): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 24
             if reldiff[idx] > eps: print "           S^0(MSHO): %6.2f  vs  %6.2f  (cal/mol/K)"%(data[idx],ref[idx])
             idx = 25
             if reldiff[idx] > eps: print "           G^0(MSHO): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 26
             if reldiff[idx] > eps: print "           Cp (MSHO): %6.2f  vs  %6.2f  (cal/mol/K)"%(data[idx],ref[idx])

             idx = 27
             if reldiff[idx] > eps: print "           U^0(E2DT): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 28
             if reldiff[idx] > eps: print "           H^0(E2DT): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 29
             if reldiff[idx] > eps: print "           S^0(E2DT): %6.2f  vs  %6.2f  (cal/mol/K)"%(data[idx],ref[idx])
             idx = 30
             if reldiff[idx] > eps: print "           G^0(E2DT): %6.2f  vs  %6.2f  (kcal/mol)"%(data[idx],ref[idx])
             idx = 31
             if reldiff[idx] > eps: print "           Cp (E2DT): %6.2f  vs  %6.2f  (cal/mol/K)"%(data[idx],ref[idx])
          else:
             print "%sOK"%descrip


def check_data(data,ref,SXX):
    stringexe = "        Did you run Q2DTor with --%s?"

    descrip  = ["        Minimum energy in numerical PES ................."]
    descrip += ["        PES relative Energy at (phi1,phi2)=(0,0) ........"]
    descrip += ["        Minimum relative energy in numerical PES ........"]
    descrip += ["        Average relative energy in numerical PES ........"]
    descrip += ["        Maximum relative energy in numerical PES ........"]
    descrip += ["        Constant term in Fourier series ................."]
    descrip += ["        Information about stationary points ............."]
    descrip += ["        First eight 2D-NS eigenvalues ..................."]
    descrip += ["        Zero-Point Energy for 1WHO ......................"]
    descrip += ["        Zero-Point Energy for MSHO ......................"]
    descrip += ["        Zero-Point Energy for 2D-NS ....................."]
    descrip += ["        Zero-Point Energy for EHR ......................."]
    descrip += ["        Zero-Point Energy for E2DT ......................"]
    descrip += ["        Partition/Thermodynamic functions at  100K ......"]
    descrip += ["        Partition/Thermodynamic functions at  700K ......"]
    descrip += ["        Partition/Thermodynamic functions at 2500K ......"]

    #--------
    idx = 0
    if data[idx] is None: print "%sNot found...\n%s"%(descrip[idx],stringexe%"pes"); return
    else:
       diff = abs(data[idx]-ref[idx])*2625.5
       if diff > 1.0: print "%sdiffers more than 1 kj/mol, %.6f vs %.6f (hartree)"%(descrip[idx],data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 1
    if data[idx] is None: print "%sNot found...\n%s"%(descrip[idx],stringexe%"pes"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5: print "%sdiffers more than 0.5 cm^-1, %.6f vs %.6f (cm^-1)"%(descrip[idx],data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]
    #--------


    #--------
    idx = 2
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"fourier"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5: print "%sdiffers more than 0.5 cm^-1, %.6f vs %.6f (cm^-1)"%(descrip[idx],data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 3
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"fourier"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5: print "%sdiffers more than 0.5 cm^-1, %.6f vs %.6f (cm^-1)"%(descrip[idx],data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 4
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"fourier"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5: print "%sdiffers more than 0.5 cm^-1, %.6f vs %.6f (cm^-1)"%(descrip[idx],data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 5
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"fourier"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5: print "%sdiffers more than 0.5 cm^-1, %.6f vs %.6f (cm^-1)"%(descrip[idx],data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]
    #--------


    #--------
    idx = 6
    if data[idx] == {}  : print "%sNot found\n\n%s"%(descrip[idx],stringexe%"optsp"); return
    else:
       string = compare_sps(data[idx],ref[idx],SXX)
       if string == "": print "%sOK"%descrip[idx]
       else           : print "%sWARNING"%descrip[idx]; print string
    #--------


    #--------
    idx = 7
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"tor2dns"); return
    else:
       diff = get_diff(data[idx],ref[idx])
       norm = get_norm(diff)
       if norm > 0.5:
          print "%sWARNING"%(descrip[idx])
          for idx2 in range(len(ref[idx])):
              print "          %.2f vs %.2f"%(data[idx][idx2],ref[idx][idx2])
       else         : print "%sOK"%descrip[idx]
    #--------


    #--------
    idx = 8
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5:
          print "%sWARNING"%(descrip[idx])
          print "%.6f vs %.6f (cm^-1)"%(data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 9
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5:
          print "%sWARNING"%(descrip[idx])
          print "%.6f vs %.6f (cm^-1)"%(data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 10
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5:
          print "%sWARNING"%(descrip[idx])
          print "%.6f vs %.6f (cm^-1)"%(data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 11
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5:
          print "%sWARNING"%(descrip[idx])
          print "%.6f vs %.6f (cm^-1)"%(data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 12
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else:
       diff = abs(data[idx]-ref[idx])
       if diff > 0.5:
          print "%sWARNING"%(descrip[idx])
          print "%.6f vs %.6f (cm^-1)"%(data[idx],ref[idx])
       else         : print "%sOK"%descrip[idx]

    idx = 13
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else: check_rovibpf(data[idx],ref[idx],descrip[idx])

    idx = 14
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else: check_rovibpf(data[idx],ref[idx],descrip[idx])

    idx = 15
    if data[idx] is None: print "%sNot found\n\n%s"%(descrip[idx],stringexe%"rovibpf"); return
    else: check_rovibpf(data[idx],ref[idx],descrip[idx])
    #--------

def compare_sps(data,ref,SXX):
    not_in_ref  = list(set(data.keys()).difference( set( ref.keys()) ))
    not_in_data = list(set( ref.keys()).difference( set(data.keys()) ))

    # Correlate points #
    ref2data = {}
    for sp1 in data.keys():
           #phi1, phi2, relE, freq_first, freq_last, ZPE, I1, I2, L12, ZPEnt = data[sp1]
        for sp2 in ref.keys():
            if sp1[0:2] != sp2[0:2]: continue
            diff = (data[sp1][0]-ref[sp2][0])**2 + (data[sp1][1]-ref[sp2][1])**2
            if diff < 0.1**2:
               ref2data[sp2] = sp1

    # Differences #
    not_in_ref  = []
    not_in_data = []
    for sp in data.keys():
        if sp not in ref2data.values(): not_in_ref.append(sp)
    for sp in ref.keys():
        if sp not in ref2data.keys() : not_in_data.append(sp)
    string = ""
    for sp in not_in_ref:
       if "_" in sp: sptype = "TS"
       else        : sptype = sp[0:3]
       string += "             unexpected   : %3s at (%.1f,%.1f) \n"%(sptype,data[sp][0],data[sp][1])
    for sp in not_in_data:
       phi1 = int(round(ref[sp][0]))
       phi2 = int(round(ref[sp][1]))
       name = "%s_%03i_%03i"%(SXX,phi1,phi2)
       if "_" in sp: sptype = "TS"
       else        : sptype = sp[0:3]
       string += "             not presented  : %3s at (%.1f,%.1f) \n"%(sptype, ref[sp][0], ref[sp][1])
       if sptype == "MIN": string += "             line for splist: ' 0   %3i.00   %3i.00   -   NO   %s' \n"%(phi1,phi2,name)
       if sptype == "TS" : string += "             line for splist: ' 1   %3i.00   %3i.00   -   NO   %s' \n"%(phi1,phi2,name)
       if sptype == "MAX": string += "             line for splist: ' 2   %3i.00   %3i.00   -   NO   %s' \n"%(phi1,phi2,name)

    # Compare assigned #
    for sp2,sp1 in ref2data.items():

        diff = abs(data[sp1][2] - ref[sp2][2])
        if diff > 0.1: string += "             %s: relative energy differs (%7.2f vs %7.2f) \n"%(sp, data[sp1][2],ref[sp2][2])

        diff = abs(data[sp1][3] - ref[sp2][3])
        if diff > 1.0: string += "             %s: first frequency differs (%7.2f vs %7.2f) \n"%(sp, data[sp1][3],ref[sp2][3])

        diff = abs(data[sp1][4] - ref[sp2][4])
        if diff > 1.0: string += "             %s: last  frequency differs (%7.2f vs %7.2f) \n"%(sp, data[sp1][4],ref[sp2][4])
    
        diff = abs(data[sp1][5] - ref[sp2][5])
        if diff > 1.0: string += "             %s: ZPE differs significantly (%7.2f vs %7.2f) \n"%(sp, data[sp1][5],ref[sp2][5])
    
        diff = abs(data[sp1][6] - ref[sp2][6])
        if diff > 1.0: string += "             %s: first reduced moment of inertia differs (%7.2f vs %7.2f) \n"%(sp, data[sp1][6],ref[sp2][6])
    
        diff = abs(data[sp1][7] - ref[sp2][7])
        if diff > 1.0: string += "             %s: second reduced moment of inertia differs (%7.2f vs %7.2f) \n"%(sp, data[sp1][7],ref[sp2][7])
    
        diff = abs(data[sp1][8] - ref[sp2][8])
        if diff > 1.0: string += "             %s: coupling factor between torsions differs (%7.2f vs %7.2f) \n"%(sp, data[sp1][8],ref[sp2][8])
    
        diff = abs(data[sp1][9] - ref[sp2][9])
        if diff > 1.0: string += "             %s: non-torsional ZPE differs (%7.2f vs %7.2f) \n"%(sp, data[sp1][9],ref[sp2][9])
    
    return string
    


#----------------------------------------#
# OPTIONS of this script                 #
#----------------------------------------#
def option_createINPs():

    Tlist = "    100.0  150.0  200.0      "
    Tlist+= "    250.0  300.0  400.0      "
    Tlist+= "    500.0  700.0 1000.0      "
    Tlist+= "   1500.0 2000.0 2500.0      "
    Tlist = [float(T) for T in Tlist.split()]

    dict_SXX = get_geometries()

    for software in SOFTWARES.keys():
        folder = software.upper()+"/"

        print "    - Creating test files in: %s"%folder
        # Create folders
        if not os.path.exists(folder): os.mkdir(folder)

        # Go system by system
        for SXX in sorted(dict_SXX.keys()):

            print "      * %s"%SXX
            # Create SXX folder
            if os.path.exists(folder+SXX):
               print "        folder '%s' already exists..."%(folder+SXX)
               continue
            if not os.path.exists(folder+SXX): os.mkdir(folder+SXX)

            # Get defaults
            torsion1     , torsion2       , tsigma1     , tsigma2,  \
            level        , charge         , multiplicity, ttype  ,  \
            t1step       , t2step         , symmetry    ,           \
            weight       , ignore         ,                         \
            cos1         , cos2           , cos1cos2    , sin1sin2, \
            sin1         , sin2           , cos1sin2    , sin1cos2, \
            tolerance    , freqscal       ,                         \
            dijvar       , kmax           , maxeigen    ,           \
            interpolation, integrationstep                            = get_defaults()

            # Get data for SXX
            xyz_string, torsion1, torsion2, tsigma1, tsigma2, symmetry, ttype, charge, multiplicity = dict_SXX[SXX]
            charge  = str(charge)
            multiplicity = str(multiplicity)
            level, freqscal = SOFTWARES[software]

            if tsigma1 == 2:
               cos1 = "2 4 6 8"
               sin1 = "2 4 6 8"
               cos1cos2 = "2 4 6 , %s"
               cos1sin2 = "2 4 6 , %s"
               sin1cos2 = "2 4 6 , %s"
               sin1sin2 = "2 4 6 , %s"
            elif tsigma1 == 3:
               cos1 = "3 6 9"
               sin1 = "3 6 9"
               cos1cos2 = "3 6 , %s"
               cos1sin2 = "3 6 , %s"
               sin1cos2 = "3 6 , %s"
               sin1sin2 = "3 6 , %s"
            else:
               cos1 = "1-9"
               sin1 = "1-9"
               cos1cos2 = "1-7 , %s"
               cos1sin2 = "1-7 , %s"
               sin1cos2 = "1-7 , %s"
               sin1sin2 = "1-7 , %s"

            if tsigma2 == 2:
               cos2 = "2 4 6 8"
               sin2 = "2 4 6 8"
               cos1cos2 = cos1cos2%"2 4 6"
               cos1sin2 = cos1sin2%"2 4 6"
               sin1cos2 = sin1cos2%"2 4 6"
               sin1sin2 = sin1sin2%"2 4 6"
            elif tsigma2 == 3:
               cos2 = "3 6 9"
               sin2 = "3 6 9"
               cos1cos2 = cos1cos2%"3 6"
               cos1sin2 = cos1sin2%"3 6"
               sin1cos2 = sin1cos2%"3 6"
               sin1sin2 = sin1sin2%"3 6"
            else:
               cos2 = "1-9"
               sin2 = "1-9"
               cos1cos2 = cos1cos2%"1-7"
               cos1sin2 = cos1sin2%"1-7"
               sin1cos2 = sin1cos2%"1-7"
               sin1sin2 = sin1sin2%"1-7"


            if "b" in symmetry:
               sin1     = "none"
               sin2     = "none"
               cos1sin2 = "none"
               sin1cos2 = "none"

            # Get string
            inp_string = string_inp(torsion1     , torsion2       , tsigma1     , tsigma2,  \
                                    level        , charge         , multiplicity, ttype  ,  \
                                    t1step       , t2step         , symmetry    ,           \
                                    weight       , ignore         ,                         \
                                    cos1         , cos2           , cos1cos2    , sin1sin2, \
                                    sin1         , sin2           , cos1sin2    , sin1cos2, \
                                    tolerance    , freqscal       ,                         \
                                    dijvar       , kmax           , maxeigen    ,           \
                                    interpolation, integrationstep,                         \
                                    Tlist                                                   )

            # Write xyz file
            thefile = open(folder+SXX+"/"+SXX+".xyz",'w')
            thefile.write(xyz_string)
            thefile.close()

            # Write input file
            thefile = open(folder+SXX+"/"+SXX+".inp",'w')
            thefile.write(inp_string)
            thefile.close()
        print
#----------------------------------------#
def option_check():

    dict_SXX = get_geometries()

    #- - - - - - - - - - - #
    # Select what to check #
    #- - - - - - - - - - - #
    print "        Which folder do you want to check?"
    print "        (0) GAUSSIAN/ and ORCA/"
    print "        (1) only GAUSSIAN/"
    print "        (2) only ORCA/"
    answer = raw_input("        >> ")
    print
    answer = check_answer(answer,maxint=2,minint=0)
    if answer == -1: print "    Your choice is not valid"; exit()
    if answer == 0: folders = ["GAUSSIAN/","ORCA/"    ]
    if answer == 1: folders = ["GAUSSIAN/"]
    if answer == 2: folders = ["ORCA/"    ]

    print "        Which system do you want to check?"
    print "        ( 0) all the systems: S01 to S20"
    print "        ( 1) only S01 "
    print "        ( 2) only S02 "
    print "        ...         "
    print "        (19) only S19"
    print "        (20) only S20"
    answer = raw_input("        >> ")
    print
    answer = check_answer(answer,maxint=20,minint=0)
    if answer == -1: print "    Your choice is not valid"; exit()
    if answer == 0: list_SXX = sorted(dict_SXX.keys())
    else          : list_SXX = ["S%02i"%answer]


    #- - - - - - - - - - - - #
    # Check selected systems #
    #- - - - - - - - - - - - #
    for folder in folders:
        print "    Checking results in ./%s"%folder
        print "   ================================="
        print
        if not os.path.exists(folder):
           print "      folder was not found..."
           print
           print
           continue
        # Reference data
        if folder == "GAUSSIAN/": DATA = DATA_GAUSSIAN.split("\n")
        if folder == "ORCA/"    : DATA = DATA_ORCA.split("\n")
        if folder == "MPWB1K/"  : DATA = DATA_MPWB1K.split("\n")
        dict_REF = {}
        for SXX in list_SXX:
            start = "start_%s"%SXX
            end   = "end_%s"%SXX
            lines = select_lines(DATA,start,end)
            dict_REF[SXX] = string_to_data("\n".join(lines))

        # Get user data and check it
        for SXX in list_SXX:
            print "        ---------------"
            print "        - system: %3s -"%SXX
            print "        ---------------"
            out = folder+SXX+"/"+SXX+".out"
            if not os.path.exists(out):
               print "        file not found: %s"%out
               print
               continue
            dataSXX = data_for_SXX(out,SXX,arg="--optsp")
            check_data(dataSXX,dict_REF[SXX],SXX)
            print
        print
        print
#----------------------------------------#
def option_backup():
    dict_SXX = get_geometries()
    folders = ["SXX-hf_321g-Gaussian/","SXX-hf_321g-Orca/","SXX-MPWB1K-Gaussian/"]
    for folder in folders:
        folder = "/home/david/TEST/01-TESTS_Q2DTor/" + folder
        print folder
        BACKUP = "#======================================#\n"
        for SXX in sorted(dict_SXX.keys()):
            file_out = folder+SXX+"/"+SXX+".out"
            # string
            data = tuple([SXX] + list(data_for_SXX(file_out,SXX)))
            the_string = data_to_string(*data)
            BACKUP += the_string
            BACKUP += "#======================================#\n"
        # Save backup
        if folder == "SXX-hf_321g-Gaussian/": name = "GAUSSIAN"
        if folder == "SXX-hf_321g-Orca/"    : name = "ORCA"
        if folder == "SXX-MPWB1K-Gaussian/" : name = "MPWB1K"
        the_file = open("DATA_%s"%name,"w")
        the_file.write(BACKUP)
        the_file.close()
#----------------------------------------#


#===================#
# THE MAIN FUNCTION #
#===================#
def main():
    print
    print " << Test creator for Q2DTor >>"
    print
    print "    (1) Create input files"
    print "    (2) Check results"
    print
    answer = raw_input("    your choice: ")
    print
    # backup (not for users)
    if answer == "backup": option_backup(); exit()
    # Check answer
    answer = check_answer(answer,maxint=2,minint=1)
    # Act according to answer
    if answer == -1: print "    Your choice is not valid"; exit()
    if answer ==  1: option_createINPs()
    if answer ==  2: option_check()
#===================#





DATA_GAUSSIAN = """
#======================================#
start_S01
  minE   -299.779185
  Vmin      0.00
  Vave   3916.26
  Vmax   7301.80
  V(0,0) 3797.16
  const  3940.69
  SP  MIN01  180.00  180.00      0.00  +175.00  +3870.00  10182.93  6.336  0.729  -0.069   9792.73
  SP  MIN02    0.00  180.00    151.22  +162.00  +3877.00  10166.06  6.475  0.744  +0.382   9778.98
  SP  MIN03  180.00    0.00    304.63  +229.00  +3849.00  10227.38  7.010  0.747  -0.344   9794.05
  SP  MIN04    0.00    0.00   3797.13  +138.00  +3916.00   9957.72  6.544  0.675  -0.189   9712.95
  SP  TS_01   92.38  179.83   2592.00  -156.80  +3851.00   9936.48  7.700  0.746  +0.199   9634.63
  SP  TS_02  173.98   90.17   3577.22  -548.10  +3866.00   9746.50  6.541  0.728  -0.161   9638.45
  SP  TS_03    6.67   72.99   4574.95  -408.70  +3854.00   9697.55  6.530  0.719  +0.030   9604.69
  SP  TS_04   82.33   13.32   5557.98  -165.70  +3895.00   9775.74  7.666  0.690  -0.232   9581.99
  SP  MAX01   88.04   82.07   6396.37  -440.00  +3845.00   9482.74  7.798  0.734  +0.113   9463.50
  SP  MAX02  270.81   77.40   7307.63  -497.30  +3853.00   9467.93  7.734  0.728  -0.287   9444.67
  TOR2DNS   372.81  519.98  534.20  672.35  682.77  692.41  821.33  847.40
  zpe1WHO  29.114
  zpeMSHO  29.114
  zpe2DNS   1.066
  zpeEHR   27.999
  zpeE2DT  29.065
  PF100    1.475E+04 1.705E+04 1.793E+04 1.258E+00 5.640E-02 8.041E+02
  PF100    3.472E-60 4.013E-60 5.420E-60 5.889E-03 5.640E-02 5.191E-59
  PF100    1.652E+06 1.000E+00
  PF100    5.734E-54 6.629E-54 8.953E-54
  PF100    6.560E-01 8.540E-01 5.607E+01 -4.753E+00 9.848E+00
  PF100    7.110E-01 9.100E-01 5.692E+01 -4.781E+00 1.095E+01
  PF100    7.180E-01 9.170E-01 5.708E+01 -4.791E+00 1.113E+01
  PF700    1.105E+07 2.405E+07 3.127E+07 1.227E+01 6.022E+00 1.534E+07
  PF700    8.984E-03 1.956E-02 2.636E-02 5.705E+00 6.022E+00 2.783E-02
  PF700    2.141E+08 1.000E+00
  PF700    1.924E+06 4.188E+06 5.645E+06
  PF700    1.082E+01 1.221E+01 8.779E+01 -4.924E+01 2.593E+01
  PF700    1.112E+01 1.251E+01 8.977E+01 -5.033E+01 2.610E+01
  PF700    1.152E+01 1.291E+01 9.086E+01 -5.069E+01 2.727E+01
  PF2500   1.184E+12 3.382E+12 6.927E+12 1.552E+02 1.252E+02 5.590E+12
  PF2500   3.374E+09 9.638E+09 1.994E+10 1.252E+02 1.252E+02 1.994E+10
  PF2500   5.162E+09 1.000E+00
  PF2500   1.741E+19 4.975E+19 1.029E+20
  PF2500   6.598E+01 7.095E+01 1.281E+02 -2.492E+02 3.568E+01
  PF2500   6.735E+01 7.232E+01 1.307E+02 -2.544E+02 3.646E+01
  PF2500   6.839E+01 7.336E+01 1.325E+02 -2.580E+02 3.553E+01
end_S01
#======================================#
start_S02
  minE   -226.398436
  Vmin      0.00
  Vave   1881.82
  Vmax   4546.20
  V(0,0) 3716.41
  const  1847.42
  SP  MIN01    3.31  214.80      0.00  +128.00  +3818.00  14084.03  3.903  0.860  -0.226  13915.42
  SP  MIN02  149.01  156.24    468.58   +87.00  +3819.00  13987.44  7.221  0.857  +0.197  13857.77
  SP  TS_01    0.00  180.00     21.73   -93.70  +3826.00  14022.18  3.873  0.852  -0.196  13920.04
  SP  TS_02  180.00  180.00    578.10  -100.80  +3826.00  13860.75  7.053  0.849  +0.266  13838.36
  SP  TS_03   74.30  173.55   1372.59  -146.30  +3803.00  13916.12  6.566  0.865  +0.035  13834.52
  SP  TS_04    0.00    0.00   3724.05  -546.40  +3816.00  13941.08  3.679  0.804  -0.584  13887.11
  SP  TS_05  114.45    7.27   3808.32  -424.70  +3744.00  13862.36  7.394  0.804  -0.547  13820.92
  SP  MAX01  320.32   19.22   3815.79  -442.80  +3788.00  13852.30  4.815  0.814  -0.542  13832.17
  SP  MAX02   85.14    5.96   3864.51  -386.80  +3742.00  13814.75  7.052  0.822  -0.596  13810.10
  SP  MAX03  180.00    0.00   4553.88  -473.50  +3757.00  13786.03  7.031  0.758  -0.494  13782.44
  TOR2DNS   128.42  195.76  308.83  312.18  380.67  437.30  488.73  494.53
  zpe1WHO  40.268
  zpeMSHO  40.268
  zpe2DNS   0.367
  zpeEHR   39.786
  zpeE2DT  40.153
  PF100    8.621E+03 1.733E+04 1.085E+04 1.586E+00 3.658E-01 2.503E+03
  PF100    8.534E-85 1.716E-84 1.915E-84 2.499E-01 3.658E-01 2.804E-84
  PF100    1.207E+06 1.000E+00
  PF100    1.030E-78 2.070E-78 2.311E-78
  PF100    7.050E-01 9.040E-01 5.487E+01 -4.584E+00 1.085E+01
  PF100    7.110E-01 9.090E-01 5.632E+01 -4.722E+00 1.116E+01
  PF100    7.420E-01 9.410E-01 5.570E+01 -4.629E+00 1.127E+01
  PF700    1.079E+07 3.700E+07 2.978E+07 3.146E+01 2.433E+01 2.303E+07
  PF700    2.891E-06 9.912E-06 8.664E-06 2.416E+01 2.433E+01 8.724E-06
  PF700    1.564E+08 1.000E+00
  PF700    4.521E+02 1.550E+03 1.355E+03
  PF700    1.124E+01 1.263E+01 8.772E+01 -4.877E+01 2.740E+01
  PF700    1.175E+01 1.314E+01 9.090E+01 -5.049E+01 2.782E+01
  PF700    1.173E+01 1.312E+01 9.044E+01 -5.019E+01 2.717E+01
  PF2500   3.075E+12 1.459E+13 8.945E+12 3.309E+02 3.073E+02 8.308E+12
  PF2500   9.282E+08 4.405E+09 2.764E+09 3.073E+02 3.073E+02 2.764E+09
  PF2500   3.770E+09 1.000E+00
  PF2500   3.500E+18 1.661E+19 1.042E+19
  PF2500   7.284E+01 7.781E+01 1.321E+02 -2.524E+02 4.054E+01
  PF2500   7.359E+01 7.856E+01 1.355E+02 -2.601E+02 4.058E+01
  PF2500   7.151E+01 7.648E+01 1.337E+02 -2.577E+02 3.909E+01
end_S02
#======================================#
start_S03
  minE   -226.484421
  Vmin      0.00
  Vave   1961.03
  Vmax   3442.66
  V(0,0)    0.00
  const  1970.78
  SP  MIN01    0.00    0.00      0.00  +183.00  +3827.00  14370.85  4.478  0.785  -0.512  14072.21
  SP  MIN02  180.00  180.00    630.33   +79.00  +3899.00  14182.67  5.935  0.747  +0.311  14010.79
  SP  MIN03  197.44   81.21    700.78  +117.00  +3874.00  14232.38  6.023  0.738  +0.006  14002.14
  SP  MIN04    0.00  180.00   2185.97  +159.00  +3881.00  14220.17  3.929  0.744  -0.115  14011.77
  SP  TS_01  187.43  125.18    923.11  -264.00  +3908.00  14062.85  5.993  0.743  +0.199  13998.09
  SP  TS_02  281.41   64.52   1936.21  -151.20  +3855.00  14284.25  6.122  0.772  -0.125  14078.32
  SP  TS_03  180.00    0.00   2090.72  -402.70  +3876.00  13999.45  5.958  0.691  -0.305  13982.63
  SP  TS_04    6.11  222.40   2337.84  -286.20  +3901.00  14109.00  3.947  0.749  -0.149  14011.32
  SP  TS_05   60.27  191.24   2950.62  -161.40  +3876.00  14175.39  5.456  0.760  +0.069  14058.48
  SP  TS_06   68.76   66.22   2986.39  -171.00  +3857.00  14203.97  5.680  0.774  -0.472  14041.57
  SP  MAX01   61.99  219.37   3001.10  -245.10  +3897.00  14071.81  5.489  0.764  +0.097  14060.66
  SP  MAX02   78.98   13.29   3373.98  -376.40  +3876.00  14086.75  5.905  0.734  -0.481  14078.44
  SP  MAX03   60.02  124.43   3449.48  -317.40  +3892.00  14056.10  5.526  0.758  -0.221  14037.08
  TOR2DNS   273.49  451.80  617.10  631.11  787.43  794.38  808.46  870.86
  zpe1WHO  41.088
  zpeMSHO  41.088
  zpe2DNS   0.782
  zpeEHR   40.234
  zpeE2DT  41.016
  PF100    7.522E+03 7.544E+03 7.671E+03 1.092E+00 6.987E-02 4.907E+02
  PF100    1.202E-86 1.205E-86 1.760E-86 2.135E-02 6.987E-02 5.758E-86
  PF100    1.207E+06 1.000E+00
  PF100    1.450E-80 1.454E-80 2.123E-80
  PF100    6.500E-01 8.490E-01 5.405E+01 -4.556E+00 9.751E+00
  PF100    6.540E-01 8.530E-01 5.410E+01 -4.557E+00 1.005E+01
  PF100    6.590E-01 8.580E-01 5.418E+01 -4.560E+00 1.020E+01
  PF700    4.183E+06 1.372E+07 1.623E+07 1.869E+01 1.077E+01 9.351E+06
  PF700    6.215E-07 2.039E-06 2.540E-06 1.065E+01 1.077E+01 2.567E-06
  PF700    1.564E+08 1.000E+00
  PF700    9.721E+01 3.189E+02 3.972E+02
  PF700    1.057E+01 1.196E+01 8.488E+01 -4.745E+01 2.659E+01
  PF700    1.185E+01 1.324E+01 8.907E+01 -4.911E+01 2.766E+01
  PF700    1.190E+01 1.329E+01 8.948E+01 -4.934E+01 2.738E+01
  PF2500   7.643E+11 5.771E+12 5.472E+12 2.747E+02 2.347E+02 4.675E+12
  PF2500   1.956E+08 1.477E+09 1.421E+09 2.347E+02 2.347E+02 1.421E+09
  PF2500   3.770E+09 1.000E+00
  PF2500   7.376E+17 5.569E+18 5.358E+18
  PF2500   7.188E+01 7.685E+01 1.289E+02 -2.455E+02 4.056E+01
  PF2500   7.390E+01 7.887E+01 1.338E+02 -2.555E+02 4.073E+01
  PF2500   7.168E+01 7.665E+01 1.328E+02 -2.553E+02 3.898E+01
end_S03
#======================================#
start_S04
  minE   -153.222817
  Vmin      0.00
  Vave    903.06
  Vmax   2056.96
  V(0,0) 2056.96
  const   874.58
  SP  MIN01   62.63   55.04      0.00  +269.00  +3858.00  18727.77  0.743  2.613  -0.151  18422.95
  SP  MIN02  180.00   59.36     52.45  +254.00  +3866.00  18728.99  0.743  2.618  +0.034  18443.11
  SP  TS_01  121.59   55.83    510.06  -324.00  +3896.00  18577.45  0.741  2.619  -0.028  18430.65
  SP  TS_02    0.00   60.14    577.66  -346.00  +3889.00  18611.22  0.710  2.611  -0.213  18456.18
  SP  TS_03  180.00    0.00   1205.14  -257.70  +3864.00  18634.16  0.746  2.605  +0.032  18473.73
  SP  TS_04  292.50    4.73   1249.47  -267.00  +3858.00  18631.59  0.746  2.593  -0.144  18453.80
  SP  MAX01  238.28    3.48   1689.74  -329.50  +3895.00  18475.82  0.744  2.603  -0.031  18461.85
  SP  MAX02    0.00    0.00   2079.52  -394.60  +3902.00  18502.41  0.702  2.573  -0.216  18495.93
  TOR2DNS   281.66  281.66  281.66  282.87  282.87  282.87  312.51  312.51
  zpe1WHO  53.545
  zpeMSHO  53.545
  zpe2DNS   0.805
  zpeEHR   52.674
  zpeE2DT  53.479
  PF100    3.487E+03 8.555E+03 8.730E+03 2.770E+00 1.476E-01 4.651E+02
  PF100    3.324E-114 8.154E-114 1.161E-113 4.815E-02 1.476E-01 3.559E-113
  PF100    8.106E+05 1.000E+00
  PF100    2.694E-108 6.610E-108 9.411E-108
  PF100    6.230E-01 8.220E-01 5.146E+01 -4.325E+00 9.136E+00
  PF100    6.520E-01 8.510E-01 5.572E+01 -4.721E+00 9.355E+00
  PF100    6.640E-01 8.630E-01 5.370E+01 -4.507E+00 9.795E+00
  PF700    1.264E+06 3.710E+06 4.698E+06 1.842E+01 1.045E+01 2.666E+06
  PF700    2.424E-11 7.115E-11 9.448E-11 1.032E+01 1.045E+01 9.565E-11
  PF700    1.051E+08 1.000E+00
  PF700    2.547E-03 7.477E-03 9.929E-03
  PF700    1.036E+01 1.175E+01 8.141E+01 -4.524E+01 2.724E+01
  PF700    1.041E+01 1.180E+01 8.581E+01 -4.826E+01 2.723E+01
  PF700    1.037E+01 1.177E+01 8.404E+01 -4.706E+01 2.645E+01
  PF2500   3.893E+11 1.172E+12 1.017E+12 1.088E+02 9.252E+01 8.647E+11
  PF2500   8.118E+06 2.444E+07 2.149E+07 9.252E+01 9.252E+01 2.149E+07
  PF2500   2.533E+09 1.000E+00
  PF2500   2.057E+16 6.191E+16 5.444E+16
  PF2500   7.759E+01 8.256E+01 1.291E+02 -2.402E+02 4.534E+01
  PF2500   7.764E+01 8.261E+01 1.335E+02 -2.511E+02 4.534E+01
  PF2500   7.481E+01 7.978E+01 1.299E+02 -2.449E+02 4.351E+01
end_S04
#======================================#
start_S05
  minE   -264.168822
  Vmin      0.00
  Vave   4438.69
  Vmax   8086.16
  V(0,0)    0.00
  const  4559.96
  SP  MIN01    0.00    0.00      0.00  +159.00  +3881.00  15935.02  8.497  0.763  +0.433  15544.68
  SP  MIN02  180.00    0.00    237.69  +157.00  +3886.00  15922.99  8.284  0.745  -0.028  15548.41
  SP  MIN03    0.00  180.00   2958.96  +136.00  +3928.00  15791.57  8.564  0.695  -0.263  15518.81
  SP  MIN04  159.48  159.20   3867.58  +138.00  +3928.00  15733.39  8.306  0.714  -0.391  15462.91
  SP  TS_01   92.48    0.21   3730.41  -175.10  +3854.00  15695.97  9.511  0.763  +0.231  15379.97
  SP  TS_02  180.00  180.00   3995.54  -135.80  +3960.00  15660.95  7.715  0.684  -0.346  15499.51
  SP  TS_03  169.40  107.10   4324.31  -374.50  +3859.00  15506.57  8.442  0.759  -0.302  15395.42
  SP  TS_04    3.77  258.14   4330.67  -469.70  +3855.00  15471.89  8.567  0.741  +0.038  15380.88
  SP  TS_05   90.15  182.63   5734.21  -152.70  +3915.00  15608.29  9.632  0.721  -0.356  15366.38
  SP  MAX01  260.31   95.60   7979.00  -536.80  +3859.00  15223.88  9.562  0.742  +0.096  15206.48
  SP  MAX02   90.47   97.09   8094.00  -537.10  +3857.00  15223.31  9.627  0.745  -0.291  15205.40
  TOR2DNS   370.38  522.47  585.20  672.97  727.67  821.92  869.14  945.34
  zpe1WHO  45.561
  zpeMSHO  45.561
  zpe2DNS   1.059
  zpeEHR   44.444
  zpeE2DT  45.503
  PF100    1.553E+04 1.613E+04 1.686E+04 1.179E+00 5.479E-02 7.836E+02
  PF100    4.177E-96 4.337E-96 6.045E-96 5.718E-03 5.479E-02 5.792E-95
  PF100    1.586E+06 1.000E+00
  PF100    6.624E-90 6.878E-90 9.586E-90
  PF100    6.600E-01 8.590E-01 5.614E+01 -4.755E+00 9.894E+00
  PF100    6.840E-01 8.830E-01 5.645E+01 -4.762E+00 1.063E+01
  PF100    6.880E-01 8.870E-01 5.658E+01 -4.771E+00 1.074E+01
  PF700    1.893E+07 3.072E+07 3.692E+07 9.522E+00 4.696E+00 1.821E+07
  PF700    1.129E-07 1.832E-07 2.295E-07 4.447E+00 4.696E+00 2.423E-07
  PF700    2.056E+08 1.000E+00
  PF700    2.321E+01 3.768E+01 4.718E+01
  PF700    1.226E+01 1.365E+01 9.083E+01 -4.994E+01 3.116E+01
  PF700    1.254E+01 1.394E+01 9.220E+01 -5.061E+01 3.158E+01
  PF700    1.282E+01 1.421E+01 9.297E+01 -5.086E+01 3.245E+01
  PF2500   3.409E+13 9.170E+13 1.485E+14 1.314E+02 1.062E+02 1.200E+14
  PF2500   3.546E+09 9.541E+09 1.563E+10 1.062E+02 1.062E+02 1.563E+10
  PF2500   4.956E+09 1.000E+00
  PF2500   1.758E+19 4.728E+19 7.744E+19
  PF2500   8.316E+01 8.813E+01 1.415E+02 -2.657E+02 4.623E+01
  PF2500   8.640E+01 9.137E+01 1.448E+02 -2.706E+02 4.786E+01
  PF2500   8.664E+01 9.161E+01 1.459E+02 -2.730E+02 4.653E+01
end_S05
#======================================#
start_S06
  minE   -264.149214
  Vmin      0.00
  Vave   6120.75
  Vmax   8844.24
  V(0,0)    0.00
  const  6214.25
  SP  MIN01    0.00    0.00      0.00  +293.00  +3606.00  16135.95  6.527  0.803  -0.462  15549.67
  SP  MIN02  180.00  180.00   4155.09  +182.00  +3930.00  15706.66  6.891  0.722  +0.276  15495.96
  SP  MIN03  180.00    0.00   4328.04  +166.00  +3923.00  15774.64  6.908  0.683  -0.006  15494.39
  SP  MIN04    0.00  180.00   5551.39  +150.00  +3918.00  15758.86  5.525  0.717  +0.092  15529.84
  SP  TS_01  182.82   74.45   5371.86  -356.50  +3855.00  15493.75  6.848  0.758  +0.139  15367.49
  SP  TS_02   11.16  246.82   6244.71  -385.10  +3850.00  15501.81  5.882  0.756  -0.006  15422.67
  SP  TS_03  261.79    9.93   7284.58  -202.40  +3888.00  15534.20  8.509  0.722  -0.174  15338.04
  SP  TS_04   87.84  201.77   7732.53  -170.70  +3915.00  15410.02  8.666  0.735  +0.335  15315.99
  SP  MAX01  266.31   84.39   8045.06  -278.60  +3842.00  15262.92  8.771  0.778  +0.244  15231.21
  SP  MAX02   88.96   77.11   8849.88  -376.90  +3847.00  15246.99  8.631  0.770  -0.285  15200.61
  TOR2DNS   468.05  710.24  950.63 1149.00 1189.15 1388.18 1425.73 1625.44
  zpe1WHO  46.135
  zpeMSHO  46.135
  zpe2DNS   1.338
  zpeEHR   44.459
  zpeE2DT  45.797
  PF100    1.328E+04 1.328E+04 1.405E+04 1.032E+00 2.811E-02 3.828E+02
  PF100    1.984E-97 1.984E-97 1.150E-96 1.227E-03 2.811E-02 2.633E-95
  PF100    1.586E+06 1.000E+00
  PF100    3.146E-91 3.146E-91 1.823E-90
  PF100    6.220E-01 8.200E-01 5.544E+01 -4.724E+00 9.116E+00
  PF100    6.220E-01 8.200E-01 5.544E+01 -4.724E+00 9.116E+00
  PF100    6.270E-01 8.250E-01 5.560E+01 -4.735E+00 9.351E+00
  PF700    8.274E+06 8.305E+06 1.107E+07 3.459E+00 1.439E+00 4.606E+06
  PF700    3.266E-08 3.278E-08 5.573E-08 1.322E+00 1.439E+00 6.067E-08
  PF700    2.056E+08 1.000E+00
  PF700    6.715E+00 6.740E+00 1.146E+01
  PF700    1.164E+01 1.303E+01 8.830E+01 -4.878E+01 3.052E+01
  PF700    1.168E+01 1.307E+01 8.838E+01 -4.879E+01 3.109E+01
  PF700    1.207E+01 1.346E+01 8.950E+01 -4.919E+01 3.206E+01
  PF2500   1.014E+13 3.375E+13 6.830E+13 5.053E+01 3.860E+01 5.217E+13
  PF2500   9.395E+08 3.128E+09 6.775E+09 3.860E+01 3.860E+01 6.775E+09
  PF2500   4.956E+09 1.000E+00
  PF2500   4.656E+18 1.550E+19 3.358E+19
  PF2500   8.239E+01 8.736E+01 1.388E+02 -2.597E+02 4.627E+01
  PF2500   9.140E+01 9.636E+01 1.448E+02 -2.657E+02 4.920E+01
  PF2500   9.180E+01 9.677E+01 1.464E+02 -2.692E+02 4.760E+01
end_S06
#======================================#
start_S07
  minE   -289.184021
  Vmin      0.00
  Vave   2473.92
  Vmax   3952.26
  V(0,0)    0.00
  const  2531.49
  SP  MIN01    0.00    0.00      0.00  +138.00  +3904.00  18384.01  0.745  11.407  +0.510  18108.27
  SP  MIN02  199.72  151.42    932.55  +100.00  +3882.00  18395.73  0.783  13.450  -0.483  18087.53
  SP  MIN03    7.11  240.72   1730.78   +78.00  +3900.00  18365.61  0.738  13.512  +0.227  18125.03
  SP  TS_01  180.00  180.00   1100.67   -84.10  +3892.00  18340.83  0.778  12.707  -0.545  18077.60
  SP  TS_02    2.53  285.16   2003.36   -90.60  +3900.00  18343.80  0.746  13.430  +0.372  18131.37
  SP  TS_03  180.00    0.00   2058.45  -223.80  +3902.00  18154.51  0.721  11.441  -0.316  18076.89
  SP  TS_04   73.64  236.77   2568.51  -365.00  +3871.00  18101.00  0.774  13.667  +0.129  18035.04
  SP  TS_05    0.00  180.00   2610.65  -115.30  +3897.00  18273.98  0.731  11.748  +0.008  18053.76
  SP  TS_06  242.04   64.51   3226.94  -115.90  +3833.00  18128.54  0.762  13.252  -0.047  18024.60
  SP  TS_07   85.79  125.48   3530.69  -341.00  +3839.00  18062.46  0.769  13.326  -0.338  18008.02
  SP  MAX01  271.25   67.19   3270.61  -233.30  +3835.00  18031.56  0.772  13.358  +0.167  18021.17
  SP  MAX02  100.79   73.42   3937.81  -233.60  +3832.00  18028.03  0.774  13.428  -0.305  18018.13
  SP  MAX03   76.15  171.69   3961.30  -399.30  +3856.00  17988.27  0.768  11.900  -0.226  17972.28
  TOR2DNS   253.69  378.46  501.92  622.26  624.01  744.40  747.85  863.21
  zpe1WHO  52.563
  zpeMSHO  52.563
  zpe2DNS   0.725
  zpeEHR   51.774
  zpeE2DT  52.499
  PF100    1.938E+04 1.938E+04 2.044E+04 1.206E+00 9.839E-02 1.667E+03
  PF100    2.597E-111 2.597E-111 3.761E-111 3.135E-02 9.839E-02 1.180E-110
  PF100    1.720E+06 1.000E+00
  PF100    4.468E-105 4.468E-105 6.470E-105
  PF100    6.770E-01 8.760E-01 5.691E+01 -4.815E+00 1.032E+01
  PF100    6.770E-01 8.760E-01 5.691E+01 -4.815E+00 1.033E+01
  PF100    6.880E-01 8.860E-01 5.712E+01 -4.826E+00 1.054E+01
  PF700    5.129E+07 7.325E+07 1.057E+08 1.360E+01 8.223E+00 6.393E+07
  PF700    1.994E-09 2.847E-09 4.301E-09 8.075E+00 8.223E+00 4.380E-09
  PF700    2.230E+08 1.000E+00
  PF700    4.446E-01 6.349E-01 9.591E-01
  PF700    1.345E+01 1.484E+01 9.468E+01 -5.144E+01 3.422E+01
  PF700    1.443E+01 1.582E+01 9.678E+01 -5.193E+01 3.683E+01
  PF700    1.492E+01 1.631E+01 9.821E+01 -5.244E+01 3.742E+01
  PF2500   5.272E+14 1.958E+15 3.037E+15 3.052E+02 2.637E+02 2.624E+15
  PF2500   1.340E+10 4.975E+10 7.816E+10 2.637E+02 2.637E+02 7.816E+10
  PF2500   5.376E+09 1.000E+00
  PF2500   7.203E+19 2.675E+20 4.202E+20
  PF2500   9.264E+01 9.760E+01 1.509E+02 -2.797E+02 5.161E+01
  PF2500   9.548E+01 1.004E+02 1.547E+02 -2.862E+02 5.193E+01
  PF2500   9.387E+01 9.884E+01 1.549E+02 -2.884E+02 5.017E+01
end_S07
#======================================#
start_S08
  minE   -190.855295
  Vmin      0.00
  Vave   1413.28
  Vmax   2211.99
  V(0,0) 1103.34
  const  1429.10
  SP  MIN01    0.00   60.81      0.00  +226.00  +3431.00  20149.05  4.074  2.939  -1.498  19909.42
  SP  MIN02  149.74   54.23   1140.83   +67.00  +3433.00  19933.41  7.041  2.555  -0.609  19808.93
  SP  TS_01    0.00    0.00   1103.74  -233.00  +3432.00  20000.83  3.802  2.907  -1.509  19920.36
  SP  TS_02  180.00   61.04   1191.31   -61.30  +3435.00  19882.53  6.983  2.453  -0.592  19799.54
  SP  TS_03  180.00    0.00   1523.81  -153.40  +3436.00  19854.67  6.875  2.428  -0.683  19812.38
  SP  TS_04   75.11   63.02   1663.40  -129.10  +3422.00  19883.59  6.249  2.887  -1.071  19794.99
  SP  MAX01   90.10    2.64   2215.38  -187.80  +3424.00  19791.64  6.583  2.842  -0.976  19787.39
  TOR2DNS   221.40  221.40  221.40  426.26  426.26  426.26  446.06  446.06
  zpe1WHO  57.609
  zpeMSHO  57.609
  zpe2DNS   0.633
  zpeEHR   56.924
  zpeE2DT  57.557
  PF100    7.760E+03 7.760E+03 8.102E+03 1.100E+00 9.508E-02 7.003E+02
  PF100    9.730E-123 9.730E-123 1.320E-122 4.549E-02 9.508E-02 2.759E-122
  PF100    1.147E+06 1.000E+00
  PF100    1.116E-116 1.116E-116 1.515E-116
  PF100    6.520E-01 8.500E-01 5.403E+01 -4.553E+00 1.003E+01
  PF100    6.520E-01 8.500E-01 5.621E+01 -4.771E+00 1.003E+01
  PF100    6.660E-01 8.640E-01 5.425E+01 -4.561E+00 1.038E+01
  PF700    7.792E+06 1.598E+07 1.986E+07 1.544E+01 9.844E+00 1.266E+07
  PF700    8.048E-12 1.651E-11 2.129E-11 9.797E+00 9.844E+00 2.139E-11
  PF700    1.487E+08 1.000E+00
  PF700    1.197E-03 2.456E-03 3.167E-03
  PF700    1.195E+01 1.334E+01 8.798E+01 -4.825E+01 3.139E+01
  PF700    1.349E+01 1.488E+01 9.380E+01 -5.078E+01 3.381E+01
  PF700    1.320E+01 1.459E+01 9.163E+01 -4.955E+01 3.206E+01
  PF2500   2.192E+13 1.348E+14 8.562E+13 2.130E+02 1.875E+02 7.538E+13
  PF2500   2.017E+08 1.241E+09 7.963E+08 1.875E+02 1.875E+02 7.963E+08
  PF2500   3.585E+09 1.000E+00
  PF2500   7.233E+17 4.449E+18 2.855E+18
  PF2500   8.883E+01 9.380E+01 1.423E+02 -2.619E+02 5.113E+01
  PF2500   9.147E+01 9.644E+01 1.491E+02 -2.764E+02 5.127E+01
  PF2500   8.776E+01 9.273E+01 1.446E+02 -2.687E+02 4.932E+01
end_S08
#======================================#
start_S09
  minE   -190.854256
  Vmin      0.00
  Vave   1319.63
  Vmax   2495.70
  V(0,0) 1559.84
  const  1333.27
  SP  MIN01   11.16   67.99      0.00  +202.00  +3860.00  20068.70  5.276  0.791  -0.469  19766.77
  SP  MIN02  227.80   60.91    390.66  +111.00  +3859.00  20007.82  7.353  0.748  -0.153  19757.97
  SP  MIN03    0.00  180.00    498.21  +153.00  +3885.00  19961.01  5.069  0.746  -0.092  19760.92
  SP  MIN04  139.89  175.97   1043.16  +103.00  +3870.00  19924.33  7.247  0.765  +0.277  19732.55
  SP  MIN05  128.11   63.63   1092.76   +84.00  +3859.00  19930.92  7.284  0.769  -0.364  19744.96
  SP  TS_01    2.26  225.49    640.65  -261.10  +3902.00  19876.09  5.100  0.756  -0.163  19764.70
  SP  TS_02  161.87   72.10   1214.35  -102.90  +3862.00  19849.57  7.205  0.759  -0.203  19688.60
  SP  TS_03  180.00  180.00   1258.25   -99.80  +3881.00  19818.65  7.117  0.759  +0.318  19683.28
  SP  TS_04  118.79   14.40   1333.09  -308.80  +3885.00  19842.41  7.290  0.741  -0.483  19787.97
  SP  TS_05  216.62  133.54   1336.82  -321.80  +3896.00  19790.83  7.249  0.755  +0.247  19729.95
  SP  TS_06  137.99  118.58   1406.61  -285.10  +3894.00  19786.26  7.332  0.763  -0.001  19724.83
  SP  TS_07   77.47   49.69   1413.42  -121.40  +3862.00  19924.19  6.996  0.780  -0.556  19781.06
  SP  TS_08    0.00    0.00   1562.22  -472.90  +3937.00  19834.33  4.737  0.728  -0.524  19774.28
  SP  TS_09  307.73   50.83   1739.99  -176.00  +3883.00  19879.73  6.097  0.758  -0.324  19762.40
  SP  TS_10   69.09  178.96   2016.75  -154.60  +3848.00  19922.51  6.810  0.770  +0.089  19762.17
  SP  MAX01  174.65  126.43   1547.08  -293.30  +3900.00  19695.08  7.186  0.754  +0.159  19680.81
  SP  MAX02  313.70   28.45   1772.26  -251.50  +3902.00  19787.31  5.935  0.747  -0.443  19777.11
  SP  MAX03  180.00    0.00   2073.38  -401.10  +3904.00  19722.63  7.123  0.708  -0.399  19714.78
  SP  MAX04   65.81  127.82   2410.71  -350.40  +3886.00  19764.09  6.787  0.759  -0.186  19749.25
  SP  MAX05   69.00  235.75   2507.94  -353.50  +3891.00  19766.08  6.634  0.773  +0.088  19748.81
  TOR2DNS   278.16  278.16  460.91  460.92  610.85  612.83  625.55  625.55
  zpe1WHO  57.379
  zpeMSHO  57.379
  zpe2DNS   0.795
  zpeEHR   56.516
  zpeE2DT  57.311
  PF100    8.378E+03 1.696E+04 1.735E+04 2.197E+00 1.349E-01 1.066E+03
  PF100    3.338E-122 6.756E-122 9.732E-122 4.015E-02 1.349E-01 3.271E-121
  PF100    1.147E+06 1.000E+00
  PF100    3.830E-116 7.752E-116 1.117E-115
  PF100    6.440E-01 8.420E-01 5.410E+01 -4.568E+00 9.683E+00
  PF100    6.560E-01 8.550E-01 5.562E+01 -4.708E+00 1.030E+01
  PF100    6.670E-01 8.660E-01 5.578E+01 -4.713E+00 1.061E+01
  PF700    7.882E+06 4.722E+07 5.848E+07 3.718E+01 2.121E+01 3.336E+07
  PF700    9.603E-12 5.753E-11 7.482E-11 2.099E+01 2.121E+01 7.559E-11
  PF700    1.487E+08 1.000E+00
  PF700    1.428E-03 8.558E-03 1.113E-02
  PF700    1.212E+01 1.351E+01 8.825E+01 -4.827E+01 3.192E+01
  PF700    1.334E+01 1.473E+01 9.356E+01 -5.076E+01 3.324E+01
  PF700    1.320E+01 1.460E+01 9.378E+01 -5.105E+01 3.194E+01
  PF2500   2.591E+13 3.523E+14 2.549E+14 4.050E+02 3.451E+02 2.172E+14
  PF2500   2.498E+08 3.396E+09 2.491E+09 3.451E+02 3.451E+02 2.491E+09
  PF2500   3.585E+09 1.000E+00
  PF2500   8.955E+17 1.218E+19 8.930E+18
  PF2500   8.920E+01 9.417E+01 1.428E+02 -2.627E+02 5.111E+01
  PF2500   9.116E+01 9.613E+01 1.487E+02 -2.757E+02 5.123E+01
  PF2500   8.778E+01 9.275E+01 1.467E+02 -2.741E+02 4.929E+01
end_S09
#======================================#
start_S10
  minE   -190.869302
  Vmin      0.00
  Vave   1481.14
  Vmax   2492.77
  V(0,0)    0.00
  const  1523.90
  SP  MIN01    0.00    0.00      0.00  +182.00  +3893.00  19934.88  0.729  2.932  +0.068  19616.82
  SP  MIN02  180.00    0.00   1128.10  +133.00  +3917.00  19793.03  0.728  2.924  -0.103  19607.17
  SP  TS_01    0.00   60.17    657.77  -185.60  +3893.00  19847.07  0.730  2.927  +0.068  19612.95
  SP  TS_02   88.62    1.23   1735.82  -333.60  +3844.00  19625.06  0.765  2.930  -0.015  19513.51
  SP  TS_03  166.58   60.60   2104.54  -227.50  +3921.00  19641.64  0.725  2.909  -0.100  19591.44
  SP  MAX01  180.00   60.02   2105.86  -229.50  +3930.00  19604.81  0.722  2.908  -0.103  19597.35
  SP  MAX02   89.68   62.15   2495.43  -317.40  +3838.00  19519.93  0.768  2.924  -0.018  19509.31
  TOR2DNS   294.44  294.44  294.44  459.64  459.64  459.69  610.14  610.97
  zpe1WHO  56.997
  zpeMSHO  56.997
  zpe2DNS   0.842
  zpeEHR   56.087
  zpeE2DT  56.929
  PF100    8.557E+03 8.557E+03 8.933E+03 1.109E+00 6.547E-02 5.275E+02
  PF100    2.338E-121 2.338E-121 3.428E-121 1.603E-02 6.547E-02 1.400E-120
  PF100    1.147E+06 1.000E+00
  PF100    2.682E-115 2.682E-115 3.933E-115
  PF100    6.430E-01 8.420E-01 5.414E+01 -4.572E+00 9.504E+00
  PF100    6.430E-01 8.420E-01 5.632E+01 -4.790E+00 9.505E+00
  PF100    6.540E-01 8.530E-01 5.433E+01 -4.581E+00 9.823E+00
  PF700    9.087E+06 1.171E+07 1.538E+07 8.376E+00 4.665E+00 8.564E+06
  PF700    1.458E-11 1.879E-11 2.589E-11 4.573E+00 4.665E+00 2.641E-11
  PF700    1.487E+08 1.000E+00
  PF700    2.168E-03 2.794E-03 3.851E-03
  PF700    1.238E+01 1.377E+01 8.890E+01 -4.846E+01 3.224E+01
  PF700    1.308E+01 1.447E+01 9.260E+01 -5.034E+01 3.403E+01
  PF700    1.305E+01 1.444E+01 9.091E+01 -4.920E+01 3.303E+01
  PF2500   3.546E+13 8.851E+13 7.525E+13 8.361E+01 7.058E+01 6.352E+13
  PF2500   3.692E+08 9.214E+08 7.941E+08 7.058E+01 7.058E+01 7.941E+08
  PF2500   3.585E+09 1.000E+00
  PF2500   1.324E+18 3.304E+18 2.847E+18
  PF2500   8.959E+01 9.455E+01 1.435E+02 -2.643E+02 5.112E+01
  PF2500   9.151E+01 9.648E+01 1.483E+02 -2.743E+02 5.132E+01
  PF2500   8.842E+01 9.339E+01 1.446E+02 -2.680E+02 4.942E+01
end_S10
#======================================#
start_S11
  minE   -190.862329
  Vmin      0.00
  Vave    610.39
  Vmax   1156.79
  V(0,0)    0.00
  const   629.76
  SP  MIN01    0.00    0.00      0.00   +87.00  +3920.00  19970.45  0.714  2.887  -0.106  19757.03
  SP  MIN02  208.97    3.71    392.64  +146.00  +3897.00  19859.59  0.742  2.928  +0.073  19708.11
  SP  TS_01    0.00   60.73    169.21  -100.00  +3905.00  19902.83  0.723  2.901  -0.099  19729.25
  SP  TS_02  180.00    0.00    413.27  -118.40  +3928.00  19800.50  0.729  2.926  +0.086  19724.01
  SP  TS_03   77.11    8.29    715.05  -255.50  +3847.00  19745.91  0.768  2.930  -0.026  19646.61
  SP  TS_04  180.00   59.49    759.38  -140.30  +3927.00  19750.84  0.726  2.919  +0.087  19714.84
  SP  MAX01   81.37   67.24   1159.05  -291.30  +3835.00  19655.95  0.769  2.918  -0.018  19640.14
  TOR2DNS   198.91  199.30  199.30  272.98  272.98  279.19  323.29  344.53
  zpe1WHO  57.098
  zpeMSHO  57.098
  zpe2DNS   0.569
  zpeEHR   56.488
  zpeE2DT  57.057
  PF100    1.168E+04 1.204E+04 1.507E+04 1.578E+00 1.887E-01 1.802E+03
  PF100    1.912E-121 1.971E-121 3.042E-121 9.023E-02 1.887E-01 6.362E-121
  PF100    1.147E+06 1.000E+00
  PF100    2.194E-115 2.262E-115 3.490E-115
  PF100    7.170E-01 9.160E-01 5.549E+01 -4.634E+00 1.070E+01
  PF100    7.410E-01 9.400E-01 5.798E+01 -4.858E+00 1.168E+01
  PF100    7.520E-01 9.510E-01 5.636E+01 -4.685E+00 1.137E+01
  PF700    2.264E+07 5.278E+07 4.732E+07 2.145E+01 1.433E+01 3.161E+07
  PF700    3.376E-11 7.869E-11 7.269E-11 1.425E+01 1.433E+01 7.307E-11
  PF700    1.487E+08 1.000E+00
  PF700    5.021E-03 1.171E-02 1.081E-02
  PF700    1.246E+01 1.385E+01 9.083E+01 -4.973E+01 3.197E+01
  PF700    1.302E+01 1.441E+01 9.550E+01 -5.244E+01 3.227E+01
  PF700    1.238E+01 1.378E+01 9.219E+01 -5.076E+01 3.065E+01
  PF2500   8.877E+13 2.886E+14 1.130E+14 1.256E+02 1.120E+02 1.007E+14
  PF2500   9.054E+08 2.943E+09 1.162E+09 1.120E+02 1.120E+02 1.162E+09
  PF2500   3.585E+09 1.000E+00
  PF2500   3.246E+18 1.055E+19 4.166E+18
  PF2500   8.952E+01 9.449E+01 1.453E+02 -2.689E+02 5.110E+01
  PF2500   9.026E+01 9.523E+01 1.502E+02 -2.802E+02 5.113E+01
  PF2500   8.626E+01 9.123E+01 1.445E+02 -2.701E+02 4.918E+01
end_S11
#======================================#
start_S12
  minE   -228.516006
  Vmin      0.00
  Vave   2089.17
  Vmax   4296.80
  V(0,0)  670.77
  const  2129.79
  SP  MIN01  180.00    0.00      0.00  +139.00  +3396.00  21229.19  6.778  3.028  -0.785  21056.28
  SP  MIN02    0.00    0.00    670.93  +157.00  +3411.00  21195.03  6.819  2.955  +0.395  21034.95
  SP  TS_01  180.00   60.20    373.33  -129.40  +3391.00  21171.91  6.643  3.030  -0.780  21070.23
  SP  TS_02    0.00   60.18   1127.22  -148.30  +3409.00  21140.92  6.882  2.943  +0.378  21059.37
  SP  TS_03   91.77    2.27   3571.73  -173.20  +3389.00  20988.39  8.061  3.023  -0.227  20888.28
  SP  MAX01   93.13   61.67   4312.02  -198.00  +3385.00  20901.97  8.021  3.007  -0.244  20896.74
  TOR2DNS   151.53  151.54  151.54  274.30  274.30  274.84  322.91  322.91
  zpe1WHO  60.697
  zpeMSHO  60.697
  zpe2DNS   0.433
  zpeEHR   60.203
  zpeE2DT  60.636
  PF100    1.951E+04 1.951E+04 2.112E+04 1.338E+00 2.162E-01 3.412E+03
  PF100    4.357E-129 4.357E-129 6.415E-129 1.512E-01 2.162E-01 9.172E-129
  PF100    1.521E+06 1.000E+00
  PF100    6.627E-123 6.627E-123 9.758E-123
  PF100    7.100E-01 9.090E-01 5.701E+01 -4.792E+00 1.126E+01
  PF100    7.110E-01 9.090E-01 5.920E+01 -5.010E+00 1.128E+01
  PF100    7.390E-01 9.380E-01 5.745E+01 -4.808E+00 1.181E+01
  PF700    9.296E+07 1.200E+08 1.531E+08 1.937E+01 1.425E+01 1.126E+08
  PF700    1.043E-11 1.346E-11 1.794E-11 1.419E+01 1.425E+01 1.802E-11
  PF700    1.972E+08 1.000E+00
  PF700    2.056E-03 2.653E-03 3.538E-03
  PF700    1.403E+01 1.542E+01 9.645E+01 -5.209E+01 3.599E+01
  PF700    1.445E+01 1.584E+01 9.973E+01 -5.397E+01 3.661E+01
  PF700    1.438E+01 1.577E+01 9.794E+01 -5.278E+01 3.641E+01
  PF2500   3.311E+15 5.836E+15 6.490E+15 2.005E+02 1.837E+02 5.948E+15
  PF2500   1.636E+10 2.885E+10 3.248E+10 1.837E+02 1.837E+02 3.248E+10
  PF2500   4.753E+09 1.000E+00
  PF2500   7.778E+19 1.371E+20 1.544E+20
  PF2500   1.009E+02 1.058E+02 1.576E+02 -2.882E+02 5.699E+01
  PF2500   1.017E+02 1.067E+02 1.613E+02 -2.965E+02 5.706E+01
  PF2500   1.006E+02 1.055E+02 1.588E+02 -2.916E+02 5.584E+01
end_S12
#======================================#
start_S13
  minE   -228.498044
  Vmin      0.00
  Vave   1513.19
  Vmax   3051.53
  V(0,0) 1304.50
  const  1545.20
  SP  MIN01  180.00    0.00      0.00  +101.00  +3923.00  21325.06  8.837  0.705  -0.066  21113.21
  SP  MIN02   43.50    8.86    109.30  +136.00  +3870.00  21413.05  9.112  0.763  -0.464  21087.25
  SP  MIN03  179.58  165.02    156.27   +87.00  +3921.00  21221.81  8.722  0.738  +0.328  21099.45
  SP  MIN04    0.00  180.00   1049.97   +81.00  +3924.00  21283.87  7.071  0.728  +0.128  21154.62
  SP  TS_01  180.00  180.00    158.02   -63.60  +3931.00  21183.42  8.721  0.733  +0.336  21103.45
  SP  TS_02  185.40   74.20    621.77  -268.10  +3845.00  21116.58  8.721  0.774  +0.135  21017.98
  SP  TS_03  249.24    7.93   1028.46  -108.50  +3882.00  21256.10  10.503  0.748  -0.186  21041.41
  SP  TS_04    0.00    0.00   1304.56  -262.80  +3980.00  21326.76  6.409  0.721  -0.443  21176.99
  SP  TS_05   14.33   95.01   1372.59  -256.50  +3822.00  21153.49  7.560  0.787  -0.217  21068.26
  SP  TS_06   90.23  188.07   2594.19  -107.40  +3921.00  21062.33  10.785  0.740  +0.348  21001.02
  SP  MAX01   84.43  258.98   2816.08  -268.10  +3843.00  20942.75  10.831  0.781  +0.295  20933.90
  SP  MAX02   97.89   91.25   3054.65  -316.30  +3838.00  20932.75  10.794  0.780  -0.217  20913.97
  TOR2DNS   200.72  269.46  304.94  381.60  389.20  389.20  406.71  407.42
  zpe1WHO  60.971
  zpeMSHO  60.971
  zpe2DNS   0.574
  zpeEHR   60.366
  zpeE2DT  60.940
  PF100    2.152E+04 4.757E+04 3.827E+04 2.017E+00 2.268E-01 4.302E+03
  PF100    1.210E-129 2.675E-129 2.525E-129 1.123E-01 2.268E-01 5.098E-129
  PF100    1.521E+06 1.000E+00
  PF100    1.840E-123 4.068E-123 3.841E-123
  PF100    7.260E-01 9.250E-01 5.736E+01 -4.811E+00 1.128E+01
  PF100    8.550E-01 1.054E+00 6.023E+01 -4.969E+00 1.265E+01
  PF100    8.350E-01 1.034E+00 5.959E+01 -4.926E+00 1.267E+01
  PF700    1.149E+08 6.473E+08 5.044E+08 4.400E+01 2.943E+01 3.374E+08
  PF700    1.058E-11 5.962E-11 4.754E-11 2.913E+01 2.943E+01 4.803E-11
  PF700    1.972E+08 1.000E+00
  PF700    2.087E-03 1.176E-02 9.373E-03
  PF700    1.428E+01 1.567E+01 9.722E+01 -5.239E+01 3.664E+01
  PF700    1.468E+01 1.607E+01 1.012E+02 -5.479E+01 3.698E+01
  PF700    1.470E+01 1.609E+01 1.008E+02 -5.444E+01 3.679E+01
  PF2500   4.926E+15 3.653E+16 2.317E+16 4.323E+02 3.851E+02 2.064E+16
  PF2500   2.304E+10 1.709E+11 1.091E+11 3.851E+02 3.851E+02 1.091E+11
  PF2500   4.753E+09 1.000E+00
  PF2500   1.095E+20 8.122E+20 5.184E+20
  PF2500   1.012E+02 1.061E+02 1.585E+02 -2.902E+02 5.685E+01
  PF2500   1.019E+02 1.068E+02 1.628E+02 -3.002E+02 5.692E+01
  PF2500   9.986E+01 1.048E+02 1.611E+02 -2.979E+02 5.528E+01
end_S13
#======================================#
start_S14
  minE   -228.504086
  Vmin      0.00
  Vave   2484.54
  Vmax   4376.82
  V(0,0) 1374.31
  const  2512.13
  SP  MIN01  180.00    0.00      0.00  +156.00  +3910.00  21312.26  8.562  0.727  -0.033  21037.40
  SP  MIN02  167.48  129.99   1041.63  +163.00  +3853.00  21208.36  8.617  0.778  -0.346  20971.04
  SP  MIN03   28.56    1.57   1293.14  +100.00  +3903.00  21203.82  8.307  0.748  +0.373  20939.17
  SP  MIN04   40.82  212.74   2069.43  +130.00  +3880.00  21123.24  8.780  0.746  -0.137  20914.19
  SP  TS_01  175.57   86.42   1203.16  -261.20  +3847.00  21059.00  8.618  0.784  -0.214  20951.09
  SP  TS_02    0.00    0.00   1374.35   -80.60  +3908.00  21160.99  7.740  0.748  +0.393  20949.62
  SP  TS_03  180.00  180.00   1527.10  -298.00  +3957.00  21125.01  8.022  0.720  -0.336  21034.89
  SP  TS_04   29.85  281.08   2645.33  -339.30  +3855.00  20922.58  8.372  0.771  +0.222  20850.89
  SP  TS_05    0.00  180.00   2668.15  -185.20  +3930.00  20995.20  7.749  0.723  -0.218  20948.02
  SP  TS_06  100.49    0.71   2835.17  -151.90  +3881.00  21121.50  9.957  0.750  +0.184  20897.63
  SP  TS_07  100.13  185.06   3113.91  -142.40  +3913.00  21036.23  9.997  0.748  -0.344  20902.97
  SP  MAX01   99.32   81.47   4348.89  -370.80  +3845.00  20827.37  9.959  0.773  -0.225  20811.99
  SP  MAX02  254.26   80.68   4403.10  -383.60  +3849.00  20824.92  9.804  0.778  +0.121  20809.29
  TOR2DNS   248.31  390.86  532.10  586.63  672.19  732.45  811.25  876.81
  zpe1WHO  60.935
  zpeMSHO  60.935
  zpe2DNS   0.710
  zpeEHR   60.149
  zpeE2DT  60.859
  PF100    1.674E+04 1.674E+04 1.769E+04 1.157E+00 9.325E-02 1.426E+03
  PF100    1.131E-129 1.131E-129 1.751E-129 3.249E-02 9.325E-02 5.027E-129
  PF100    1.521E+06 1.000E+00
  PF100    1.721E-123 1.721E-123 2.664E-123
  PF100    6.650E-01 8.640E-01 5.625E+01 -4.761E+00 1.013E+01
  PF100    6.650E-01 8.640E-01 5.625E+01 -4.761E+00 1.013E+01
  PF100    6.770E-01 8.750E-01 5.648E+01 -4.772E+00 1.041E+01
  PF700    5.090E+07 8.879E+07 1.194E+08 1.272E+01 7.741E+00 7.264E+07
  PF700    4.813E-12 8.396E-12 1.192E-11 7.636E+00 7.741E+00 1.208E-11
  PF700    1.972E+08 1.000E+00
  PF700    9.490E-04 1.656E-03 2.351E-03
  PF700    1.419E+01 1.558E+01 9.547E+01 -5.125E+01 3.688E+01
  PF700    1.564E+01 1.703E+01 9.865E+01 -5.203E+01 4.015E+01
  PF700    1.572E+01 1.711E+01 9.935E+01 -5.244E+01 3.972E+01
  PF2500   2.151E+15 1.293E+16 1.404E+16 2.613E+02 2.265E+02 1.217E+16
  PF2500   1.014E+10 6.094E+10 6.716E+10 2.265E+02 2.265E+02 6.716E+10
  PF2500   4.753E+09 1.000E+00
  PF2500   4.818E+19 2.896E+20 3.192E+20
  PF2500   1.012E+02 1.061E+02 1.569E+02 -2.861E+02 5.685E+01
  PF2500   1.046E+02 1.095E+02 1.618E+02 -2.950E+02 5.716E+01
  PF2500   1.025E+02 1.074E+02 1.611E+02 -2.954E+02 5.546E+01
end_S14
#======================================#
start_S15
  minE   -228.513668
  Vmin      0.00
  Vave   1847.79
  Vmax   3644.27
  V(0,0)    0.00
  const  1899.12
  SP  MIN01    0.00    0.00      0.00  +159.00  +3365.00  21386.01  5.897  3.034  -0.694  21217.31
  SP  MIN02  180.00    0.00    740.07  +120.00  +3365.00  21272.57  7.290  2.964  +0.596  21132.25
  SP  TS_01    0.00   58.76    340.84   -19.40  +3362.00  21270.30  6.342  3.014  -0.700  21176.66
  SP  TS_02  180.00   60.40   1009.58  -102.00  +3372.00  21197.60  7.123  2.978  +0.613  21117.81
  SP  TS_03  265.73    1.21   3232.20  -175.60  +3341.00  21042.22  8.740  3.025  +0.042  20963.30
  SP  MAX01   92.88   58.96   3654.25  -179.30  +3351.00  20957.36  8.807  3.024  +0.043  20948.86
  TOR2DNS   156.10  156.11  156.11  290.35  290.35  290.85  319.35  319.35
  zpe1WHO  61.146
  zpeMSHO  61.146
  zpe2DNS   0.446
  zpeEHR   60.663
  zpeE2DT  61.110
  PF100    1.948E+04 1.948E+04 2.269E+04 1.320E+00 1.995E-01 3.429E+03
  PF100    4.557E-130 4.558E-130 6.363E-130 1.397E-01 1.995E-01 9.084E-130
  PF100    1.521E+06 1.000E+00
  PF100    6.931E-124 6.932E-124 9.678E-124
  PF100    7.120E-01 9.110E-01 5.702E+01 -4.792E+00 1.137E+01
  PF100    7.120E-01 9.110E-01 5.921E+01 -5.010E+00 1.139E+01
  PF100    7.510E-01 9.500E-01 5.772E+01 -4.822E+00 1.213E+01
  PF700    8.566E+07 1.246E+08 1.793E+08 2.067E+01 1.506E+01 1.307E+08
  PF700    6.961E-12 1.013E-11 1.495E-11 1.500E+01 1.506E+01 1.502E-11
  PF700    1.972E+08 1.000E+00
  PF700    1.373E-03 1.997E-03 2.948E-03
  PF700    1.384E+01 1.523E+01 9.600E+01 -5.198E+01 3.568E+01
  PF700    1.446E+01 1.585E+01 9.983E+01 -5.403E+01 3.661E+01
  PF700    1.438E+01 1.577E+01 9.825E+01 -5.301E+01 3.639E+01
  PF2500   2.633E+15 6.082E+15 7.321E+15 2.189E+02 2.001E+02 6.692E+15
  PF2500   1.189E+10 2.747E+10 3.330E+10 2.001E+02 2.001E+02 3.330E+10
  PF2500   4.753E+09 1.000E+00
  PF2500   5.652E+19 1.305E+20 1.583E+20
  PF2500   1.005E+02 1.055E+02 1.570E+02 -2.871E+02 5.696E+01
  PF2500   1.017E+02 1.067E+02 1.613E+02 -2.967E+02 5.705E+01
  PF2500   1.002E+02 1.052E+02 1.589E+02 -2.922E+02 5.562E+01
end_S15
#======================================#
start_S16
  minE   -155.243239
  Vmin      0.00
  Vave   1030.07
  Vmax   1939.43
  V(0,0) 1743.23
  const  1013.69
  SP  MIN01  116.65   58.76      0.00  +105.00  +3384.00  25610.69  7.763  2.841  -0.406  25430.72
  SP  MIN02    0.00   60.28    289.93  +151.00  +3395.00  25652.90  5.330  2.953  -1.278  25432.14
  SP  TS_01  180.00   59.98    751.70  -122.70  +3382.00  25513.36  7.570  2.675  -0.126  25380.35
  SP  TS_02   50.90   58.63    812.93  -133.40  +3386.00  25573.31  6.552  2.945  -1.034  25438.73
  SP  TS_03  119.35  358.63   1149.83  -248.50  +3385.00  25540.23  7.784  2.806  -0.398  25481.43
  SP  TS_04    0.00    0.00   1748.77  -274.00  +3400.00  25566.43  5.072  2.912  -1.259  25529.45
  SP  MAX01  180.00    0.00   1927.43  -255.10  +3382.00  25449.53  7.623  2.640  -0.172  25448.01
  SP  MAX02  303.80    0.21   1951.13  -253.70  +3383.00  25485.47  6.747  2.923  -1.002  25480.43
  TOR2DNS   170.94  170.94  170.94  170.94  170.94  170.94  274.16  274.16
  zpe1WHO  73.225
  zpeMSHO  73.225
  zpe2DNS   0.489
  zpeEHR   72.710
  zpeE2DT  73.199
  PF100    1.144E+04 2.297E+04 2.347E+04 2.704E+00 3.901E-01 3.386E+03
  PF100    1.070E-156 2.149E-156 2.500E-156 2.312E-01 3.901E-01 4.220E-156
  PF100    1.089E+06 1.000E+00
  PF100    1.166E-150 2.341E-150 2.723E-150
  PF100    7.090E-01 9.080E-01 5.527E+01 -4.619E+00 1.078E+01
  PF100    7.130E-01 9.110E-01 5.888E+01 -4.976E+00 1.094E+01
  PF100    7.220E-01 9.210E-01 5.683E+01 -4.762E+00 1.119E+01
  PF700    3.033E+07 7.115E+07 9.366E+07 3.873E+01 2.750E+01 6.649E+07
  PF700    4.175E-16 9.793E-16 1.313E-15 2.726E+01 2.750E+01 1.325E-15
  PF700    1.412E+08 1.000E+00
  PF700    5.895E-08 1.383E-07 1.854E-07
  PF700    1.357E+01 1.496E+01 9.290E+01 -5.007E+01 3.655E+01
  PF700    1.369E+01 1.509E+01 9.695E+01 -5.278E+01 3.664E+01
  PF700    1.378E+01 1.517E+01 9.544E+01 -5.164E+01 3.598E+01
  PF2500   1.606E+15 4.072E+15 3.968E+15 3.046E+02 2.760E+02 3.596E+15
  PF2500   6.377E+08 1.617E+09 1.584E+09 2.760E+02 2.760E+02 1.584E+09
  PF2500   3.404E+09 1.000E+00
  PF2500   2.170E+18 5.503E+18 5.390E+18
  PF2500   1.064E+02 1.114E+02 1.577E+02 -2.830E+02 6.181E+01
  PF2500   1.066E+02 1.115E+02 1.618E+02 -2.931E+02 6.182E+01
  PF2500   1.039E+02 1.089E+02 1.585E+02 -2.875E+02 5.998E+01
end_S16
#======================================#
start_S17
  minE   -155.243948
  Vmin      0.00
  Vave    224.48
  Vmax    414.82
  V(0,0)    0.00
  const   230.35
  SP  MIN01    0.00    0.00      0.00  +112.00  +3327.00  25516.61  2.917  2.917  -0.051  25395.50
  SP  TS_01   60.30    0.00    228.69   -91.00  +3338.00  25422.02  2.925  2.931  -0.034  25366.31
  SP  MAX01   59.88   59.88    414.81   -82.50  +3348.00  25341.37  2.935  2.935  -0.021  25340.85
  TOR2DNS   109.51  109.65  109.65  109.65  109.65  109.80  109.80  109.80
  zpe1WHO  72.956
  zpeMSHO  72.956
  zpe2DNS   0.313
  zpeEHR   72.609
  zpeE2DT  72.922
  PF100    6.840E+03 6.840E+03 1.351E+04 1.930E+00 4.521E-01 3.165E+03
  PF100    2.478E-156 2.478E-156 5.781E-156 3.992E-01 4.521E-01 6.546E-156
  PF100    1.089E+06 1.000E+00
  PF100    2.698E-150 2.698E-150 6.297E-150
  PF100    7.560E-01 9.550E-01 5.472E+01 -4.517E+00 1.165E+01
  PF100    7.560E-01 9.550E-01 5.909E+01 -4.954E+00 1.165E+01
  PF100    8.750E-01 1.073E+00 5.726E+01 -4.652E+00 1.206E+01
  PF700    3.172E+07 3.172E+07 5.968E+07 2.349E+01 1.878E+01 4.773E+07
  PF700    5.297E-16 5.297E-16 1.021E-15 1.876E+01 1.878E+01 1.022E-15
  PF700    1.412E+08 1.000E+00
  PF700    7.480E-08 7.480E-08 1.441E-07
  PF700    1.384E+01 1.523E+01 9.337E+01 -5.013E+01 3.653E+01
  PF700    1.384E+01 1.523E+01 9.773E+01 -5.319E+01 3.653E+01
  PF700    1.303E+01 1.443E+01 9.348E+01 -5.101E+01 3.465E+01
  PF2500   1.914E+15 1.914E+15 1.395E+15 9.838E+01 9.237E+01 1.310E+15
  PF2500   8.021E+08 8.021E+08 5.885E+08 9.237E+01 9.237E+01 5.885E+08
  PF2500   3.404E+09 1.000E+00
  PF2500   2.730E+18 2.730E+18 2.003E+18
  PF2500   1.066E+02 1.116E+02 1.582E+02 -2.839E+02 6.181E+01
  PF2500   1.066E+02 1.116E+02 1.626E+02 -2.948E+02 6.181E+01
  PF2500   1.023E+02 1.073E+02 1.558E+02 -2.823E+02 5.984E+01
end_S17
#======================================#
start_S18
  minE   -266.151290
  Vmin      0.00
  Vave   2975.48
  Vmax   5645.15
  V(0,0) 1173.86
  const  3023.27
  SP  MIN01  180.00    0.00      0.00  +141.00  +3403.00  22697.94  7.619  4.822  -2.615  22497.12
  SP  MIN02  180.00  180.00    221.23  +106.00  +3403.00  22593.15  8.728  8.169  -4.505  22418.22
  SP  MIN03   32.42    0.08   1041.63   +89.00  +3404.00  22650.65  8.508  7.089  -5.056  22474.07
  SP  MIN04   32.79  179.45   1213.48  +129.00  +3405.00  22563.59  5.652  6.805  -2.213  22398.60
  SP  TS_01    0.00    0.00   1173.75   -86.80  +3405.00  22633.35  8.471  7.047  -5.221  22516.13
  SP  TS_02    0.00  180.00   1355.69  -110.90  +3407.00  22535.46  5.056  6.646  -1.806  22441.12
  SP  TS_03  100.56  359.63   2244.57  -160.90  +3398.00  22491.02  8.222  6.257  -3.622  22397.42
  SP  TS_04  100.57  180.05   2383.93  -161.90  +3398.00  22410.83  8.129  7.599  -3.724  22319.74
  SP  TS_05  180.14   92.58   3756.53  -227.30  +3401.00  22347.81  8.403  7.208  -3.780  22259.18
  SP  TS_06   38.98   92.09   4624.55  -222.00  +3401.00  22327.79  6.662  6.545  -2.835  22245.28
  SP  TS_07   39.22  266.74   4662.96  -216.90  +3400.00  22317.98  8.096  7.547  -4.160  22243.25
  SP  MAX01    0.96  267.32   4929.18  -248.20  +3401.00  22291.00  7.072  7.031  -3.388  22272.72
  SP  MAX02  257.67   92.84   5636.77  -239.80  +3395.00  22199.12  9.134  7.880  -4.496  22187.97
  SP  MAX03  102.62   92.72   5659.81  -240.50  +3396.00  22200.47  7.380  6.233  -2.668  22188.86
  TOR2DNS   188.72  320.09  385.35  433.24  450.33  498.17  562.15  579.45
  zpe1WHO  64.897
  zpeMSHO  64.897
  zpe2DNS   0.540
  zpeEHR   64.323
  zpeE2DT  64.862
  PF100    3.901E+04 4.554E+04 4.695E+04 1.292E+00 1.552E-01 5.641E+03
  PF100    5.790E-138 6.758E-138 8.292E-138 8.552E-02 1.552E-01 1.505E-137
  PF100    1.928E+06 1.000E+00
  PF100    1.116E-131 1.303E-131 1.599E-131
  PF100    7.500E-01 9.490E-01 5.926E+01 -4.977E+00 1.238E+01
  PF100    8.010E-01 9.990E-01 6.007E+01 -5.007E+00 1.317E+01
  PF100    8.030E-01 1.002E+00 6.015E+01 -5.014E+00 1.334E+01
  PF700    6.128E+08 1.537E+09 1.838E+09 2.747E+01 1.879E+01 1.257E+09
  PF700    3.358E-12 8.421E-12 1.032E-11 1.864E+01 1.879E+01 1.041E-11
  PF700    2.500E+08 1.000E+00
  PF700    8.394E-04 2.105E-03 2.581E-03
  PF700    1.574E+01 1.713E+01 1.031E+02 -5.504E+01 4.034E+01
  PF700    1.670E+01 1.809E+01 1.063E+02 -5.632E+01 4.208E+01
  PF700    1.695E+01 1.835E+01 1.070E+02 -5.657E+01 4.259E+01
  PF2500   2.154E+17 1.166E+18 1.625E+18 5.411E+02 4.854E+02 1.458E+18
  PF2500   4.571E+11 2.475E+12 3.473E+12 4.854E+02 4.854E+02 3.473E+12
  PF2500   6.025E+09 1.000E+00
  PF2500   2.754E+21 1.491E+22 2.092E+22
  PF2500   1.122E+02 1.172E+02 1.709E+02 -3.102E+02 6.272E+01
  PF2500   1.143E+02 1.192E+02 1.751E+02 -3.185E+02 6.290E+01
  PF2500   1.141E+02 1.190E+02 1.757E+02 -3.202E+02 6.179E+01
end_S18
#======================================#
start_S19
  minE   -342.672656
  Vmin      0.00
  Vave   1249.88
  Vmax   2314.69
  V(0,0) 2314.69
  const  1248.95
  SP  MIN01   25.09   65.84      0.00   +91.00  +3859.00  31564.72  16.607  0.811  -0.675  31315.92
  SP  MIN02    0.00  180.00    713.29   +50.00  +3891.00  31410.80  15.866  0.776  +0.442  31284.24
  SP  TS_01    2.66  222.53    821.27  -243.70  +3903.00  31333.12  16.001  0.779  +0.231  31280.56
  SP  TS_02   88.19    0.00   1321.90  -275.80  +3890.00  31350.95  17.740  0.796  -1.113  31332.48
  SP  TS_03  130.43   66.81   1584.83  -104.30  +3883.00  31400.52  16.744  0.794  -0.611  31273.98
  SP  TS_04   88.63  180.00   1871.24   -71.90  +3849.00  31440.16  17.551  0.794  +0.536  31291.86
  SP  MAX01   89.93  128.31   2225.03  -337.80  +3887.00  31297.24  17.623  0.790  +0.205  31277.72
  SP  MAX02    0.00    0.00   2331.92  -522.10  +3960.00  31299.29  15.324  0.752  -1.072  31286.86
  TOR2DNS   241.97  241.97  241.97  241.97  329.03  329.03  329.06  329.06
  zpe1WHO  90.248
  zpeMSHO  90.248
  zpe2DNS   0.692
  zpeEHR   89.537
  zpeE2DT  90.229
  PF100    8.504E+04 1.701E+05 1.758E+05 2.829E+00 2.690E-01 1.671E+04
  PF100    4.975E-193 9.953E-193 1.135E-192 8.705E-02 2.690E-01 3.507E-192
  PF100    2.915E+06 1.000E+00
  PF100    1.450E-186 2.901E-186 3.307E-186
  PF100    7.490E-01 9.480E-01 6.162E+01 -5.214E+00 1.184E+01
  PF100    7.500E-01 9.480E-01 6.438E+01 -5.489E+00 1.187E+01
  PF100    7.540E-01 9.520E-01 6.310E+01 -5.358E+00 1.204E+01
  PF700    4.728E+09 1.378E+10 1.806E+10 3.638E+01 2.240E+01 1.112E+10
  PF700    3.152E-19 9.184E-19 1.221E-18 2.213E+01 2.240E+01 1.236E-18
  PF700    3.778E+08 1.000E+00
  PF700    1.191E-10 3.470E-10 4.614E-10
  PF700    1.977E+01 2.116E+01 1.137E+02 -5.846E+01 5.407E+01
  PF700    2.037E+01 2.176E+01 1.181E+02 -6.091E+01 5.491E+01
  PF700    2.045E+01 2.184E+01 1.174E+02 -6.033E+01 5.414E+01
  PF2500   1.789E+21 8.053E+21 7.567E+21 3.385E+02 2.945E+02 6.584E+21
  PF2500   2.308E+13 1.039E+14 9.802E+13 2.945E+02 2.945E+02 9.802E+13
  PF2500   9.108E+09 1.000E+00
  PF2500   2.102E+23 9.464E+23 8.928E+23
  PF2500   1.514E+02 1.563E+02 2.054E+02 -3.570E+02 8.479E+01
  PF2500   1.525E+02 1.575E+02 2.102E+02 -3.680E+02 8.487E+01
  PF2500   1.497E+02 1.547E+02 2.076E+02 -3.642E+02 8.300E+01
end_S19
#======================================#
start_S20
  minE   -342.683138
  Vmin      0.00
  Vave    682.62
  Vmax   1040.98
  V(0,0)  993.74
  const   677.78
  SP  MIN01  180.00   59.40      0.00  +129.00  +3919.00  31290.67  0.749  3.050  +0.046  31068.36
  SP  MIN02    0.00   60.70    447.51  +135.00  +3941.00  31235.35  0.752  3.041  -0.029  31051.81
  SP  TS_01  180.00    0.00    283.34  -107.60  +3920.00  31219.87  0.752  3.053  +0.046  31075.98
  SP  TS_02   92.64   49.23    773.65  -191.20  +3835.00  31043.88  0.787  3.049  +0.011  30981.50
  SP  TS_03    0.00    0.00    993.78  -173.10  +3964.00  31141.82  0.742  3.023  -0.030  31078.00
  SP  MAX01  309.03   15.06   1046.89  -171.70  +3882.00  31037.62  0.772  3.035  -0.016  31014.31
  TOR2DNS   203.28  203.34  203.34  308.37  308.37  310.19  382.44  395.55
  zpe1WHO  89.465
  zpeMSHO  89.465
  zpe2DNS   0.581
  zpeEHR   88.829
  zpeE2DT  89.410
  PF100    6.948E+04 6.973E+04 8.012E+04 1.356E+00 1.425E-01 8.422E+03
  PF100    2.096E-191 2.104E-191 3.179E-191 7.277E-02 1.425E-01 6.224E-191
  PF100    2.915E+06 1.000E+00
  PF100    6.110E-185 6.132E-185 9.265E-185
  PF100    7.310E-01 9.290E-01 6.103E+01 -5.173E+00 1.207E+01
  PF100    7.350E-01 9.340E-01 6.326E+01 -5.393E+00 1.231E+01
  PF100    7.780E-01 9.760E-01 6.178E+01 -5.202E+00 1.320E+01
  PF700    6.059E+09 9.465E+09 1.523E+10 2.028E+01 1.339E+01 1.006E+10
  PF700    7.095E-19 1.108E-18 1.854E-18 1.335E+01 1.339E+01 1.860E-18
  PF700    3.778E+08 1.000E+00
  PF700    2.681E-10 4.187E-10 7.006E-10
  PF700    2.030E+01 2.169E+01 1.150E+02 -5.881E+01 5.446E+01
  PF700    2.074E+01 2.213E+01 1.187E+02 -6.095E+01 5.484E+01
  PF700    2.040E+01 2.179E+01 1.170E+02 -6.009E+01 5.316E+01
  PF2500   3.160E+21 6.551E+21 5.281E+21 1.270E+02 1.130E+02 4.698E+21
  PF2500   4.773E+13 9.896E+13 8.066E+13 1.130E+02 1.130E+02 8.066E+13
  PF2500   9.108E+09 1.000E+00
  PF2500   4.347E+23 9.013E+23 7.347E+23
  PF2500   1.521E+02 1.571E+02 2.068E+02 -3.599E+02 8.481E+01
  PF2500   1.527E+02 1.577E+02 2.107E+02 -3.689E+02 8.485E+01
  PF2500   1.490E+02 1.540E+02 2.066E+02 -3.624E+02 8.289E+01
end_S20
#======================================#
"""

DATA_ORCA = """
#======================================#
start_S01
  minE   -299.779185
  Vmin      0.00
  Vave   3916.28
  Vmax   7301.80
  V(0,0) 3797.17
  const  3940.71
  SP  MIN01  179.92  180.01      0.00  +175.00  +3870.00  10182.99  6.336  0.729  -0.069   9792.78
  SP  MIN02    0.00  180.00    151.00  +162.00  +3877.00  10166.10  6.476  0.744  +0.382   9779.05
  SP  MIN03  180.00    0.00    304.41  +229.00  +3849.00  10227.56  7.010  0.747  -0.344   9794.19
  SP  MIN04    0.00    0.00   3796.91  +138.00  +3916.00   9957.91  6.547  0.675  -0.189   9713.01
  SP  TS_01   92.38  179.83   2591.78  -156.80  +3851.00   9936.43  7.700  0.746  +0.199   9634.60
  SP  TS_02  173.99   90.17   3577.00  -548.10  +3866.00   9746.57  6.541  0.728  -0.161   9638.47
  SP  TS_03    6.67   72.98   4574.73  -408.80  +3854.00   9697.65  6.530  0.719  +0.029   9604.74
  SP  TS_04   82.33   13.32   5557.76  -165.80  +3895.00   9775.70  7.666  0.690  -0.232   9581.98
  SP  MAX01   88.00   82.00   6396.15  -439.70  +3845.00   9483.06  7.803  0.734  +0.113   9463.90
  SP  MAX02  271.00   77.00   7307.19  -495.40  +3853.00   9468.01  7.730  0.728  -0.288   9444.89
  TOR2DNS   372.75  519.92  534.14  672.28  682.70  692.35  821.26  847.34
  zpe1WHO  29.115
  zpeMSHO  29.115
  zpe2DNS   1.066
  zpeEHR   27.999
  zpeE2DT  29.065
  PF100    1.475E+04 1.706E+04 1.793E+04 1.258E+00 5.645E-02 8.049E+02
  PF100    3.468E-60 4.011E-60 5.421E-60 5.894E-03 5.645E-02 5.192E-59
  PF100    1.652E+06 1.000E+00
  PF100    5.729E-54 6.626E-54 8.953E-54
  PF100    6.560E-01 8.540E-01 5.607E+01 -4.753E+00 9.848E+00
  PF100    7.120E-01 9.100E-01 5.692E+01 -4.782E+00 1.095E+01
  PF100    7.180E-01 9.170E-01 5.708E+01 -4.791E+00 1.113E+01
  PF700    1.104E+07 2.404E+07 3.127E+07 1.227E+01 6.022E+00 1.534E+07
  PF700    8.980E-03 1.955E-02 2.636E-02 5.705E+00 6.022E+00 2.782E-02
  PF700    2.141E+08 1.000E+00
  PF700    1.923E+06 4.187E+06 5.644E+06
  PF700    1.082E+01 1.221E+01 8.779E+01 -4.924E+01 2.593E+01
  PF700    1.112E+01 1.251E+01 8.977E+01 -5.033E+01 2.610E+01
  PF700    1.152E+01 1.291E+01 9.086E+01 -5.069E+01 2.727E+01
  PF2500   1.183E+12 3.380E+12 6.925E+12 1.552E+02 1.252E+02 5.588E+12
  PF2500   3.372E+09 9.632E+09 1.994E+10 1.252E+02 1.252E+02 1.994E+10
  PF2500   5.162E+09 1.000E+00
  PF2500   1.741E+19 4.972E+19 1.029E+20
  PF2500   6.598E+01 7.095E+01 1.281E+02 -2.492E+02 3.568E+01
  PF2500   6.735E+01 7.232E+01 1.307E+02 -2.544E+02 3.645E+01
  PF2500   6.839E+01 7.336E+01 1.325E+02 -2.580E+02 3.553E+01
end_S01
#======================================#
start_S02
  minE   -226.398436
  Vmin      0.00
  Vave   1881.82
  Vmax   4546.19
  V(0,0) 3716.40
  const  1847.42
  SP  MIN01    3.36  214.96      0.00  +128.00  +3818.00  14083.96  3.904  0.860  -0.226  13915.22
  SP  MIN02  149.25  156.16    468.58   +86.00  +3819.00  13986.83  7.218  0.857  +0.197  13857.65
  SP  TS_01    0.00  180.00     21.73   -93.40  +3826.00  14022.17  3.873  0.852  -0.196  13920.04
  SP  TS_02  180.00  180.00    578.10  -100.50  +3826.00  13860.09  7.053  0.849  +0.266  13838.36
  SP  TS_03   74.30  173.54   1372.59  -146.40  +3803.00  13916.07  6.566  0.865  +0.035  13834.50
  SP  TS_04    0.00    0.00   3724.05  -546.30  +3816.00  13941.15  3.679  0.804  -0.584  13887.08
  SP  TS_05  114.45    7.27   3808.32  -424.90  +3744.00  13862.40  7.394  0.804  -0.547  13820.94
  SP  MAX01  320.00   19.00   3815.79  -440.90  +3787.00  13851.65  4.832  0.814  -0.543  13831.50
  SP  MAX02   85.00    6.00   3864.51  -386.40  +3742.00  13814.67  7.049  0.823  -0.597  13810.00
  SP  MAX03  180.00    0.00   4553.88  -473.50  +3757.00  13786.05  7.031  0.758  -0.494  13782.45
  TOR2DNS   128.43  195.77  308.84  312.19  380.68  437.31  488.73  494.54
  zpe1WHO  40.268
  zpeMSHO  40.268
  zpe2DNS   0.367
  zpeEHR   39.786
  zpeE2DT  40.153
  PF100    8.617E+03 1.732E+04 1.084E+04 1.586E+00 3.657E-01 2.501E+03
  PF100    8.539E-85 1.717E-84 1.919E-84 2.499E-01 3.657E-01 2.809E-84
  PF100    1.207E+06 1.000E+00
  PF100    1.030E-78 2.071E-78 2.315E-78
  PF100    7.050E-01 9.040E-01 5.487E+01 -4.583E+00 1.085E+01
  PF100    7.110E-01 9.090E-01 5.631E+01 -4.722E+00 1.116E+01
  PF100    7.420E-01 9.410E-01 5.570E+01 -4.629E+00 1.127E+01
  PF700    1.077E+07 3.710E+07 2.977E+07 3.146E+01 2.433E+01 2.302E+07
  PF700    2.887E-06 9.939E-06 8.666E-06 2.416E+01 2.433E+01 8.726E-06
  PF700    1.564E+08 1.000E+00
  PF700    4.515E+02 1.555E+03 1.355E+03
  PF700    1.124E+01 1.263E+01 8.771E+01 -4.877E+01 2.740E+01
  PF700    1.175E+01 1.314E+01 9.091E+01 -5.049E+01 2.782E+01
  PF700    1.173E+01 1.312E+01 9.044E+01 -5.019E+01 2.717E+01
  PF2500   3.070E+12 1.465E+13 8.946E+12 3.309E+02 3.073E+02 8.308E+12
  PF2500   9.268E+08 4.424E+09 2.764E+09 3.073E+02 3.073E+02 2.764E+09
  PF2500   3.770E+09 1.000E+00
  PF2500   3.494E+18 1.668E+19 1.042E+19
  PF2500   7.284E+01 7.781E+01 1.321E+02 -2.524E+02 4.054E+01
  PF2500   7.360E+01 7.856E+01 1.355E+02 -2.602E+02 4.058E+01
  PF2500   7.151E+01 7.648E+01 1.337E+02 -2.577E+02 3.909E+01
end_S02
#======================================#
start_S03
  minE   -226.484421
  Vmin      0.00
  Vave   1961.05
  Vmax   3442.68
  V(0,0)    0.00
  const  1970.80
  SP  MIN01    0.00    0.00      0.00  +183.00  +3827.00  14370.70  4.478  0.785  -0.512  14072.17
  SP  MIN02  179.87  180.01    630.33   +79.00  +3899.00  14182.67  5.936  0.747  +0.311  14010.83
  SP  MIN03  197.43   81.11    700.78  +117.00  +3873.00  14232.53  6.022  0.738  +0.005  14002.06
  SP  MIN04    0.00  180.00   2185.75  +158.00  +3881.00  14220.08  3.929  0.744  -0.115  14011.79
  SP  TS_01  187.41  125.16    923.11  -263.80  +3908.00  14062.75  5.993  0.743  +0.199  13998.05
  SP  TS_02  281.41   64.53   1936.21  -151.10  +3855.00  14284.08  6.123  0.772  -0.125  14078.30
  SP  TS_03  179.97    0.00   2090.72  -402.70  +3876.00  13998.39  5.958  0.691  -0.305  13982.69
  SP  TS_04    6.11  222.39   2337.84  -286.00  +3901.00  14108.99  3.947  0.749  -0.149  14011.33
  SP  TS_05   60.30  191.25   2950.62  -161.80  +3876.00  14175.41  5.456  0.760  +0.069  14058.52
  SP  TS_06   68.77   66.23   2986.39  -171.20  +3857.00  14203.98  5.681  0.774  -0.472  14041.59
  SP  MAX01   62.00  220.00   3001.10  -251.30  +3897.00  14072.14  5.487  0.764  +0.096  14060.73
  SP  MAX02   79.00   13.00   3373.98  -378.00  +3877.00  14087.14  5.909  0.734  -0.481  14078.83
  SP  MAX03   60.00  125.00   3449.26  -317.40  +3892.00  14056.25  5.522  0.758  -0.218  14037.33
  TOR2DNS   273.50  451.81  617.12  631.12  787.45  794.39  808.48  870.88
  zpe1WHO  41.088
  zpeMSHO  41.088
  zpe2DNS   0.782
  zpeEHR   40.234
  zpeE2DT  41.016
  PF100    7.524E+03 7.546E+03 7.671E+03 1.092E+00 6.986E-02 4.906E+02
  PF100    1.205E-86 1.208E-86 1.760E-86 2.135E-02 6.986E-02 5.761E-86
  PF100    1.207E+06 1.000E+00
  PF100    1.453E-80 1.458E-80 2.124E-80
  PF100    6.500E-01 8.490E-01 5.405E+01 -4.556E+00 9.752E+00
  PF100    6.540E-01 8.530E-01 5.410E+01 -4.557E+00 1.005E+01
  PF100    6.590E-01 8.580E-01 5.418E+01 -4.560E+00 1.020E+01
  PF700    4.187E+06 1.373E+07 1.623E+07 1.869E+01 1.076E+01 9.350E+06
  PF700    6.223E-07 2.040E-06 2.539E-06 1.065E+01 1.076E+01 2.566E-06
  PF700    1.564E+08 1.000E+00
  PF700    9.734E+01 3.191E+02 3.972E+02
  PF700    1.057E+01 1.196E+01 8.488E+01 -4.746E+01 2.659E+01
  PF700    1.185E+01 1.324E+01 8.907E+01 -4.911E+01 2.766E+01
  PF700    1.190E+01 1.329E+01 8.948E+01 -4.934E+01 2.738E+01
  PF2500   7.652E+11 5.772E+12 5.472E+12 2.747E+02 2.347E+02 4.675E+12
  PF2500   1.959E+08 1.477E+09 1.421E+09 2.347E+02 2.347E+02 1.421E+09
  PF2500   3.770E+09 1.000E+00
  PF2500   7.385E+17 5.571E+18 5.358E+18
  PF2500   7.188E+01 7.685E+01 1.289E+02 -2.455E+02 4.056E+01
  PF2500   7.391E+01 7.887E+01 1.338E+02 -2.555E+02 4.073E+01
  PF2500   7.168E+01 7.665E+01 1.328E+02 -2.553E+02 3.898E+01
end_S03
#======================================#
start_S04
  minE   -153.222817
  Vmin      0.00
  Vave    903.06
  Vmax   2056.96
  V(0,0) 2056.96
  const   874.59
  SP  MIN01   62.45   54.97      0.00  +269.00  +3858.00  18727.87  0.743  2.613  -0.151  18422.94
  SP  MIN02  180.03   59.83     52.45  +254.00  +3866.00  18728.98  0.743  2.619  +0.034  18443.09
  SP  TS_01  121.59   55.84    509.84  -323.40  +3896.00  18577.39  0.741  2.619  -0.028  18430.59
  SP  TS_02    0.00   60.14    577.44  -346.10  +3889.00  18611.20  0.710  2.611  -0.213  18456.19
  SP  TS_03  180.00    0.00   1204.92  -257.90  +3864.00  18634.19  0.746  2.605  +0.032  18473.76
  SP  TS_04  292.51    4.73   1249.25  -267.00  +3858.00  18631.77  0.746  2.593  -0.144  18453.76
  SP  MAX01  238.00    3.00   1689.30  -328.80  +3894.00  18475.99  0.744  2.603  -0.030  18462.03
  SP  MAX02    0.00    0.00   2079.30  -394.90  +3902.00  18502.39  0.702  2.573  -0.216  18495.91
  TOR2DNS   281.45  281.45  281.45  282.67  282.67  282.67  312.30  312.30
  zpe1WHO  53.546
  zpeMSHO  53.546
  zpe2DNS   0.805
  zpeEHR   52.674
  zpeE2DT  53.479
  PF100    3.487E+03 8.558E+03 8.731E+03 2.770E+00 1.481E-01 4.666E+02
  PF100    3.319E-114 8.144E-114 1.165E-113 4.829E-02 1.481E-01 3.571E-113
  PF100    8.106E+05 1.000E+00
  PF100    2.690E-108 6.602E-108 9.441E-108
  PF100    6.230E-01 8.220E-01 5.146E+01 -4.325E+00 9.135E+00
  PF100    6.520E-01 8.510E-01 5.572E+01 -4.721E+00 9.354E+00
  PF100    6.640E-01 8.630E-01 5.370E+01 -4.507E+00 9.795E+00
  PF700    1.263E+06 3.710E+06 4.698E+06 1.842E+01 1.046E+01 2.667E+06
  PF700    2.422E-11 7.112E-11 9.452E-11 1.033E+01 1.046E+01 9.569E-11
  PF700    1.051E+08 1.000E+00
  PF700    2.545E-03 7.474E-03 9.933E-03
  PF700    1.036E+01 1.175E+01 8.141E+01 -4.524E+01 2.724E+01
  PF700    1.041E+01 1.180E+01 8.581E+01 -4.826E+01 2.723E+01
  PF700    1.037E+01 1.177E+01 8.404E+01 -4.706E+01 2.645E+01
  PF2500   3.891E+11 1.172E+12 1.017E+12 1.088E+02 9.253E+01 8.648E+11
  PF2500   8.112E+06 2.443E+07 2.149E+07 9.253E+01 9.253E+01 2.149E+07
  PF2500   2.533E+09 1.000E+00
  PF2500   2.055E+16 6.188E+16 5.444E+16
  PF2500   7.759E+01 8.256E+01 1.291E+02 -2.402E+02 4.534E+01
  PF2500   7.764E+01 8.261E+01 1.335E+02 -2.511E+02 4.534E+01
  PF2500   7.481E+01 7.978E+01 1.299E+02 -2.449E+02 4.351E+01
end_S04
#======================================#
start_S05
  minE   -264.168822
  Vmin      0.00
  Vave   4438.71
  Vmax   8086.15
  V(0,0)    0.00
  const  4559.98
  SP  MIN01    0.00    0.00      0.00  +159.00  +3881.00  15934.94  8.496  0.763  +0.433  15544.76
  SP  MIN02  179.99    0.00    237.69  +157.00  +3886.00  15922.90  8.284  0.745  -0.028  15548.46
  SP  MIN03    0.00  180.01   2958.96  +136.00  +3928.00  15791.47  8.564  0.695  -0.263  15518.84
  SP  MIN04  159.79  159.84   3867.80  +135.00  +3930.00  15732.14  8.281  0.712  -0.389  15463.92
  SP  TS_01   92.48    0.21   3730.41  -175.20  +3854.00  15695.86  9.511  0.763  +0.231  15379.92
  SP  TS_02  180.00  180.00   3995.54  -136.10  +3960.00  15661.00  7.715  0.684  -0.346  15499.55
  SP  TS_03  169.42  107.10   4324.31  -374.60  +3859.00  15506.44  8.442  0.759  -0.302  15395.43
  SP  TS_04    3.77  258.14   4330.67  -469.80  +3855.00  15471.80  8.567  0.741  +0.038  15380.91
  SP  TS_05   90.15  182.63   5734.21  -152.80  +3915.00  15608.26  9.631  0.721  -0.356  15366.33
  SP  MAX01  261.00   96.00   7978.12  -534.60  +3859.00  15223.92  9.565  0.742  +0.097  15206.76
  SP  MAX02   90.00   97.00   8093.57  -537.00  +3856.00  15223.37  9.623  0.745  -0.290  15205.39
  TOR2DNS   370.39  522.48  585.23  672.97  727.69  821.93  869.17  945.40
  zpe1WHO  45.560
  zpeMSHO  45.560
  zpe2DNS   1.059
  zpeEHR   44.445
  zpeE2DT  45.504
  PF100    1.553E+04 1.613E+04 1.686E+04 1.179E+00 5.478E-02 7.833E+02
  PF100    4.183E-96 4.343E-96 6.035E-96 5.717E-03 5.478E-02 5.783E-95
  PF100    1.586E+06 1.000E+00
  PF100    6.633E-90 6.888E-90 9.572E-90
  PF100    6.610E-01 8.590E-01 5.614E+01 -4.755E+00 9.895E+00
  PF100    6.840E-01 8.830E-01 5.645E+01 -4.762E+00 1.063E+01
  PF100    6.880E-01 8.870E-01 5.658E+01 -4.771E+00 1.074E+01
  PF700    1.894E+07 3.074E+07 3.692E+07 9.522E+00 4.695E+00 1.821E+07
  PF700    1.130E-07 1.834E-07 2.294E-07 4.447E+00 4.695E+00 2.422E-07
  PF700    2.056E+08 1.000E+00
  PF700    2.324E+01 3.772E+01 4.717E+01
  PF700    1.226E+01 1.365E+01 9.083E+01 -4.994E+01 3.116E+01
  PF700    1.254E+01 1.394E+01 9.221E+01 -5.061E+01 3.158E+01
  PF700    1.282E+01 1.421E+01 9.297E+01 -5.086E+01 3.245E+01
  PF2500   3.412E+13 9.220E+13 1.482E+14 1.314E+02 1.062E+02 1.197E+14
  PF2500   3.550E+09 9.592E+09 1.560E+10 1.062E+02 1.062E+02 1.560E+10
  PF2500   4.956E+09 1.000E+00
  PF2500   1.759E+19 4.754E+19 7.729E+19
  PF2500   8.316E+01 8.813E+01 1.415E+02 -2.657E+02 4.623E+01
  PF2500   8.643E+01 9.140E+01 1.448E+02 -2.707E+02 4.788E+01
  PF2500   8.662E+01 9.159E+01 1.458E+02 -2.730E+02 4.653E+01
end_S05
#======================================#
start_S06
  minE   -264.149214
  Vmin      0.00
  Vave   6120.78
  Vmax   8844.25
  V(0,0)    0.00
  const  6214.28
  SP  MIN01    0.00    0.00      0.00  +293.00  +3607.00  16135.99  6.527  0.803  -0.462  15549.77
  SP  MIN02  180.00  180.00   4155.09  +182.00  +3930.00  15706.57  6.891  0.722  +0.276  15495.94
  SP  MIN03  180.00    0.00   4328.04  +166.00  +3923.00  15774.52  6.908  0.683  -0.006  15494.38
  SP  MIN04    0.00  180.00   5551.39  +150.00  +3918.00  15758.79  5.525  0.717  +0.092  15529.82
  SP  TS_01  182.83   74.43   5371.86  -356.80  +3855.00  15493.60  6.847  0.758  +0.139  15367.51
  SP  TS_02   11.17  246.83   6244.49  -385.20  +3850.00  15501.67  5.883  0.756  -0.006  15422.61
  SP  TS_03  261.79    9.93   7284.58  -202.40  +3888.00  15534.14  8.509  0.722  -0.174  15337.98
  SP  TS_04   87.84  201.76   7732.53  -170.90  +3915.00  15410.15  8.666  0.735  +0.335  15315.95
  SP  MAX01  266.00   84.00   8045.06  -279.70  +3842.00  15263.11  8.768  0.777  +0.242  15231.06
  SP  MAX02   89.00   77.00   8849.88  -375.50  +3846.00  15246.72  8.640  0.770  -0.285  15200.47
  TOR2DNS   468.07  710.27  950.67 1149.03 1189.21 1388.23 1425.80 1625.50
  zpe1WHO  46.135
  zpeMSHO  46.135
  zpe2DNS   1.338
  zpeEHR   44.459
  zpeE2DT  45.797
  PF100    1.329E+04 1.329E+04 1.405E+04 1.032E+00 2.811E-02 3.828E+02
  PF100    1.983E-97 1.983E-97 1.148E-96 1.227E-03 2.811E-02 2.630E-95
  PF100    1.586E+06 1.000E+00
  PF100    3.144E-91 3.144E-91 1.820E-90
  PF100    6.220E-01 8.200E-01 5.544E+01 -4.724E+00 9.115E+00
  PF100    6.220E-01 8.200E-01 5.544E+01 -4.724E+00 9.115E+00
  PF100    6.260E-01 8.250E-01 5.560E+01 -4.735E+00 9.351E+00
  PF700    8.273E+06 8.305E+06 1.107E+07 3.459E+00 1.439E+00 4.605E+06
  PF700    3.265E-08 3.278E-08 5.571E-08 1.322E+00 1.439E+00 6.065E-08
  PF700    2.056E+08 1.000E+00
  PF700    6.714E+00 6.739E+00 1.145E+01
  PF700    1.164E+01 1.303E+01 8.830E+01 -4.878E+01 3.052E+01
  PF700    1.168E+01 1.307E+01 8.838E+01 -4.879E+01 3.109E+01
  PF700    1.207E+01 1.346E+01 8.950E+01 -4.919E+01 3.206E+01
  PF2500   1.014E+13 3.377E+13 6.829E+13 5.053E+01 3.860E+01 5.217E+13
  PF2500   9.394E+08 3.130E+09 6.774E+09 3.860E+01 3.860E+01 6.774E+09
  PF2500   4.956E+09 1.000E+00
  PF2500   4.655E+18 1.551E+19 3.357E+19
  PF2500   8.239E+01 8.736E+01 1.388E+02 -2.597E+02 4.627E+01
  PF2500   9.140E+01 9.637E+01 1.448E+02 -2.657E+02 4.920E+01
  PF2500   9.180E+01 9.676E+01 1.464E+02 -2.692E+02 4.760E+01
end_S06
#======================================#
start_S07
  minE   -289.184021
  Vmin      0.00
  Vave   2473.93
  Vmax   3952.27
  V(0,0)    0.00
  const  2531.49
  SP  MIN01    0.00    0.00      0.00  +138.00  +3904.00  18383.85  0.745  11.408  +0.510  18108.19
  SP  MIN02  200.33  150.58    932.99  +103.00  +3882.00  18397.04  0.783  13.485  -0.479  18088.23
  SP  MIN03    7.10  240.54   1731.00   +78.00  +3900.00  18365.48  0.738  13.509  +0.226  18124.84
  SP  TS_01  180.00  180.00   1100.67   -84.40  +3892.00  18340.82  0.778  12.706  -0.545  18077.62
  SP  TS_02    2.53  285.16   2003.36   -90.60  +3900.00  18343.81  0.746  13.430  +0.372  18131.38
  SP  TS_03  180.00    0.00   2058.45  -223.90  +3902.00  18154.45  0.721  11.441  -0.316  18076.95
  SP  TS_04   73.64  236.78   2568.51  -364.70  +3871.00  18101.00  0.774  13.667  +0.129  18035.06
  SP  TS_05    0.00  180.00   2610.65  -115.50  +3897.00  18274.01  0.731  11.748  +0.008  18053.80
  SP  TS_06  242.09   64.48   3226.94  -115.90  +3833.00  18128.64  0.763  13.252  -0.046  18024.48
  SP  TS_07   85.79  125.48   3530.69  -340.80  +3839.00  18062.46  0.769  13.326  -0.338  18008.02
  SP  MAX01  271.00   67.00   3270.61  -229.40  +3835.00  18031.03  0.773  13.361  +0.165  18020.79
  SP  MAX02  101.00   73.00   3937.81  -231.20  +3832.00  18027.53  0.774  13.417  -0.305  18017.69
  SP  MAX03   76.00  172.00   3961.08  -399.20  +3856.00  17988.33  0.768  11.893  -0.224  17972.34
  TOR2DNS   253.68  378.45  501.91  622.24  623.99  744.38  747.84  863.19
  zpe1WHO  52.562
  zpeMSHO  52.562
  zpe2DNS   0.725
  zpeEHR   51.774
  zpeE2DT  52.499
  PF100    1.939E+04 1.939E+04 2.044E+04 1.206E+00 9.839E-02 1.667E+03
  PF100    2.604E-111 2.604E-111 3.766E-111 3.135E-02 9.839E-02 1.182E-110
  PF100    1.720E+06 1.000E+00
  PF100    4.480E-105 4.480E-105 6.479E-105
  PF100    6.770E-01 8.760E-01 5.691E+01 -4.815E+00 1.033E+01
  PF100    6.770E-01 8.760E-01 5.691E+01 -4.815E+00 1.033E+01
  PF100    6.880E-01 8.860E-01 5.712E+01 -4.826E+00 1.054E+01
  PF700    5.134E+07 7.296E+07 1.058E+08 1.360E+01 8.223E+00 6.392E+07
  PF700    1.996E-09 2.837E-09 4.302E-09 8.077E+00 8.223E+00 4.380E-09
  PF700    2.230E+08 1.000E+00
  PF700    4.452E-01 6.326E-01 9.594E-01
  PF700    1.345E+01 1.484E+01 9.468E+01 -5.144E+01 3.422E+01
  PF700    1.442E+01 1.581E+01 9.677E+01 -5.192E+01 3.683E+01
  PF700    1.492E+01 1.631E+01 9.822E+01 -5.244E+01 3.743E+01
  PF2500   5.278E+14 1.945E+15 3.037E+15 3.052E+02 2.637E+02 2.624E+15
  PF2500   1.342E+10 4.942E+10 7.817E+10 2.637E+02 2.637E+02 7.817E+10
  PF2500   5.376E+09 1.000E+00
  PF2500   7.212E+19 2.657E+20 4.202E+20
  PF2500   9.264E+01 9.761E+01 1.509E+02 -2.797E+02 5.161E+01
  PF2500   9.548E+01 1.004E+02 1.547E+02 -2.862E+02 5.193E+01
  PF2500   9.387E+01 9.884E+01 1.549E+02 -2.884E+02 5.017E+01
end_S07
#======================================#
start_S08
  minE   -190.855295
  Vmin      0.00
  Vave   1413.28
  Vmax   2211.99
  V(0,0) 1103.34
  const  1429.10
  SP  MIN01    0.22   60.39      0.00  +226.00  +3431.00  20149.04  4.074  2.938  -1.498  19909.50
  SP  MIN02  150.06   54.20   1140.83   +67.00  +3433.00  19933.22  7.037  2.553  -0.610  19809.04
  SP  TS_01    0.00    0.00   1103.74  -232.90  +3432.00  20000.75  3.802  2.907  -1.509  19920.37
  SP  TS_02  180.01   61.04   1191.31   -61.20  +3435.00  19882.49  6.982  2.453  -0.592  19799.57
  SP  TS_03  180.00    0.00   1523.81  -153.50  +3436.00  19854.71  6.875  2.428  -0.683  19812.44
  SP  TS_04   75.12   63.02   1663.40  -129.00  +3422.00  19883.49  6.249  2.887  -1.071  19794.91
  SP  MAX01   90.00    2.00   2215.16  -187.90  +3424.00  19791.86  6.580  2.843  -0.977  19787.57
  TOR2DNS   221.41  221.41  221.41  426.27  426.27  426.27  446.06  446.06
  zpe1WHO  57.609
  zpeMSHO  57.609
  zpe2DNS   0.633
  zpeEHR   56.924
  zpeE2DT  57.557
  PF100    7.761E+03 7.761E+03 8.104E+03 1.100E+00 9.508E-02 7.005E+02
  PF100    9.732E-123 9.732E-123 1.319E-122 4.549E-02 9.508E-02 2.757E-122
  PF100    1.147E+06 1.000E+00
  PF100    1.117E-116 1.117E-116 1.514E-116
  PF100    6.520E-01 8.500E-01 5.403E+01 -4.553E+00 1.003E+01
  PF100    6.520E-01 8.500E-01 5.622E+01 -4.771E+00 1.003E+01
  PF100    6.650E-01 8.640E-01 5.425E+01 -4.561E+00 1.038E+01
  PF700    7.796E+06 1.602E+07 1.986E+07 1.544E+01 9.843E+00 1.266E+07
  PF700    8.052E-12 1.655E-11 2.129E-11 9.796E+00 9.843E+00 2.139E-11
  PF700    1.487E+08 1.000E+00
  PF700    1.198E-03 2.461E-03 3.167E-03
  PF700    1.195E+01 1.334E+01 8.798E+01 -4.825E+01 3.139E+01
  PF700    1.349E+01 1.489E+01 9.381E+01 -5.078E+01 3.381E+01
  PF700    1.320E+01 1.459E+01 9.163E+01 -4.955E+01 3.206E+01
  PF2500   2.193E+13 1.353E+14 8.564E+13 2.130E+02 1.875E+02 7.539E+13
  PF2500   2.018E+08 1.245E+09 7.964E+08 1.875E+02 1.875E+02 7.964E+08
  PF2500   3.585E+09 1.000E+00
  PF2500   7.237E+17 4.464E+18 2.855E+18
  PF2500   8.883E+01 9.380E+01 1.423E+02 -2.619E+02 5.113E+01
  PF2500   9.147E+01 9.644E+01 1.491E+02 -2.764E+02 5.127E+01
  PF2500   8.776E+01 9.273E+01 1.446E+02 -2.687E+02 4.932E+01
end_S08
#======================================#
start_S09
  minE   -190.854256
  Vmin      0.00
  Vave   1319.63
  Vmax   2495.71
  V(0,0) 1559.84
  const  1333.28
  SP  MIN01   10.36   69.50      0.00  +202.00  +3861.00  20065.81  5.262  0.791  -0.462  19766.33
  SP  MIN02  229.56   60.16    390.66  +110.00  +3859.00  20010.33  7.360  0.748  -0.158  19760.47
  SP  MIN03    0.00  180.01    496.89  +153.00  +3885.00  19960.95  5.069  0.746  -0.092  19760.90
  SP  MIN04  140.05  179.33   1044.92  +103.00  +3871.00  19922.67  7.243  0.765  +0.284  19732.49
  SP  MIN05  128.04   63.59   1091.45   +84.00  +3859.00  19931.31  7.285  0.769  -0.364  19745.01
  SP  TS_01    2.25  225.50    639.33  -260.90  +3902.00  19876.03  5.100  0.756  -0.163  19764.65
  SP  TS_02  161.88   72.09   1213.04  -102.60  +3862.00  19849.63  7.205  0.759  -0.204  19688.52
  SP  TS_03  180.00  180.00   1256.93   -99.80  +3881.00  19818.63  7.117  0.759  +0.318  19683.25
  SP  TS_04  118.79   14.40   1331.77  -308.60  +3885.00  19842.50  7.290  0.741  -0.483  19787.90
  SP  TS_05  216.63  133.54   1335.50  -321.60  +3896.00  19791.00  7.249  0.755  +0.247  19729.88
  SP  TS_06  138.00  118.58   1405.30  -284.80  +3894.00  19786.37  7.332  0.763  -0.001  19724.74
  SP  TS_07   77.48   49.70   1412.10  -121.20  +3862.00  19924.24  6.997  0.780  -0.556  19780.98
  SP  TS_08    0.00    0.00   1560.90  -472.90  +3937.00  19834.28  4.737  0.728  -0.524  19774.26
  SP  TS_09  307.72   50.84   1738.68  -175.90  +3883.00  19879.77  6.098  0.758  -0.324  19762.33
  SP  TS_10   69.09  178.96   2015.44  -154.50  +3848.00  19922.46  6.810  0.770  +0.089  19762.11
  SP  MAX01  175.00  127.00   1545.54  -292.20  +3899.00  19694.82  7.185  0.754  +0.163  19680.67
  SP  MAX02  314.00   28.00   1770.94  -256.80  +3903.00  19787.57  5.925  0.747  -0.445  19777.24
  SP  MAX03  180.00    0.00   2072.06  -401.10  +3903.00  19722.63  7.122  0.708  -0.399  19714.78
  SP  MAX04   66.00  128.00   2409.39  -350.00  +3886.00  19764.27  6.791  0.759  -0.184  19749.40
  SP  MAX05   69.00  236.00   2506.62  -353.70  +3891.00  19766.08  6.627  0.774  +0.087  19748.77
  TOR2DNS   276.85  276.85  459.60  459.62  609.55  611.53  624.24  624.24
  zpe1WHO  57.371
  zpeMSHO  57.371
  zpe2DNS   0.792
  zpeEHR   56.515
  zpeE2DT  57.306
  PF100    8.377E+03 1.694E+04 1.734E+04 2.197E+00 1.375E-01 1.086E+03
  PF100    3.479E-122 7.037E-122 9.976E-122 4.091E-02 1.375E-01 3.353E-121
  PF100    1.147E+06 1.000E+00
  PF100    3.992E-116 8.074E-116 1.145E-115
  PF100    6.440E-01 8.430E-01 5.410E+01 -4.568E+00 9.690E+00
  PF100    6.550E-01 8.540E-01 5.562E+01 -4.708E+00 1.029E+01
  PF100    6.670E-01 8.660E-01 5.578E+01 -4.712E+00 1.059E+01
  PF700    7.933E+06 4.714E+07 5.840E+07 3.718E+01 2.127E+01 3.340E+07
  PF700    9.723E-12 5.777E-11 7.498E-11 2.105E+01 2.127E+01 7.576E-11
  PF700    1.487E+08 1.000E+00
  PF700    1.446E-03 8.593E-03 1.115E-02
  PF700    1.212E+01 1.351E+01 8.827E+01 -4.827E+01 3.192E+01
  PF700    1.335E+01 1.474E+01 9.356E+01 -5.075E+01 3.325E+01
  PF700    1.321E+01 1.460E+01 9.379E+01 -5.105E+01 3.193E+01
  PF2500   2.617E+13 3.525E+14 2.549E+14 4.050E+02 3.453E+02 2.173E+14
  PF2500   2.526E+08 3.404E+09 2.493E+09 3.453E+02 3.453E+02 2.493E+09
  PF2500   3.585E+09 1.000E+00
  PF2500   9.059E+17 1.220E+19 8.939E+18
  PF2500   8.921E+01 9.418E+01 1.428E+02 -2.628E+02 5.111E+01
  PF2500   9.117E+01 9.613E+01 1.487E+02 -2.757E+02 5.123E+01
  PF2500   8.778E+01 9.275E+01 1.467E+02 -2.741E+02 4.929E+01
end_S09
#======================================#
start_S10
  minE   -190.869302
  Vmin      0.00
  Vave   1481.16
  Vmax   2492.80
  V(0,0)    0.00
  const  1523.91
  SP  MIN01    0.00    0.00      0.00  +182.00  +3893.00  19934.83  0.729  2.932  +0.068  19616.78
  SP  MIN02  180.00    0.00   1128.10  +133.00  +3917.00  19793.23  0.728  2.924  -0.103  19607.18
  SP  TS_01    0.00   60.17    657.77  -185.60  +3893.00  19847.06  0.730  2.927  +0.068  19612.96
  SP  TS_02   88.62    1.23   1735.82  -333.50  +3844.00  19625.08  0.765  2.930  -0.015  19513.52
  SP  TS_03  166.55   60.60   2104.54  -227.50  +3921.00  19641.94  0.725  2.909  -0.100  19591.41
  SP  MAX01  180.00   60.00   2105.86  -229.50  +3930.00  19604.76  0.722  2.908  -0.103  19597.31
  SP  MAX02   90.00   62.00   2495.43  -315.30  +3838.00  19519.81  0.768  2.924  -0.018  19509.29
  TOR2DNS   294.46  294.47  294.47  459.68  459.68  459.73  610.16  611.00
  zpe1WHO  56.997
  zpeMSHO  56.997
  zpe2DNS   0.842
  zpeEHR   56.087
  zpeE2DT  56.929
  PF100    8.557E+03 8.557E+03 8.932E+03 1.109E+00 6.545E-02 5.273E+02
  PF100    2.339E-121 2.339E-121 3.429E-121 1.603E-02 6.545E-02 1.400E-120
  PF100    1.147E+06 1.000E+00
  PF100    2.684E-115 2.684E-115 3.934E-115
  PF100    6.430E-01 8.420E-01 5.414E+01 -4.572E+00 9.504E+00
  PF100    6.430E-01 8.420E-01 5.632E+01 -4.790E+00 9.504E+00
  PF100    6.540E-01 8.530E-01 5.433E+01 -4.581E+00 9.823E+00
  PF700    9.086E+06 1.170E+07 1.537E+07 8.376E+00 4.665E+00 8.563E+06
  PF700    1.458E-11 1.877E-11 2.589E-11 4.573E+00 4.665E+00 2.641E-11
  PF700    1.487E+08 1.000E+00
  PF700    2.168E-03 2.792E-03 3.850E-03
  PF700    1.238E+01 1.377E+01 8.890E+01 -4.846E+01 3.224E+01
  PF700    1.308E+01 1.447E+01 9.259E+01 -5.034E+01 3.403E+01
  PF700    1.305E+01 1.444E+01 9.091E+01 -4.920E+01 3.303E+01
  PF2500   3.545E+13 8.834E+13 7.524E+13 8.361E+01 7.058E+01 6.351E+13
  PF2500   3.691E+08 9.196E+08 7.940E+08 7.058E+01 7.058E+01 7.940E+08
  PF2500   3.585E+09 1.000E+00
  PF2500   1.323E+18 3.297E+18 2.847E+18
  PF2500   8.959E+01 9.455E+01 1.435E+02 -2.643E+02 5.112E+01
  PF2500   9.151E+01 9.647E+01 1.483E+02 -2.743E+02 5.132E+01
  PF2500   8.842E+01 9.339E+01 1.446E+02 -2.680E+02 4.942E+01
end_S10
#======================================#
start_S11
  minE   -190.862329
  Vmin      0.00
  Vave    610.39
  Vmax   1156.80
  V(0,0)    0.00
  const   629.77
  SP  MIN01    0.00    0.00      0.00   +87.00  +3920.00  19970.50  0.714  2.887  -0.106  19757.15
  SP  MIN02  209.92    0.57    395.49  +144.00  +3896.00  19857.43  0.743  2.928  +0.072  19707.35
  SP  TS_01    0.00   60.72    169.21  -100.00  +3905.00  19902.86  0.723  2.901  -0.099  19729.31
  SP  TS_02  180.00    0.00    413.27  -118.50  +3928.00  19800.55  0.729  2.926  +0.086  19724.07
  SP  TS_03   77.11    8.29    715.05  -255.30  +3847.00  19745.94  0.768  2.930  -0.026  19646.65
  SP  TS_04  179.98   59.49    759.38  -140.20  +3927.00  19750.83  0.726  2.919  +0.087  19714.87
  SP  MAX01   81.00   67.00   1159.05  -291.90  +3835.00  19655.97  0.769  2.918  -0.018  19640.07
  TOR2DNS   198.93  199.32  199.32  273.00  273.00  279.21  323.31  344.55
  zpe1WHO  57.099
  zpeMSHO  57.099
  zpe2DNS   0.569
  zpeEHR   56.489
  zpeE2DT  57.057
  PF100    1.168E+04 1.204E+04 1.508E+04 1.578E+00 1.887E-01 1.802E+03
  PF100    1.911E-121 1.970E-121 3.037E-121 9.020E-02 1.887E-01 6.352E-121
  PF100    1.147E+06 1.000E+00
  PF100    2.193E-115 2.260E-115 3.484E-115
  PF100    7.170E-01 9.160E-01 5.550E+01 -4.634E+00 1.070E+01
  PF100    7.410E-01 9.400E-01 5.798E+01 -4.858E+00 1.169E+01
  PF100    7.520E-01 9.510E-01 5.636E+01 -4.685E+00 1.137E+01
  PF700    2.265E+07 5.319E+07 4.733E+07 2.145E+01 1.433E+01 3.161E+07
  PF700    3.376E-11 7.929E-11 7.269E-11 1.425E+01 1.433E+01 7.307E-11
  PF700    1.487E+08 1.000E+00
  PF700    5.022E-03 1.179E-02 1.081E-02
  PF700    1.246E+01 1.385E+01 9.083E+01 -4.973E+01 3.197E+01
  PF700    1.303E+01 1.442E+01 9.553E+01 -5.245E+01 3.227E+01
  PF700    1.238E+01 1.378E+01 9.219E+01 -5.076E+01 3.065E+01
  PF2500   8.879E+13 2.918E+14 1.130E+14 1.256E+02 1.120E+02 1.008E+14
  PF2500   9.056E+08 2.976E+09 1.162E+09 1.120E+02 1.120E+02 1.162E+09
  PF2500   3.585E+09 1.000E+00
  PF2500   3.247E+18 1.067E+19 4.166E+18
  PF2500   8.952E+01 9.449E+01 1.453E+02 -2.689E+02 5.110E+01
  PF2500   9.027E+01 9.523E+01 1.502E+02 -2.802E+02 5.113E+01
  PF2500   8.626E+01 9.123E+01 1.445E+02 -2.701E+02 4.918E+01
end_S11
#======================================#
start_S12
  minE   -228.516006
  Vmin      0.00
  Vave   2089.22
  Vmax   4296.87
  V(0,0)  670.79
  const  2129.84
  SP  MIN01  179.99    0.00      0.00  +138.00  +3395.00  21229.01  6.778  3.028  -0.785  21056.30
  SP  MIN02    0.16    0.00    671.15  +157.00  +3411.00  21194.97  6.819  2.955  +0.395  21034.99
  SP  TS_01  180.00   60.20    373.33  -129.40  +3391.00  21171.75  6.643  3.030  -0.780  21070.26
  SP  TS_02    0.00   60.18   1127.22  -148.50  +3409.00  21140.84  6.882  2.943  +0.378  21059.38
  SP  TS_03   91.77    2.27   3571.73  -173.20  +3389.00  20988.35  8.061  3.023  -0.227  20888.22
  SP  MAX01   93.00   62.00   4312.02  -197.80  +3385.00  20901.85  8.015  3.007  -0.240  20896.67
  TOR2DNS   151.53  151.55  151.55  274.32  274.32  274.85  322.93  322.93
  zpe1WHO  60.697
  zpeMSHO  60.697
  zpe2DNS   0.433
  zpeEHR   60.203
  zpeE2DT  60.636
  PF100    1.952E+04 1.952E+04 2.112E+04 1.338E+00 2.161E-01 3.412E+03
  PF100    4.370E-129 4.371E-129 6.412E-129 1.512E-01 2.161E-01 9.168E-129
  PF100    1.521E+06 1.000E+00
  PF100    6.647E-123 6.648E-123 9.753E-123
  PF100    7.110E-01 9.090E-01 5.701E+01 -4.792E+00 1.126E+01
  PF100    7.110E-01 9.100E-01 5.920E+01 -5.010E+00 1.128E+01
  PF100    7.390E-01 9.380E-01 5.745E+01 -4.808E+00 1.181E+01
  PF700    9.313E+07 1.201E+08 1.531E+08 1.937E+01 1.425E+01 1.126E+08
  PF700    1.045E-11 1.348E-11 1.794E-11 1.419E+01 1.425E+01 1.802E-11
  PF700    1.972E+08 1.000E+00
  PF700    2.060E-03 2.658E-03 3.537E-03
  PF700    1.403E+01 1.542E+01 9.645E+01 -5.209E+01 3.599E+01
  PF700    1.445E+01 1.584E+01 9.973E+01 -5.398E+01 3.661E+01
  PF700    1.438E+01 1.577E+01 9.794E+01 -5.278E+01 3.641E+01
  PF2500   3.318E+15 5.844E+15 6.488E+15 2.005E+02 1.837E+02 5.947E+15
  PF2500   1.640E+10 2.889E+10 3.247E+10 1.837E+02 1.837E+02 3.247E+10
  PF2500   4.753E+09 1.000E+00
  PF2500   7.795E+19 1.373E+20 1.543E+20
  PF2500   1.009E+02 1.058E+02 1.576E+02 -2.882E+02 5.699E+01
  PF2500   1.017E+02 1.067E+02 1.613E+02 -2.965E+02 5.706E+01
  PF2500   1.006E+02 1.055E+02 1.588E+02 -2.916E+02 5.584E+01
end_S12
#======================================#
start_S13
  minE   -228.498044
  Vmin      0.00
  Vave   1513.20
  Vmax   3051.53
  V(0,0) 1304.57
  const  1545.22
  SP  MIN01  180.00    0.00      0.00  +102.00  +3923.00  21325.11  8.837  0.705  -0.066  21113.23
  SP  MIN02   43.47    9.05    109.30  +136.00  +3870.00  21413.14  9.107  0.764  -0.464  21087.31
  SP  MIN03  179.73  160.17    157.36  +124.00  +3914.00  21236.61  8.724  0.741  +0.323  21096.30
  SP  MIN04    0.00  180.02   1049.97   +82.00  +3924.00  21284.21  7.072  0.728  +0.127  21154.62
  SP  TS_01  180.00  180.00    158.02   -63.60  +3931.00  21183.58  8.721  0.733  +0.336  21103.47
  SP  TS_02  185.40   74.19    621.77  -267.80  +3845.00  21116.70  8.721  0.774  +0.135  21017.99
  SP  TS_03  249.24    7.93   1028.46  -108.20  +3882.00  21255.98  10.502  0.748  -0.186  21041.31
  SP  TS_04    0.00    0.00   1304.56  -262.80  +3980.00  21326.72  6.409  0.721  -0.443  21176.98
  SP  TS_05   14.33   95.01   1372.59  -256.30  +3822.00  21153.46  7.560  0.787  -0.217  21068.25
  SP  TS_06   90.23  188.08   2594.19  -106.90  +3921.00  21062.39  10.785  0.740  +0.348  21000.90
  SP  MAX01   85.00  259.00   2815.86  -267.60  +3843.00  20942.47  10.840  0.781  +0.296  20933.61
  SP  MAX02   97.00   92.00   3053.99  -312.70  +3838.00  20932.66  10.806  0.780  -0.217  20914.41
  TOR2DNS   200.72  269.47  304.94  381.61  389.21  389.21  406.71  407.43
  zpe1WHO  60.972
  zpeMSHO  60.972
  zpe2DNS   0.574
  zpeEHR   60.366
  zpeE2DT  60.940
  PF100    2.151E+04 3.999E+04 3.853E+04 2.017E+00 2.268E-01 4.331E+03
  PF100    1.209E-129 2.247E-129 2.542E-129 1.123E-01 2.268E-01 5.131E-129
  PF100    1.521E+06 1.000E+00
  PF100    1.838E-123 3.417E-123 3.866E-123
  PF100    7.260E-01 9.250E-01 5.736E+01 -4.811E+00 1.128E+01
  PF100    8.470E-01 1.046E+00 5.980E+01 -4.934E+00 1.264E+01
  PF100    8.350E-01 1.034E+00 5.961E+01 -4.927E+00 1.264E+01
  PF700    1.148E+08 5.209E+08 5.044E+08 4.400E+01 2.943E+01 3.374E+08
  PF700    1.057E-11 4.797E-11 4.753E-11 2.913E+01 2.943E+01 4.802E-11
  PF700    1.972E+08 1.000E+00
  PF700    2.085E-03 9.460E-03 9.371E-03
  PF700    1.428E+01 1.567E+01 9.722E+01 -5.238E+01 3.664E+01
  PF700    1.469E+01 1.608E+01 1.008E+02 -5.449E+01 3.703E+01
  PF700    1.470E+01 1.609E+01 1.008E+02 -5.444E+01 3.679E+01
  PF2500   4.921E+15 2.979E+16 2.314E+16 4.323E+02 3.851E+02 2.062E+16
  PF2500   2.302E+10 1.393E+11 1.089E+11 3.851E+02 3.851E+02 1.089E+11
  PF2500   4.753E+09 1.000E+00
  PF2500   1.094E+20 6.622E+20 5.178E+20
  PF2500   1.012E+02 1.061E+02 1.585E+02 -2.902E+02 5.685E+01
  PF2500   1.019E+02 1.069E+02 1.624E+02 -2.991E+02 5.694E+01
  PF2500   9.986E+01 1.048E+02 1.611E+02 -2.979E+02 5.528E+01
end_S13
#======================================#
start_S14
  minE   -228.504086
  Vmin      0.00
  Vave   2484.56
  Vmax   4376.86
  V(0,0) 1374.30
  const  2512.15
  SP  MIN01  180.00    0.00      0.00  +157.00  +3910.00  21312.35  8.563  0.727  -0.033  21037.39
  SP  MIN02  167.47  129.97   1041.63  +163.00  +3853.00  21208.13  8.613  0.778  -0.345  20971.08
  SP  MIN03   29.59    0.67   1294.24  +104.00  +3903.00  21205.44  8.350  0.749  +0.374  20938.96
  SP  MIN04   40.39  210.51   2070.96  +130.00  +3885.00  21123.57  8.761  0.744  -0.147  20917.23
  SP  TS_01  175.57   86.48   1203.16  -261.00  +3847.00  21059.10  8.617  0.784  -0.215  20951.14
  SP  TS_02    0.00    0.00   1374.35   -80.50  +3908.00  21160.95  7.740  0.748  +0.393  20949.59
  SP  TS_03  180.00  180.00   1527.10  -298.00  +3957.00  21125.10  8.022  0.720  -0.336  21034.93
  SP  TS_04   29.85  281.08   2645.33  -339.20  +3855.00  20922.59  8.372  0.771  +0.222  20850.87
  SP  TS_05    0.00  180.00   2668.15  -185.20  +3930.00  20995.15  7.749  0.723  -0.218  20948.08
  SP  TS_06  100.49    0.71   2835.17  -151.80  +3881.00  21121.46  9.957  0.750  +0.184  20897.58
  SP  TS_07  100.14  185.05   3113.91  -142.40  +3913.00  21036.19  9.997  0.748  -0.344  20902.92
  SP  MAX01   99.00   82.00   4348.67  -369.30  +3845.00  20827.36  9.963  0.773  -0.227  20812.12
  SP  MAX02  255.00   81.00   4402.22  -382.20  +3849.00  20824.80  9.811  0.779  +0.123  20809.47
  TOR2DNS   248.32  390.87  532.11  586.63  672.21  732.46  811.26  876.82
  zpe1WHO  60.935
  zpeMSHO  60.935
  zpe2DNS   0.710
  zpeEHR   60.149
  zpeE2DT  60.859
  PF100    1.673E+04 1.673E+04 1.769E+04 1.157E+00 9.325E-02 1.426E+03
  PF100    1.129E-129 1.129E-129 1.752E-129 3.249E-02 9.325E-02 5.029E-129
  PF100    1.521E+06 1.000E+00
  PF100    1.718E-123 1.718E-123 2.665E-123
  PF100    6.650E-01 8.640E-01 5.625E+01 -4.761E+00 1.012E+01
  PF100    6.650E-01 8.640E-01 5.625E+01 -4.761E+00 1.013E+01
  PF100    6.760E-01 8.750E-01 5.648E+01 -4.772E+00 1.041E+01
  PF700    5.083E+07 8.828E+07 1.194E+08 1.272E+01 7.741E+00 7.266E+07
  PF700    4.806E-12 8.346E-12 1.192E-11 7.636E+00 7.741E+00 1.209E-11
  PF700    1.972E+08 1.000E+00
  PF700    9.476E-04 1.646E-03 2.351E-03
  PF700    1.419E+01 1.558E+01 9.547E+01 -5.125E+01 3.688E+01
  PF700    1.563E+01 1.702E+01 9.863E+01 -5.202E+01 4.015E+01
  PF700    1.572E+01 1.711E+01 9.935E+01 -5.244E+01 3.972E+01
  PF2500   2.148E+15 1.285E+16 1.404E+16 2.613E+02 2.265E+02 1.217E+16
  PF2500   1.012E+10 6.054E+10 6.716E+10 2.265E+02 2.265E+02 6.716E+10
  PF2500   4.753E+09 1.000E+00
  PF2500   4.811E+19 2.878E+20 3.192E+20
  PF2500   1.012E+02 1.061E+02 1.569E+02 -2.861E+02 5.685E+01
  PF2500   1.046E+02 1.095E+02 1.618E+02 -2.950E+02 5.716E+01
  PF2500   1.025E+02 1.074E+02 1.611E+02 -2.954E+02 5.546E+01
end_S14
#======================================#
start_S15
  minE   -228.513667
  Vmin      0.00
  Vave   1847.73
  Vmax   3644.39
  V(0,0)    0.00
  const  1899.06
  SP  MIN01  359.58    0.00      0.00  +159.00  +3365.00  21385.94  5.898  3.034  -0.694  21217.26
  SP  MIN02  180.00    0.00    738.31  +119.00  +3365.00  21272.36  7.290  2.964  +0.596  21132.29
  SP  TS_01  359.85   59.64    339.09   -20.80  +3362.00  21270.29  6.342  3.014  -0.700  21176.65
  SP  TS_02  180.00   60.40   1007.83  -102.00  +3372.00  21197.36  7.123  2.978  +0.613  21117.76
  SP  TS_03  265.72    1.21   3230.45  -175.80  +3341.00  21042.33  8.740  3.025  +0.042  20963.38
  SP  MAX01   93.00   60.00   3652.28  -180.00  +3351.00  20957.57  8.804  3.024  +0.044  20949.00
  TOR2DNS   155.97  155.98  155.98  290.17  290.17  290.68  319.18  319.18
  zpe1WHO  61.145
  zpeMSHO  61.145
  zpe2DNS   0.446
  zpeEHR   60.663
  zpeE2DT  61.109
  PF100    1.948E+04 1.948E+04 2.397E+04 1.320E+00 1.998E-01 3.628E+03
  PF100    4.561E-130 4.562E-130 6.739E-130 1.400E-01 1.998E-01 9.620E-130
  PF100    1.521E+06 1.000E+00
  PF100    6.938E-124 6.939E-124 1.025E-123
  PF100    7.120E-01 9.110E-01 5.702E+01 -4.792E+00 1.137E+01
  PF100    7.120E-01 9.110E-01 5.921E+01 -5.010E+00 1.139E+01
  PF100    7.360E-01 9.340E-01 5.767E+01 -4.833E+00 1.221E+01
  PF700    8.567E+07 1.249E+08 1.802E+08 2.067E+01 1.506E+01 1.313E+08
  PF700    6.962E-12 1.015E-11 1.503E-11 1.500E+01 1.506E+01 1.510E-11
  PF700    1.972E+08 1.000E+00
  PF700    1.373E-03 2.001E-03 2.964E-03
  PF700    1.384E+01 1.523E+01 9.600E+01 -5.198E+01 3.568E+01
  PF700    1.446E+01 1.585E+01 9.983E+01 -5.403E+01 3.661E+01
  PF700    1.437E+01 1.576E+01 9.825E+01 -5.301E+01 3.639E+01
  PF2500   2.633E+15 6.096E+15 7.333E+15 2.189E+02 2.001E+02 6.703E+15
  PF2500   1.189E+10 2.753E+10 3.336E+10 2.001E+02 2.001E+02 3.336E+10
  PF2500   4.753E+09 1.000E+00
  PF2500   5.653E+19 1.309E+20 1.586E+20
  PF2500   1.005E+02 1.055E+02 1.570E+02 -2.871E+02 5.696E+01
  PF2500   1.017E+02 1.066E+02 1.613E+02 -2.967E+02 5.705E+01
  PF2500   1.002E+02 1.051E+02 1.589E+02 -2.922E+02 5.562E+01
end_S15
#======================================#
start_S16
  minE   -155.243239
  Vmin      0.00
  Vave   1030.08
  Vmax   1939.45
  V(0,0) 1743.23
  const  1013.71
  SP  MIN01  119.31   59.67      0.00  +105.00  +3385.00  25608.12  7.755  2.831  -0.385  25428.50
  SP  MIN02    0.00   60.11    286.41  +151.00  +3395.00  25652.89  5.330  2.953  -1.278  25432.10
  SP  TS_01  180.00   59.98    748.19  -122.40  +3382.00  25513.31  7.570  2.675  -0.126  25380.35
  SP  TS_02   50.90   58.63    809.42  -133.30  +3386.00  25573.25  6.552  2.945  -1.034  25438.67
  SP  TS_03  119.36  358.63   1146.32  -248.60  +3385.00  25540.32  7.784  2.806  -0.398  25481.36
  SP  TS_04    0.00    0.00   1745.26  -274.10  +3400.00  25566.60  5.072  2.912  -1.259  25529.43
  SP  MAX01  180.00    0.00   1923.91  -255.30  +3382.00  25449.53  7.623  2.640  -0.172  25448.01
  SP  MAX02  305.00    1.00   1946.74  -253.60  +3383.00  25486.12  6.700  2.923  -1.013  25480.96
  TOR2DNS   167.44  167.44  167.44  167.44  167.44  167.44  270.66  270.66
  zpe1WHO  73.217
  zpeMSHO  73.217
  zpe2DNS   0.479
  zpeEHR   72.704
  zpeE2DT  73.182
  PF100    1.142E+04 2.292E+04 2.342E+04 2.704E+00 4.103E-01 3.553E+03
  PF100    1.108E-156 2.225E-156 2.709E-156 2.431E-01 4.103E-01 4.572E-156
  PF100    1.089E+06 1.000E+00
  PF100    1.207E-150 2.424E-150 2.950E-150
  PF100    7.100E-01 9.080E-01 5.527E+01 -4.619E+00 1.079E+01
  PF100    7.130E-01 9.120E-01 5.888E+01 -4.976E+00 1.095E+01
  PF100    7.220E-01 9.210E-01 5.683E+01 -4.762E+00 1.118E+01
  PF700    3.042E+07 7.134E+07 9.323E+07 3.873E+01 2.769E+01 6.667E+07
  PF700    4.209E-16 9.871E-16 1.323E-15 2.745E+01 2.769E+01 1.334E-15
  PF700    1.412E+08 1.000E+00
  PF700    5.943E-08 1.394E-07 1.868E-07
  PF700    1.357E+01 1.496E+01 9.291E+01 -5.007E+01 3.655E+01
  PF700    1.370E+01 1.509E+01 9.696E+01 -5.279E+01 3.664E+01
  PF700    1.378E+01 1.517E+01 9.543E+01 -5.163E+01 3.599E+01
  PF2500   1.614E+15 4.087E+15 3.954E+15 3.046E+02 2.766E+02 3.591E+15
  PF2500   6.416E+08 1.625E+09 1.583E+09 2.766E+02 2.766E+02 1.583E+09
  PF2500   3.404E+09 1.000E+00
  PF2500   2.184E+18 5.531E+18 5.389E+18
  PF2500   1.064E+02 1.114E+02 1.578E+02 -2.830E+02 6.181E+01
  PF2500   1.066E+02 1.115E+02 1.619E+02 -2.931E+02 6.182E+01
  PF2500   1.039E+02 1.089E+02 1.585E+02 -2.875E+02 5.998E+01
end_S16
#======================================#
start_S17
  minE   -155.243948
  Vmin      0.00
  Vave    224.49
  Vmax    414.82
  V(0,0)    0.00
  const   230.35
  SP  MIN01    0.00    0.00      0.00  +111.00  +3327.00  25516.47  2.918  2.918  -0.051  25395.56
  SP  TS_01   60.30    0.00    228.69   -91.00  +3338.00  25422.06  2.925  2.931  -0.034  25366.35
  SP  MAX01   60.00   60.00    414.81   -82.60  +3348.00  25341.39  2.935  2.935  -0.021  25340.86
  TOR2DNS   109.52  109.66  109.66  109.66  109.66  109.81  109.81  109.81
  zpe1WHO  72.955
  zpeMSHO  72.955
  zpe2DNS   0.313
  zpeEHR   72.610
  zpeE2DT  72.923
  PF100    6.848E+03 6.848E+03 1.364E+04 1.930E+00 4.520E-01 3.195E+03
  PF100    2.485E-156 2.485E-156 5.833E-156 3.992E-01 4.520E-01 6.604E-156
  PF100    1.089E+06 1.000E+00
  PF100    2.707E-150 2.707E-150 6.353E-150
  PF100    7.560E-01 9.550E-01 5.472E+01 -4.517E+00 1.165E+01
  PF100    7.560E-01 9.550E-01 5.909E+01 -4.954E+00 1.165E+01
  PF100    8.760E-01 1.074E+00 5.729E+01 -4.654E+00 1.205E+01
  PF700    3.180E+07 3.180E+07 6.046E+07 2.349E+01 1.878E+01 4.835E+07
  PF700    5.312E-16 5.312E-16 1.034E-15 1.876E+01 1.878E+01 1.035E-15
  PF700    1.412E+08 1.000E+00
  PF700    7.501E-08 7.501E-08 1.460E-07
  PF700    1.384E+01 1.523E+01 9.337E+01 -5.013E+01 3.653E+01
  PF700    1.384E+01 1.523E+01 9.774E+01 -5.319E+01 3.653E+01
  PF700    1.304E+01 1.443E+01 9.350E+01 -5.103E+01 3.465E+01
  PF2500   1.919E+15 1.919E+15 1.413E+15 9.838E+01 9.237E+01 1.327E+15
  PF2500   8.044E+08 8.044E+08 5.963E+08 9.237E+01 9.237E+01 5.963E+08
  PF2500   3.404E+09 1.000E+00
  PF2500   2.738E+18 2.738E+18 2.030E+18
  PF2500   1.066E+02 1.116E+02 1.582E+02 -2.839E+02 6.181E+01
  PF2500   1.066E+02 1.116E+02 1.626E+02 -2.948E+02 6.181E+01
  PF2500   1.023E+02 1.073E+02 1.559E+02 -2.823E+02 5.984E+01
end_S17
#======================================#
start_S18
  minE   -266.151289
  Vmin      0.00
  Vave   2975.30
  Vmax   5644.94
  V(0,0) 1173.65
  const  3023.08
  SP  MIN01  179.99    0.00      0.00  +141.00  +3403.00  22697.82  7.619  4.822  -2.615  22497.08
  SP  MIN02  180.00  180.00    221.23  +106.00  +3403.00  22593.00  8.729  8.169  -4.505  22418.21
  SP  MIN03   30.67    0.11   1042.94   +85.00  +3404.00  22650.93  8.504  7.086  -5.074  22476.81
  SP  MIN04   32.75  179.42   1213.48  +128.00  +3405.00  22563.43  5.648  6.805  -2.211  22398.65
  SP  TS_01    0.00    0.00   1173.75   -86.70  +3405.00  22633.42  8.471  7.046  -5.221  22516.17
  SP  TS_02    0.00  180.00   1355.69  -110.90  +3407.00  22535.41  5.056  6.646  -1.806  22441.13
  SP  TS_03  100.56  359.63   2244.57  -161.10  +3398.00  22490.87  8.222  6.257  -3.622  22397.26
  SP  TS_04  100.57  180.05   2383.93  -162.00  +3398.00  22410.58  8.129  7.599  -3.724  22319.58
  SP  TS_05  180.15   92.58   3756.53  -227.30  +3401.00  22347.82  8.403  7.208  -3.780  22259.14
  SP  TS_06   38.98   92.08   4624.55  -222.10  +3401.00  22327.67  6.662  6.545  -2.835  22245.21
  SP  TS_07   39.22  266.74   4662.96  -217.00  +3400.00  22317.82  8.096  7.547  -4.160  22243.14
  SP  MAX01    1.00  267.00   4929.18  -248.40  +3402.00  22290.80  7.066  7.027  -3.381  22272.42
  SP  MAX02  258.00   93.00   5636.77  -239.80  +3395.00  22199.04  9.136  7.876  -4.497  22187.88
  SP  MAX03  103.00   93.00   5659.59  -240.70  +3396.00  22200.15  7.388  6.240  -2.679  22188.46
  TOR2DNS   188.76  320.07  385.30  433.19  450.29  498.12  562.07  579.40
  zpe1WHO  64.896
  zpeMSHO  64.896
  zpe2DNS   0.540
  zpeEHR   64.322
  zpeE2DT  64.862
  PF100    3.902E+04 4.556E+04 4.697E+04 1.292E+00 1.551E-01 5.637E+03
  PF100    5.801E-138 6.773E-138 8.295E-138 8.549E-02 1.551E-01 1.505E-137
  PF100    1.928E+06 1.000E+00
  PF100    1.119E-131 1.306E-131 1.599E-131
  PF100    7.500E-01 9.490E-01 5.926E+01 -4.977E+00 1.238E+01
  PF100    8.010E-01 9.990E-01 6.007E+01 -5.008E+00 1.317E+01
  PF100    8.030E-01 1.002E+00 6.016E+01 -5.014E+00 1.334E+01
  PF700    6.133E+08 1.547E+09 1.837E+09 2.747E+01 1.879E+01 1.257E+09
  PF700    3.362E-12 8.480E-12 1.032E-11 1.863E+01 1.879E+01 1.041E-11
  PF700    2.500E+08 1.000E+00
  PF700    8.404E-04 2.120E-03 2.580E-03
  PF700    1.574E+01 1.713E+01 1.031E+02 -5.505E+01 4.034E+01
  PF700    1.671E+01 1.811E+01 1.063E+02 -5.633E+01 4.210E+01
  PF700    1.695E+01 1.835E+01 1.070E+02 -5.657E+01 4.259E+01
  PF2500   2.156E+17 1.183E+18 1.626E+18 5.411E+02 4.854E+02 1.459E+18
  PF2500   4.577E+11 2.511E+12 3.475E+12 4.854E+02 4.854E+02 3.475E+12
  PF2500   6.025E+09 1.000E+00
  PF2500   2.758E+21 1.513E+22 2.094E+22
  PF2500   1.122E+02 1.172E+02 1.709E+02 -3.102E+02 6.272E+01
  PF2500   1.143E+02 1.193E+02 1.752E+02 -3.186E+02 6.290E+01
  PF2500   1.141E+02 1.190E+02 1.757E+02 -3.202E+02 6.179E+01
end_S18
#======================================#
start_S19
  minE   -342.672656
  Vmin      0.00
  Vave   1249.82
  Vmax   2314.67
  V(0,0) 2314.67
  const  1248.89
  SP  MIN01   25.42   65.45      0.00   +90.00  +3859.00  31565.15  16.614  0.811  -0.680  31316.05
  SP  MIN02    0.00  180.01    713.07   +50.00  +3891.00  31410.83  15.866  0.776  +0.442  31284.20
  SP  TS_01    2.66  222.55    821.05  -243.60  +3903.00  31333.15  16.001  0.779  +0.231  31280.49
  SP  TS_02  130.43   66.81   1584.61  -104.20  +3883.00  31400.54  16.744  0.794  -0.611  31273.91
  SP  TS_03   88.20    0.00   1321.68  -275.80  +3890.00  31351.32  17.740  0.796  -1.113  31332.41
  SP  TS_04   88.63  180.00   1871.02   -71.90  +3849.00  31440.13  17.551  0.794  +0.536  31291.78
  SP  MAX01   90.00  128.00   2224.59  -338.90  +3887.00  31297.17  17.623  0.790  +0.202  31277.50
  SP  MAX02    0.00    0.00   2331.70  -522.10  +3960.00  31299.27  15.323  0.751  -1.072  31286.83
  TOR2DNS   241.78  241.78  241.79  241.79  328.84  328.84  328.87  328.87
  zpe1WHO  90.249
  zpeMSHO  90.249
  zpe2DNS   0.691
  zpeEHR   89.537
  zpeE2DT  90.228
  PF100    8.512E+04 1.703E+05 1.759E+05 2.829E+00 2.696E-01 1.676E+04
  PF100    4.949E-193 9.901E-193 1.136E-192 8.728E-02 2.696E-01 3.510E-192
  PF100    2.915E+06 1.000E+00
  PF100    1.442E-186 2.886E-186 3.312E-186
  PF100    7.490E-01 9.480E-01 6.162E+01 -5.214E+00 1.184E+01
  PF100    7.500E-01 9.490E-01 6.438E+01 -5.489E+00 1.187E+01
  PF100    7.540E-01 9.520E-01 6.310E+01 -5.358E+00 1.204E+01
  PF700    4.731E+09 1.378E+10 1.807E+10 3.638E+01 2.241E+01 1.113E+10
  PF700    3.151E-19 9.177E-19 1.222E-18 2.214E+01 2.241E+01 1.237E-18
  PF700    3.778E+08 1.000E+00
  PF700    1.191E-10 3.468E-10 4.616E-10
  PF700    1.977E+01 2.116E+01 1.137E+02 -5.846E+01 5.407E+01
  PF700    2.037E+01 2.176E+01 1.181E+02 -6.091E+01 5.491E+01
  PF700    2.045E+01 2.184E+01 1.174E+02 -6.033E+01 5.414E+01
  PF2500   1.789E+21 8.045E+21 7.569E+21 3.385E+02 2.945E+02 6.586E+21
  PF2500   2.307E+13 1.038E+14 9.806E+13 2.945E+02 2.945E+02 9.806E+13
  PF2500   9.108E+09 1.000E+00
  PF2500   2.102E+23 9.452E+23 8.931E+23
  PF2500   1.514E+02 1.563E+02 2.054E+02 -3.570E+02 8.479E+01
  PF2500   1.525E+02 1.575E+02 2.102E+02 -3.680E+02 8.487E+01
  PF2500   1.497E+02 1.547E+02 2.076E+02 -3.642E+02 8.300E+01
end_S19
#======================================#
start_S20
  minE   -342.683138
  Vmin      0.00
  Vave    682.61
  Vmax   1040.96
  V(0,0)  993.79
  const   677.77
  SP  MIN01  179.99   59.94      0.00  +129.00  +3919.00  31290.59  0.749  3.050  +0.046  31068.33
  SP  MIN02    0.05   60.16    447.73  +135.00  +3941.00  31235.27  0.752  3.041  -0.029  31051.76
  SP  TS_01  180.00    0.00    283.34  -107.50  +3920.00  31219.83  0.752  3.053  +0.046  31075.96
  SP  TS_02   92.65   49.24    773.65  -190.80  +3835.00  31043.88  0.787  3.049  +0.011  30981.48
  SP  TS_03    0.00    0.00    993.78  -173.10  +3964.00  31141.80  0.742  3.023  -0.030  31077.99
  SP  MAX01  309.00   15.00   1046.89  -171.50  +3882.00  31037.47  0.772  3.036  -0.016  31014.21
  TOR2DNS   203.30  203.36  203.36  308.38  308.38  310.20  382.45  395.56
  zpe1WHO  89.464
  zpeMSHO  89.464
  zpe2DNS   0.581
  zpeEHR   88.829
  zpeE2DT  89.410
  PF100    6.947E+04 6.972E+04 8.012E+04 1.356E+00 1.425E-01 8.419E+03
  PF100    2.098E-191 2.106E-191 3.179E-191 7.275E-02 1.425E-01 6.225E-191
  PF100    2.915E+06 1.000E+00
  PF100    6.115E-185 6.138E-185 9.266E-185
  PF100    7.310E-01 9.290E-01 6.103E+01 -5.173E+00 1.207E+01
  PF100    7.350E-01 9.330E-01 6.326E+01 -5.393E+00 1.231E+01
  PF100    7.780E-01 9.760E-01 6.178E+01 -5.202E+00 1.320E+01
  PF700    6.059E+09 9.463E+09 1.523E+10 2.028E+01 1.339E+01 1.006E+10
  PF700    7.095E-19 1.108E-18 1.854E-18 1.335E+01 1.339E+01 1.860E-18
  PF700    3.778E+08 1.000E+00
  PF700    2.681E-10 4.187E-10 7.005E-10
  PF700    2.030E+01 2.169E+01 1.150E+02 -5.881E+01 5.446E+01
  PF700    2.074E+01 2.213E+01 1.187E+02 -6.095E+01 5.484E+01
  PF700    2.040E+01 2.179E+01 1.170E+02 -6.009E+01 5.316E+01
  PF2500   3.159E+21 6.550E+21 5.280E+21 1.270E+02 1.130E+02 4.697E+21
  PF2500   4.773E+13 9.896E+13 8.065E+13 1.130E+02 1.130E+02 8.065E+13
  PF2500   9.108E+09 1.000E+00
  PF2500   4.347E+23 9.013E+23 7.345E+23
  PF2500   1.521E+02 1.571E+02 2.068E+02 -3.599E+02 8.481E+01
  PF2500   1.528E+02 1.577E+02 2.107E+02 -3.689E+02 8.485E+01
  PF2500   1.490E+02 1.540E+02 2.065E+02 -3.624E+02 8.289E+01
end_S20
#======================================#
"""

DATA_MPWB1K = """
#======================================#
start_S01
  minE   -302.955893
  Vmin      0.00
  Vave   3521.54
  Vmax   6336.91
  V(0,0) 3334.22
  const  3542.96
  SP  MIN01  180.00    0.00      0.00  +168.00  +3691.00   9574.60  6.954  0.786  -0.354   9161.09
  SP  MIN02  180.00  180.00    339.09  +100.00  +3771.00   9508.80  6.358  0.768  -0.077   9160.72
  SP  MIN03    0.00  180.00    787.69   +92.00  +3770.00   9481.37  6.549  0.783  +0.377   9135.03
  SP  MIN04    0.00    0.00   3334.48   +80.00  +3807.00   9361.91  6.494  0.719  -0.179   9093.94
  SP  TS_01   82.74  180.08   1586.14   -98.20  +3746.00   9367.68  7.840  0.784  +0.216   9081.40
  SP  TS_02   66.20    6.28   3760.70  -102.10  +3803.00   9286.16  7.611  0.720  -0.244   9063.17
  SP  TS_03  163.90   92.33   4449.19  -616.30  +3836.00   9097.39  6.610  0.741  -0.122   9019.86
  SP  TS_04   14.64   84.86   5554.90  -551.00  +3826.00   9031.04  6.690  0.734  +0.100   8983.92
  SP  MAX01   73.08   87.55   5716.88  -531.20  +3818.00   8960.27  7.815  0.746  +0.148   8937.82
  SP  MAX02  274.14   85.17   6369.81  -569.00  +3826.00   8937.36  7.919  0.741  -0.256   8912.75
  TOR2DNS   366.14  519.38  671.02  679.36  773.81  821.04  865.70  932.94
  zpe1WHO  27.375
  zpeMSHO  27.375
  zpe2DNS   1.047
  zpeEHR   26.193
  zpeE2DT  27.240
  PF100    1.426E+04 1.461E+04 1.520E+04 1.140E+00 5.363E-02 7.155E+02
  PF100    2.124E-56 2.176E-56 4.475E-56 5.874E-03 5.363E-02 4.086E-55
  PF100    1.652E+06 1.000E+00
  PF100    3.508E-50 3.594E-50 7.392E-50
  PF100    6.580E-01 8.570E-01 5.603E+01 -4.746E+00 9.903E+00
  PF100    6.780E-01 8.770E-01 5.628E+01 -4.751E+00 1.073E+01
  PF100    6.770E-01 8.760E-01 5.635E+01 -4.759E+00 1.066E+01
  PF700    1.288E+07 3.314E+07 4.182E+07 1.265E+01 6.290E+00 2.080E+07
  PF700    3.657E-02 9.411E-02 1.309E-01 5.960E+00 6.290E+00 1.381E-01
  PF700    2.141E+08 1.000E+00
  PF700    7.831E+06 2.015E+07 2.803E+07
  PF700    1.115E+01 1.254E+01 8.857E+01 -4.946E+01 2.666E+01
  PF700    1.200E+01 1.339E+01 9.166E+01 -5.077E+01 2.751E+01
  PF700    1.237E+01 1.377E+01 9.266E+01 -5.109E+01 2.807E+01
  PF2500   1.915E+12 9.522E+12 1.516E+13 1.830E+02 1.482E+02 1.228E+13
  PF2500   7.747E+09 3.852E+10 6.302E+10 1.482E+02 1.482E+02 6.302E+10
  PF2500   5.162E+09 1.000E+00
  PF2500   3.999E+19 1.988E+20 3.253E+20
  PF2500   6.709E+01 7.206E+01 1.295E+02 -2.516E+02 3.590E+01
  PF2500   6.930E+01 7.427E+01 1.335E+02 -2.596E+02 3.653E+01
  PF2500   6.922E+01 7.419E+01 1.344E+02 -2.619E+02 3.539E+01
end_S01
#======================================#
start_S02
  minE   -228.835127
  Vmin      0.00
  Vave   1772.08
  Vmax   3816.86
  V(0,0) 2149.47
  const  1768.97
  SP  MIN01    1.07  241.96      0.00  +176.00  +3789.00  13115.85  3.823  0.843  -0.296  12900.33
  SP  MIN02  199.24  109.05    535.08   +75.00  +3776.00  13015.69  7.047  0.806  +0.054  12844.56
  SP  MIN03  160.14  128.61    556.37   +65.00  +3793.00  13001.16  7.012  0.824  +0.081  12851.28
  SP  TS_01    0.00  180.00    114.57  -118.50  +3815.00  13054.30  3.768  0.826  -0.188  12928.44
  SP  TS_02  166.44  188.11    769.48  -179.40  +3813.00  12889.01  6.901  0.823  +0.249  12860.32
  SP  TS_03    0.00    0.00   2151.29  -504.20  +3775.00  12995.18  3.672  0.793  -0.586  12881.47
  SP  TS_04   72.67  115.58   2275.51  -185.90  +3780.00  12884.86  6.525  0.838  -0.283  12781.84
  SP  TS_05   75.77  241.91   2396.66  -193.70  +3794.00  12893.47  6.310  0.845  +0.041  12799.72
  SP  TS_06  175.78  121.68    567.34   -53.30  +3788.00  12970.98  6.950  0.816  +0.075  12846.86
  SP  TS_07  141.18    7.35   3012.51  -401.20  +3746.00  12878.68  7.072  0.764  -0.510  12817.48
  SP  MAX01  180.00  180.00    774.96  -175.60  +3813.00  12864.80  6.861  0.822  +0.249  12862.63
  SP  MAX02   76.81  197.00   2460.97  -187.20  +3802.00  12815.74  6.447  0.846  +0.102  12806.88
  SP  MAX03  180.00    0.00   3327.24  -429.70  +3752.00  12810.30  6.801  0.739  -0.500  12805.50
  SP  MAX04   80.32    6.73   3823.47  -358.60  +3736.00  12800.28  6.816  0.808  -0.603  12791.71
  TOR2DNS   184.64  201.33  299.07  382.16  415.63  431.97  485.68  531.13
  zpe1WHO  37.500
  zpeMSHO  37.500
  zpe2DNS   0.528
  zpeEHR   36.884
  zpeE2DT  37.412
  PF100    7.330E+03 1.473E+04 1.312E+04 2.135E+00 2.658E-01 1.634E+03
  PF100    8.137E-79 1.635E-78 2.272E-78 1.498E-01 2.658E-01 4.030E-78
  PF100    1.207E+06 1.000E+00
  PF100    9.817E-73 1.972E-72 2.741E-72
  PF100    6.690E-01 8.680E-01 5.419E+01 -4.551E+00 1.033E+01
  PF100    6.750E-01 8.740E-01 5.563E+01 -4.690E+00 1.073E+01
  PF100    7.100E-01 9.090E-01 5.576E+01 -4.667E+00 1.081E+01
  PF700    9.394E+06 5.095E+07 3.654E+07 3.037E+01 2.096E+01 2.522E+07
  PF700    1.841E-05 9.985E-05 7.631E-05 2.078E+01 2.096E+01 7.697E-05
  PF700    1.564E+08 1.000E+00
  PF700    2.880E+03 1.562E+04 1.194E+04
  PF700    1.176E+01 1.315E+01 8.819E+01 -4.858E+01 2.869E+01
  PF700    1.269E+01 1.408E+01 9.288E+01 -5.093E+01 2.925E+01
  PF700    1.235E+01 1.374E+01 9.172E+01 -5.047E+01 2.889E+01
  PF2500   4.563E+12 4.292E+13 2.104E+13 3.256E+02 2.928E+02 1.892E+13
  PF2500   2.405E+09 2.262E+10 1.129E+10 2.928E+02 2.928E+02 1.129E+10
  PF2500   3.770E+09 1.000E+00
  PF2500   9.066E+18 8.529E+19 4.255E+19
  PF2500   7.463E+01 7.960E+01 1.336E+02 -2.544E+02 4.088E+01
  PF2500   7.583E+01 8.080E+01 1.385E+02 -2.655E+02 4.092E+01
  PF2500   7.363E+01 7.859E+01 1.362E+02 -2.620E+02 3.938E+01
end_S02
#======================================#
start_S03
  minE   -228.946781
  Vmin      0.00
  Vave   2168.20
  Vmax   3509.54
  V(0,0)    0.00
  const  2187.36
  SP  MIN01    0.00    0.00      0.00  +198.00  +3733.00  13365.64  4.383  0.784  -0.511  13063.37
  SP  MIN02  180.00  180.00   1188.67   +86.00  +3872.00  13225.76  5.902  0.743  +0.285  13067.73
  SP  MIN03  193.84   79.16   1295.12   +81.00  +3858.00  13254.97  5.953  0.727  -0.042  13070.18
  SP  MIN04    0.00  180.01   2022.90  +166.00  +3866.00  13252.14  3.844  0.741  -0.128  13045.30
  SP  TS_01  178.19  121.22   1450.95  -221.90  +3894.00  13119.19  5.950  0.736  +0.133  13065.21
  SP  TS_02  281.09   62.81   2112.66  -138.20  +3828.00  13354.18  6.045  0.761  -0.140  13158.66
  SP  TS_03  180.00    0.00   2169.51  -323.30  +3832.00  13074.39  5.817  0.694  -0.306  13033.02
  SP  TS_04    6.00  225.33   2177.85  -261.20  +3893.00  13160.23  3.859  0.745  -0.166  13056.21
  SP  TS_05   66.60  193.45   3168.56  -163.50  +3862.00  13212.03  5.526  0.754  +0.087  13129.16
  SP  TS_06   67.10   68.35   3172.73  -173.10  +3844.00  13261.98  5.592  0.763  -0.462  13119.02
  SP  MAX01   68.64  212.26   3180.19  -182.80  +3878.00  13148.08  5.557  0.756  +0.112  13138.78
  SP  MAX02   60.78  125.67   3514.67  -260.90  +3887.00  13130.65  5.493  0.750  -0.211  13112.52
  SP  MAX03   78.17   15.31   3518.18  -356.40  +3852.00  13152.34  5.844  0.726  -0.480  13144.25
  TOR2DNS   289.40  487.80  644.76  683.34  841.55  874.64  955.84 1037.04
  zpe1WHO  38.214
  zpeMSHO  38.214
  zpe2DNS   0.827
  zpeEHR   37.350
  zpeE2DT  38.177
  PF100    7.187E+03 7.187E+03 7.142E+03 1.068E+00 6.065E-02 4.057E+02
  PF100    2.193E-80 2.193E-80 2.623E-80 1.660E-02 6.065E-02 9.581E-80
  PF100    1.207E+06 1.000E+00
  PF100    2.646E-74 2.646E-74 3.164E-74
  PF100    6.470E-01 8.460E-01 5.393E+01 -4.547E+00 9.763E+00
  PF100    6.470E-01 8.460E-01 5.393E+01 -4.547E+00 9.763E+00
  PF100    6.510E-01 8.490E-01 5.395E+01 -4.546E+00 9.876E+00
  PF700    5.060E+06 9.626E+06 1.126E+07 1.135E+01 6.342E+00 6.291E+06
  PF700    5.935E-06 1.129E-05 1.356E-05 6.262E+00 6.342E+00 1.373E-05
  PF700    1.564E+08 1.000E+00
  PF700    9.283E+02 1.766E+03 2.121E+03
  PF700    1.108E+01 1.247E+01 8.599E+01 -4.772E+01 2.779E+01
  PF700    1.285E+01 1.424E+01 8.980E+01 -4.861E+01 3.156E+01
  PF700    1.284E+01 1.423E+01 9.009E+01 -4.883E+01 3.052E+01
  PF2500   1.566E+12 1.188E+13 9.483E+12 2.355E+02 1.994E+02 8.028E+12
  PF2500   7.147E+08 5.424E+09 4.360E+09 1.994E+02 1.994E+02 4.360E+09
  PF2500   3.770E+09 1.000E+00
  PF2500   2.695E+18 2.045E+19 1.644E+19
  PF2500   7.370E+01 7.867E+01 1.311E+02 -2.490E+02 4.092E+01
  PF2500   7.711E+01 8.208E+01 1.365E+02 -2.591E+02 4.108E+01
  PF2500   7.431E+01 7.927E+01 1.349E+02 -2.580E+02 3.928E+01
end_S03
#======================================#
start_S04
  minE   -154.965066
  Vmin      0.00
  Vave    877.59
  Vmax   1896.16
  V(0,0) 1896.16
  const   851.47
  SP  MIN01  180.00   59.78      0.00  +241.00  +3861.00  17394.58  0.736  2.591  +0.027  17130.37
  SP  MIN02   60.68   56.80     41.70  +262.00  +3846.00  17402.98  0.732  2.594  -0.158  17127.49
  SP  TS_01  118.68   56.26    421.61  -281.40  +3895.00  17272.85  0.732  2.593  -0.040  17132.64
  SP  TS_02    0.00   60.38    502.16  -318.20  +3866.00  17281.00  0.706  2.597  -0.214  17138.70
  SP  TS_03  180.00    0.00   1152.24  -250.30  +3856.00  17325.54  0.738  2.577  +0.027  17178.57
  SP  TS_04  294.39    3.04   1299.51  -259.40  +3849.00  17321.95  0.733  2.571  -0.150  17173.91
  SP  MAX01  241.32    3.50   1624.55  -288.60  +3894.00  17190.41  0.733  2.575  -0.041  17178.24
  SP  MAX02    0.00    0.00   1896.04  -337.50  +3878.00  17200.90  0.702  2.561  -0.215  17195.71
  TOR2DNS   251.50  251.50  251.50  304.11  304.11  304.11  306.27  306.27
  zpe1WHO  49.734
  zpeMSHO  49.734
  zpe2DNS   0.719
  zpeEHR   48.978
  zpeE2DT  49.697
  PF100    3.330E+03 6.619E+03 6.730E+03 2.069E+00 1.449E-01 4.713E+02
  PF100    6.792E-106 1.350E-105 1.648E-105 5.549E-02 1.449E-01 4.303E-105
  PF100    8.106E+05 1.000E+00
  PF100    5.506E-100 1.094E-99 1.336E-99
  PF100    6.350E-01 8.340E-01 5.150E+01 -4.315E+00 9.516E+00
  PF100    7.040E-01 9.030E-01 5.573E+01 -4.670E+00 9.696E+00
  PF100    7.180E-01 9.170E-01 5.373E+01 -4.455E+00 1.007E+01
  PF700    1.859E+06 5.100E+06 5.866E+06 1.826E+01 1.100E+01 3.533E+06
  PF700    5.523E-10 1.515E-09 1.788E-09 1.089E+01 1.100E+01 1.806E-09
  PF700    1.051E+08 1.000E+00
  PF700    5.804E-02 1.592E-01 1.879E-01
  PF700    1.105E+01 1.244E+01 8.317E+01 -4.577E+01 2.884E+01
  PF700    1.114E+01 1.253E+01 8.748E+01 -4.871E+01 2.885E+01
  PF700    1.100E+01 1.239E+01 8.537E+01 -4.737E+01 2.797E+01
  PF2500   1.163E+12 3.338E+12 2.457E+12 1.083E+02 9.373E+01 2.126E+12
  PF2500   5.226E+07 1.499E+08 1.112E+08 9.373E+01 9.373E+01 1.112E+08
  PF2500   2.533E+09 1.000E+00
  PF2500   1.324E+17 3.798E+17 2.816E+17
  PF2500   8.002E+01 8.499E+01 1.322E+02 -2.456E+02 4.582E+01
  PF2500   8.011E+01 8.507E+01 1.365E+02 -2.563E+02 4.582E+01
  PF2500   7.711E+01 8.208E+01 1.326E+02 -2.493E+02 4.397E+01
end_S04
#======================================#
start_S05
  minE   -267.054074
  Vmin      0.00
  Vave   3764.06
  Vmax   7195.66
  V(0,0)    0.00
  const  3865.07
  SP  MIN01    0.00    0.00      0.00  +118.00  +3779.00  14726.61  8.511  0.797  +0.427  14364.91
  SP  MIN02  180.00    0.00     87.13  +113.00  +3787.00  14722.00  8.208  0.779  -0.044  14378.47
  SP  MIN03    0.00  180.00   2278.59   +93.00  +3835.00  14630.88  8.534  0.729  -0.257  14355.23
  SP  MIN04  157.24  171.45   2614.38   +94.00  +3830.00  14634.54  8.235  0.746  -0.394  14345.99
  SP  TS_01   90.73    0.40   2418.17  -141.80  +3758.00  14575.39  9.506  0.795  +0.221  14270.58
  SP  TS_02  180.00  180.00   2674.52   -88.40  +3840.00  14549.06  7.743  0.733  -0.364  14346.84
  SP  TS_03   87.62  182.64   3840.37  -121.80  +3817.00  14535.66  9.625  0.750  -0.355  14273.42
  SP  TS_04  170.68   99.06   4488.48  -517.50  +3840.00  14338.60  8.312  0.767  -0.268  14248.17
  SP  TS_05    4.44  263.12   4748.11  -563.10  +3839.00  14298.50  8.610  0.750  +0.068  14231.59
  SP  MAX01  261.34   93.28   7125.68  -600.40  +3840.00  14133.48  9.540  0.750  +0.103  14110.16
  SP  MAX02   90.63   94.54   7223.57  -600.90  +3839.00  14131.27  9.600  0.752  -0.278  14108.11
  TOR2DNS   357.53  420.97  474.21  529.63  590.20  638.48  705.44  747.19
  zpe1WHO  42.105
  zpeMSHO  42.105
  zpe2DNS   1.022
  zpeEHR   41.071
  zpeE2DT  42.094
  PF100    1.703E+04 2.227E+04 2.351E+04 1.738E+00 9.028E-02 1.221E+03
  PF100    1.628E-88 2.129E-88 2.386E-88 1.014E-02 9.028E-02 2.124E-87
  PF100    1.586E+06 1.000E+00
  PF100    2.581E-82 3.377E-82 3.784E-82
  PF100    6.910E-01 8.900E-01 5.663E+01 -4.773E+00 1.040E+01
  PF100    7.470E-01 9.460E-01 5.772E+01 -4.826E+00 1.091E+01
  PF100    7.450E-01 9.430E-01 5.781E+01 -4.837E+00 1.090E+01
  PF700    4.089E+07 7.778E+07 8.857E+07 1.380E+01 6.995E+00 4.488E+07
  PF700    2.924E-06 5.562E-06 6.388E-06 6.620E+00 6.995E+00 6.750E-06
  PF700    2.056E+08 1.000E+00
  PF700    6.011E+02 1.143E+03 1.313E+03
  PF700    1.311E+01 1.450E+01 9.358E+01 -5.101E+01 3.266E+01
  PF700    1.334E+01 1.473E+01 9.519E+01 -5.190E+01 3.348E+01
  PF700    1.359E+01 1.498E+01 9.580E+01 -5.208E+01 3.412E+01
  PF2500   1.547E+14 4.962E+14 6.643E+14 1.768E+02 1.439E+02 5.408E+14
  PF2500   3.227E+10 1.035E+11 1.389E+11 1.439E+02 1.439E+02 1.389E+11
  PF2500   4.956E+09 1.000E+00
  PF2500   1.599E+20 5.129E+20 6.883E+20
  PF2500   8.547E+01 9.044E+01 1.455E+02 -2.732E+02 4.662E+01
  PF2500   8.823E+01 9.320E+01 1.489E+02 -2.790E+02 4.755E+01
  PF2500   8.796E+01 9.293E+01 1.494E+02 -2.805E+02 4.641E+01
end_S05
#======================================#
start_S06
  minE   -267.038611
  Vmin      0.00
  Vave   5815.41
  Vmax   8717.21
  V(0,0)    0.00
  const  5908.15
  SP  MIN01    0.00    0.00      0.00  +265.00  +3273.00  14817.89  6.575  0.853  -0.467  14245.93
  SP  MIN02  180.00    0.00   3640.65  +152.00  +3823.00  14635.16  6.890  0.719  -0.001  14340.85
  SP  MIN03  180.00  180.00   3694.64  +152.00  +3873.00  14582.85  6.998  0.739  +0.283  14346.16
  SP  MIN04    0.00  180.00   4789.16  +142.00  +3860.00  14623.99  5.638  0.735  +0.085  14350.34
  SP  TS_01  183.10   83.75   5652.13  -457.30  +3836.00  14323.75  7.000  0.762  +0.148  14211.37
  SP  TS_02   17.95  255.19   6615.62  -521.00  +3833.00  14271.18  6.177  0.756  -0.006  14183.84
  SP  TS_03  263.05    9.99   6755.43  -209.20  +3819.00  14397.32  8.620  0.743  -0.177  14191.75
  SP  TS_04   87.43  182.53   7061.82  -179.10  +3881.00  14318.70  8.677  0.736  +0.261  14182.09
  SP  MAX01  268.96   89.98   8014.56  -364.90  +3829.00  14124.56  8.884  0.776  +0.250  14091.29
  SP  MAX02   87.65   82.09   8724.78  -425.20  +3839.00  14111.03  8.701  0.769  -0.274  14065.44
  TOR2DNS   447.98  655.56  862.31 1068.14 1124.24 1272.99 1329.59 1476.75
  zpe1WHO  42.366
  zpeMSHO  42.366
  zpe2DNS   1.281
  zpeEHR   40.731
  zpeE2DT  42.012
  PF100    1.315E+04 1.315E+04 1.383E+04 1.053E+00 3.325E-02 4.366E+02
  PF100    3.382E-89 3.382E-89 2.116E-88 1.672E-03 3.325E-02 4.208E-87
  PF100    1.586E+06 1.000E+00
  PF100    5.363E-83 5.363E-83 3.356E-82
  PF100    6.340E-01 8.330E-01 5.555E+01 -4.722E+00 9.566E+00
  PF100    6.340E-01 8.330E-01 5.555E+01 -4.722E+00 9.566E+00
  PF100    6.470E-01 8.460E-01 5.578E+01 -4.732E+00 9.866E+00
  PF700    1.482E+07 1.493E+07 2.063E+07 3.909E+00 1.691E+00 8.926E+06
  PF700    8.787E-07 8.852E-07 1.578E-06 1.556E+00 1.691E+00 1.714E-06
  PF700    2.056E+08 1.000E+00
  PF700    1.807E+02 1.820E+02 3.244E+02
  PF700    1.246E+01 1.385E+01 9.063E+01 -4.959E+01 3.203E+01
  PF700    1.254E+01 1.393E+01 9.076E+01 -4.960E+01 3.292E+01
  PF700    1.292E+01 1.431E+01 9.195E+01 -5.005E+01 3.372E+01
  PF2500   3.825E+13 1.297E+14 2.701E+14 5.887E+01 4.549E+01 2.088E+14
  PF2500   7.570E+09 2.566E+10 5.741E+10 4.549E+01 4.549E+01 5.741E+10
  PF2500   4.956E+09 1.000E+00
  PF2500   3.751E+19 1.272E+20 2.845E+20
  PF2500   8.481E+01 8.978E+01 1.424E+02 -2.663E+02 4.673E+01
  PF2500   9.287E+01 9.784E+01 1.481E+02 -2.723E+02 4.895E+01
  PF2500   9.385E+01 9.882E+01 1.499E+02 -2.760E+02 4.782E+01
end_S06
#======================================#
start_S07
  minE   -292.250951
  Vmin      0.00
  Vave   1619.38
  Vmax   3354.70
  V(0,0)    0.00
  const  1658.92
  SP  MIN01    0.00    0.00      0.00  +108.00  +3840.00  16928.68  0.761  10.871  +0.490  16685.10
  SP  MIN02  199.36  141.74     49.82  +117.00  +3805.00  16993.21  0.777  13.122  -0.447  16720.07
  SP  MIN03    5.81  242.47    439.39   +90.00  +3837.00  16974.60  0.753  13.176  +0.213  16733.49
  SP  MIN04  203.89    6.29   1272.73  +112.00  +3831.00  16844.12  0.740  10.824  -0.245  16674.17
  SP  TS_01  180.00  180.00    584.24  -119.40  +3799.00  16924.45  0.775  12.063  -0.518  16697.02
  SP  TS_02    1.28   62.79   1044.92   -99.20  +3830.00  16923.19  0.762  12.787  +0.381  16711.47
  SP  TS_03  180.00    0.00   1308.29  -165.00  +3844.00  16744.06  0.730  10.781  -0.292  16670.65
  SP  TS_04   79.63  238.95   1585.92  -383.20  +3853.00  16748.76  0.777  13.323  +0.096  16675.52
  SP  TS_05  264.23    3.03   1646.06  -304.10  +3831.00  16711.21  0.765  10.951  +0.083  16635.38
  SP  TS_06    0.00  180.00   1721.78  -136.20  +3832.00  16886.00  0.748  11.445  +0.003  16670.64
  SP  TS_07  201.26   58.52   2082.59  -112.30  +3845.00  16821.97  0.742  12.610  -0.305  16694.02
  SP  TS_08   86.41  118.23   2356.50  -388.00  +3839.00  16729.63  0.766  13.194  -0.313  16668.96
  SP  MAX01  274.49   59.13   2789.30  -357.20  +3832.00  16659.22  0.768  12.790  +0.196  16646.15
  SP  MAX02   94.79   61.83   3038.41  -344.90  +3836.00  16666.09  0.767  12.904  -0.207  16653.40
  SP  MAX03   80.59  174.76   3386.05  -427.90  +3848.00  16634.87  0.768  11.498  -0.216  16616.55
  TOR2DNS   231.85  311.53  311.53  336.40  418.54  418.54  438.92  522.30
  zpe1WHO  48.402
  zpeMSHO  48.402
  zpe2DNS   0.663
  zpeEHR   47.705
  zpeE2DT  48.368
  PF100    2.119E+04 2.888E+04 2.899E+04 2.118E+00 2.293E-01 3.139E+03
  PF100    3.524E-102 4.801E-102 5.707E-102 7.537E-02 2.293E-01 1.736E-101
  PF100    1.720E+06 1.000E+00
  PF100    6.062E-96 8.260E-96 9.817E-96
  PF100    7.070E-01 9.060E-01 5.739E+01 -4.833E+00 1.089E+01
  PF100    7.930E-01 9.920E-01 5.886E+01 -4.894E+00 1.196E+01
  PF100    8.160E-01 1.014E+00 5.909E+01 -4.895E+00 1.212E+01
  PF700    1.205E+08 4.056E+08 5.232E+08 4.193E+01 2.650E+01 3.307E+08
  PF700    9.328E-08 3.139E-07 4.148E-07 2.604E+01 2.650E+01 4.222E-07
  PF700    2.230E+08 1.000E+00
  PF700    2.080E+01 7.001E+01 9.251E+01
  PF700    1.445E+01 1.584E+01 9.780E+01 -5.262E+01 3.603E+01
  PF700    1.514E+01 1.653E+01 1.012E+02 -5.431E+01 3.698E+01
  PF700    1.527E+01 1.666E+01 1.019E+02 -5.467E+01 3.656E+01
  PF2500   3.002E+15 1.713E+16 1.735E+16 4.772E+02 4.176E+02 1.518E+16
  PF2500   1.763E+11 1.006E+12 1.026E+12 4.176E+02 4.176E+02 1.026E+12
  PF2500   5.376E+09 1.000E+00
  PF2500   9.477E+20 5.409E+21 5.515E+21
  PF2500   9.542E+01 1.004E+02 1.555E+02 -2.884E+02 5.209E+01
  PF2500   9.687E+01 1.018E+02 1.595E+02 -2.970E+02 5.226E+01
  PF2500   9.444E+01 9.941E+01 1.586E+02 -2.971E+02 5.046E+01
end_S07
#======================================#
start_S08
  minE   -193.020341
  Vmin      0.00
  Vave   1668.52
  Vmax   2900.32
  V(0,0) 1251.59
  const  1688.65
  SP  MIN01    0.00   60.66      0.00  +240.00  +3240.00  18562.67  4.216  2.979  -1.500  18310.38
  SP  MIN02  171.36   57.35    957.57   +31.00  +3243.00  18338.92  6.972  2.494  -0.513  18240.33
  SP  TS_01  180.00   61.07    958.45   -19.80  +3243.00  18324.58  6.967  2.485  -0.510  18241.06
  SP  TS_02    0.00    0.00   1251.88  -237.50  +3241.00  18438.87  3.882  2.938  -1.498  18346.14
  SP  TS_03  180.00    0.00   1323.21  -152.90  +3244.00  18297.35  6.862  2.451  -0.626  18249.40
  SP  TS_04   76.80   61.85   2318.75  -160.40  +3229.00  18296.09  6.389  2.918  -1.036  18205.62
  SP  MAX01   84.45    1.85   2918.57  -191.90  +3229.00  18218.27  6.474  2.887  -1.007  18211.52
  TOR2DNS   242.67  242.67  242.67  472.90  472.90  472.90  485.15  485.15
  zpe1WHO  53.073
  zpeMSHO  53.073
  zpe2DNS   0.694
  zpeEHR   52.352
  zpeE2DT  53.046
  PF100    7.320E+03 7.321E+03 7.460E+03 1.071E+00 7.916E-02 5.513E+02
  PF100    7.504E-113 7.505E-113 8.782E-113 3.262E-02 7.916E-02 2.131E-112
  PF100    1.147E+06 1.000E+00
  PF100    8.609E-107 8.611E-107 1.008E-106
  PF100    6.460E-01 8.450E-01 5.386E+01 -4.541E+00 9.917E+00
  PF100    6.470E-01 8.450E-01 5.605E+01 -4.759E+00 9.951E+00
  PF100    6.530E-01 8.520E-01 5.397E+01 -4.545E+00 1.012E+01
  PF700    1.011E+07 4.945E+07 2.384E+07 1.362E+01 8.328E+00 1.458E+07
  PF700    2.721E-10 1.331E-09 6.547E-10 8.271E+00 8.328E+00 6.593E-10
  PF700    1.487E+08 1.000E+00
  PF700    4.048E-02 1.980E-01 9.738E-02
  PF700    1.275E+01 1.414E+01 8.965E+01 -4.861E+01 3.333E+01
  PF700    1.480E+01 1.619E+01 9.791E+01 -5.235E+01 3.453E+01
  PF700    1.399E+01 1.538E+01 9.312E+01 -4.981E+01 3.427E+01
  PF2500   6.545E+13 1.047E+15 2.565E+14 1.950E+02 1.696E+02 2.231E+14
  PF2500   1.501E+09 2.401E+10 5.914E+09 1.696E+02 1.696E+02 5.914E+09
  PF2500   3.585E+09 1.000E+00
  PF2500   5.381E+18 8.610E+19 2.121E+19
  PF2500   9.170E+01 9.667E+01 1.456E+02 -2.673E+02 5.171E+01
  PF2500   9.420E+01 9.917E+01 1.543E+02 -2.866E+02 5.177E+01
  PF2500   9.108E+01 9.605E+01 1.481E+02 -2.741E+02 5.002E+01
end_S08
#======================================#
start_S09
  minE   -193.025820
  Vmin      0.00
  Vave    897.31
  Vmax   1881.18
  V(0,0)  806.31
  const   909.94
  SP  MIN01    6.00   60.47      0.00  +183.00  +3834.00  18523.55  5.005  0.780  -0.468  18276.64
  SP  MIN02  235.43   57.23      3.29  +118.00  +3840.00  18546.08  7.301  0.740  -0.189  18313.88
  SP  MIN03    0.00  180.00    202.14  +150.00  +3872.00  18483.79  4.902  0.743  -0.102  18285.25
  SP  MIN04  129.74  170.57    489.65  +111.00  +3858.00  18467.01  7.192  0.758  +0.227  18284.44
  SP  MIN05  122.56   64.74    576.34  +101.00  +3849.00  18461.57  7.205  0.760  -0.378  18295.39
  SP  TS_01    2.62  235.28    423.37  -245.90  +3897.00  18398.41  4.935  0.755  -0.201  18288.90
  SP  TS_02  125.87  114.72    760.04  -235.30  +3888.00  18365.79  7.278  0.754  -0.077  18296.83
  SP  TS_03  118.47   15.78    767.50  -275.40  +3858.00  18366.92  7.207  0.736  -0.483  18307.18
  SP  TS_04    0.00    0.00    817.54  -400.50  +3877.00  18357.91  4.730  0.738  -0.528  18276.40
  SP  TS_05  227.06  130.85    855.51  -301.40  +3891.00  18348.07  7.182  0.749  +0.217  18282.85
  SP  TS_06  180.00  180.00    950.54  -116.30  +3865.00  18363.40  6.980  0.752  +0.300  18233.47
  SP  TS_07  170.26   69.15    998.61  -128.00  +3852.00  18385.36  7.044  0.743  -0.202  18242.83
  SP  TS_08   67.54   49.05   1106.15  -134.40  +3842.00  18444.10  6.693  0.773  -0.565  18316.79
  SP  TS_09  306.43   45.61   1365.79  -174.60  +3855.00  18376.38  6.118  0.753  -0.351  18291.78
  SP  TS_10   63.56  181.41   1496.16  -161.50  +3844.00  18434.36  6.497  0.760  +0.066  18290.58
  SP  MAX01  178.09  122.07   1289.63  -270.80  +3897.00  18249.23  7.045  0.744  +0.126  18235.08
  SP  MAX02  307.87   29.65   1372.37  -168.40  +3856.00  18305.29  6.146  0.751  -0.432  18294.94
  SP  MAX03  180.00    0.00   1651.77  -345.80  +3870.00  18253.42  6.966  0.707  -0.394  18246.90
  SP  MAX04   59.44  126.08   1897.14  -311.90  +3888.00  18303.38  6.468  0.750  -0.214  18288.74
  SP  MAX05   65.88  239.11   1908.55  -297.50  +3892.00  18317.20  6.372  0.763  +0.054  18300.96
  TOR2DNS   228.89  228.89  232.08  232.32  346.87  346.87  384.91  404.26
  zpe1WHO  52.961
  zpeMSHO  52.961
  zpe2DNS   0.654
  zpeEHR   52.256
  zpeE2DT  52.910
  PF100    8.550E+03 3.102E+04 2.856E+04 4.759E+00 4.161E-01 2.497E+03
  PF100    1.539E-112 5.582E-112 6.662E-112 1.767E-01 4.161E-01 1.569E-111
  PF100    1.147E+06 1.000E+00
  PF100    1.766E-106 6.405E-106 7.643E-106
  PF100    6.640E-01 8.630E-01 5.435E+01 -4.572E+00 1.029E+01
  PF100    7.220E-01 9.200E-01 5.748E+01 -4.828E+00 1.089E+01
  PF100    7.510E-01 9.500E-01 5.762E+01 -4.812E+00 1.127E+01
  PF700    1.655E+07 1.340E+08 1.551E+08 6.630E+01 4.165E+01 9.741E+07
  PF700    4.828E-10 3.909E-09 4.695E-09 4.142E+01 4.165E+01 4.721E-09
  PF700    1.487E+08 1.000E+00
  PF700    7.181E-02 5.814E-01 6.983E-01
  PF700    1.311E+01 1.450E+01 9.114E+01 -4.930E+01 3.382E+01
  PF700    1.372E+01 1.511E+01 9.617E+01 -5.221E+01 3.428E+01
  PF700    1.354E+01 1.493E+01 9.621E+01 -5.241E+01 3.302E+01
  PF2500   1.347E+14 1.597E+15 1.072E+15 4.881E+02 4.279E+02 9.397E+14
  PF2500   3.159E+09 3.745E+10 2.541E+10 4.279E+02 4.279E+02 2.541E+10
  PF2500   3.585E+09 1.000E+00
  PF2500   1.133E+19 1.343E+20 9.109E+19
  PF2500   9.212E+01 9.709E+01 1.472E+02 -2.709E+02 5.163E+01
  PF2500   9.298E+01 9.795E+01 1.525E+02 -2.832E+02 5.167E+01
  PF2500   8.964E+01 9.461E+01 1.503E+02 -2.812E+02 4.975E+01
end_S09
#======================================#
start_S10
  minE   -193.042897
  Vmin      0.00
  Vave   1490.75
  Vmax   2830.52
  V(0,0)    0.00
  const  1533.78
  SP  MIN01    0.00    0.00      0.00  +180.00  +3833.00  18411.31  0.743  2.937  +0.067  18105.72
  SP  MIN02  180.00    0.00    781.11  +170.00  +3860.00  18311.20  0.737  2.930  -0.103  18108.54
  SP  TS_01    0.00   60.04    724.71  -189.90  +3827.00  18333.15  0.745  2.925  +0.067  18106.82
  SP  TS_02  180.00   60.05   1668.89  -215.10  +3871.00  18237.78  0.735  2.908  -0.102  18120.76
  SP  TS_03   89.04  359.20   1935.55  -387.90  +3844.00  18161.63  0.761  2.937  -0.015  18046.52
  SP  MAX01   88.81   59.63   2831.44  -394.70  +3842.00  18060.60  0.762  2.923  -0.015  18048.29
  TOR2DNS   292.82  292.82  292.82  463.09  463.09  463.11  621.34  621.82
  zpe1WHO  52.641
  zpeMSHO  52.641
  zpe2DNS   0.837
  zpeEHR   51.767
  zpeE2DT  52.604
  PF100    8.422E+03 8.422E+03 8.619E+03 1.100E+00 6.570E-02 5.148E+02
  PF100    7.620E-112 7.621E-112 9.372E-112 1.628E-02 6.570E-02 3.782E-111
  PF100    1.147E+06 1.000E+00
  PF100    8.743E-106 8.744E-106 1.075E-105
  PF100    6.470E-01 8.460E-01 5.415E+01 -4.569E+00 9.694E+00
  PF100    6.470E-01 8.460E-01 5.633E+01 -4.787E+00 9.706E+00
  PF100    6.540E-01 8.520E-01 5.426E+01 -4.573E+00 9.897E+00
  PF700    1.551E+07 2.213E+07 2.518E+07 8.752E+00 4.898E+00 1.409E+07
  PF700    5.699E-10 8.132E-10 9.499E-10 4.794E+00 4.898E+00 9.705E-10
  PF700    1.487E+08 1.000E+00
  PF700    8.477E-02 1.210E-01 1.413E-01
  PF700    1.332E+01 1.471E+01 9.131E+01 -4.921E+01 3.415E+01
  PF700    1.398E+01 1.537E+01 9.514E+01 -5.123E+01 3.521E+01
  PF700    1.393E+01 1.532E+01 9.315E+01 -4.988E+01 3.468E+01
  PF2500   1.463E+14 3.426E+14 2.802E+14 8.453E+01 7.142E+01 2.368E+14
  PF2500   3.659E+09 8.571E+09 7.062E+09 7.142E+01 7.142E+01 7.062E+09
  PF2500   3.585E+09 1.000E+00
  PF2500   1.312E+19 3.073E+19 2.532E+19
  PF2500   9.245E+01 9.742E+01 1.475E+02 -2.713E+02 5.163E+01
  PF2500   9.373E+01 9.870E+01 1.519E+02 -2.810E+02 5.173E+01
  PF2500   9.119E+01 9.615E+01 1.483E+02 -2.746E+02 4.997E+01
end_S10
#======================================#
start_S11
  minE   -193.036365
  Vmin      0.00
  Vave    793.89
  Vmax   1686.27
  V(0,0)    0.00
  const   818.75
  SP  MIN01    0.00    0.00      0.00   +47.00  +3839.00  18463.09  0.738  2.880  -0.105  18254.40
  SP  MIN02  180.00    0.00    209.38   +99.00  +3890.00  18413.51  0.736  2.918  +0.084  18257.19
  SP  TS_01    0.00   60.62    107.32   -97.70  +3826.00  18422.71  0.745  2.899  -0.097  18235.35
  SP  TS_02  180.00   59.43    514.01  -141.90  +3886.00  18348.13  0.733  2.914  +0.085  18222.67
  SP  TS_03   83.00    7.26   1270.32  -350.30  +3841.00  18255.08  0.763  2.925  -0.020  18165.95
  SP  MAX01   85.43   65.74   1702.03  -385.30  +3832.00  18161.04  0.763  2.913  -0.015  18144.70
  TOR2DNS   211.54  212.72  212.72  264.26  264.26  276.44  302.25  330.02
  zpe1WHO  52.789
  zpeMSHO  52.789
  zpe2DNS   0.605
  zpeEHR   52.192
  zpeE2DT  52.797
  PF100    1.651E+04 1.763E+04 1.651E+04 1.870E+00 2.303E-01 2.033E+03
  PF100    7.090E-112 7.572E-112 6.804E-112 8.913E-02 2.303E-01 1.758E-111
  PF100    1.147E+06 1.000E+00
  PF100    8.135E-106 8.688E-106 7.807E-106
  PF100    7.600E-01 9.590E-01 5.662E+01 -4.703E+00 1.101E+01
  PF100    7.880E-01 9.870E-01 5.921E+01 -4.934E+00 1.162E+01
  PF100    7.580E-01 9.570E-01 5.660E+01 -4.703E+00 1.109E+01
  PF700    6.069E+07 9.596E+07 6.517E+07 1.926E+01 1.262E+01 4.271E+07
  PF700    2.005E-09 3.170E-09 2.141E-09 1.247E+01 1.262E+01 2.167E-09
  PF700    1.487E+08 1.000E+00
  PF700    2.982E-01 4.715E-01 3.184E-01
  PF700    1.335E+01 1.474E+01 9.407E+01 -5.111E+01 3.383E+01
  PF700    1.357E+01 1.496E+01 9.748E+01 -5.327E+01 3.394E+01
  PF700    1.328E+01 1.467E+01 9.410E+01 -5.120E+01 3.301E+01
  PF2500   5.566E+14 1.003E+15 3.971E+14 1.164E+02 1.031E+02 3.516E+14
  PF2500   1.352E+10 2.436E+10 9.628E+09 1.031E+02 1.031E+02 9.628E+09
  PF2500   3.585E+09 1.000E+00
  PF2500   4.846E+19 8.733E+19 3.452E+19
  PF2500   9.232E+01 9.729E+01 1.501E+02 -2.780E+02 5.162E+01
  PF2500   9.260E+01 9.757E+01 1.536E+02 -2.864E+02 5.162E+01
  PF2500   8.945E+01 9.442E+01 1.483E+02 -2.763E+02 4.978E+01
end_S11
#======================================#
start_S12
  minE   -231.121999
  Vmin      0.00
  Vave   1999.52
  Vmax   3973.25
  V(0,0) 1157.38
  const  2024.37
  SP  MIN01  180.00    0.00      0.00  +146.00  +3204.00  19472.51  6.705  3.039  -0.791  19305.21
  SP  MIN02    0.00    0.00   1157.73  +111.00  +3218.00  19396.95  6.747  2.960  +0.388  19268.14
  SP  TS_01  180.00   60.03    409.10  -138.80  +3202.00  19415.85  6.581  3.037  -0.783  19327.03
  SP  TS_02    0.00   60.19   1587.46  -142.30  +3216.00  19339.56  6.792  2.948  +0.381  19282.10
  SP  TS_03   92.60    1.86   3213.33  -157.10  +3203.00  19276.38  8.217  3.029  -0.244  19178.71
  SP  MAX01   93.34   61.22   3990.49  -197.60  +3200.00  19191.83  8.210  3.009  -0.249  19185.04
  TOR2DNS   144.24  144.25  144.25  268.79  268.79  269.15  301.84  301.84
  zpe1WHO  55.675
  zpeMSHO  55.675
  zpe2DNS   0.412
  zpeEHR   55.196
  zpeE2DT  55.609
  PF100    1.978E+04 1.978E+04 2.143E+04 1.353E+00 2.353E-01 3.728E+03
  PF100    4.186E-118 4.186E-118 6.321E-118 1.698E-01 2.353E-01 8.759E-118
  PF100    1.521E+06 1.000E+00
  PF100    6.367E-112 6.367E-112 9.614E-112
  PF100    7.230E-01 9.220E-01 5.717E+01 -4.795E+00 1.169E+01
  PF100    7.230E-01 9.220E-01 5.935E+01 -5.013E+00 1.169E+01
  PF100    7.510E-01 9.500E-01 5.761E+01 -4.811E+00 1.217E+01
  PF700    1.865E+08 2.185E+08 2.776E+08 1.831E+01 1.366E+01 2.071E+08
  PF700    7.736E-10 9.065E-10 1.208E-09 1.361E+01 1.366E+01 1.212E-09
  PF700    1.972E+08 1.000E+00
  PF700    1.525E-01 1.787E-01 2.382E-01
  PF700    1.510E+01 1.649E+01 9.935E+01 -5.306E+01 3.818E+01
  PF700    1.557E+01 1.696E+01 1.025E+02 -5.481E+01 3.952E+01
  PF700    1.545E+01 1.684E+01 1.007E+02 -5.361E+01 3.906E+01
  PF2500   1.820E+16 3.497E+16 3.419E+16 2.017E+02 1.857E+02 3.146E+16
  PF2500   2.472E+11 4.751E+11 4.706E+11 1.857E+02 1.857E+02 4.706E+11
  PF2500   4.753E+09 1.000E+00
  PF2500   1.175E+21 2.258E+21 2.237E+21
  PF2500   1.042E+02 1.091E+02 1.623E+02 -2.967E+02 5.759E+01
  PF2500   1.057E+02 1.107E+02 1.664E+02 -3.054E+02 5.781E+01
  PF2500   1.039E+02 1.089E+02 1.635E+02 -2.998E+02 5.631E+01
end_S12
#======================================#
start_S13
  minE   -231.104744
  Vmin      0.00
  Vave   1913.44
  Vmax   3911.20
  V(0,0) 1022.60
  const  1954.23
  SP  MIN01  180.00    0.00      0.00  +102.00  +3838.00  19682.37  8.848  0.731  -0.060  19450.71
  SP  MIN02  180.00  180.00    146.83  +130.00  +3884.00  19624.91  8.867  0.745  +0.342  19439.85
  SP  MIN03   38.97    5.93    232.42  +153.00  +3745.00  19780.44  8.696  0.788  -0.456  19435.92
  SP  MIN04    0.00  180.01    916.09  +107.00  +3876.00  19737.60  6.921  0.737  +0.119  19531.36
  SP  TS_01    0.00    0.00   1021.87  -227.30  +3860.00  19677.77  6.414  0.760  -0.442  19502.08
  SP  TS_02  183.29   84.48   1418.46  -380.10  +3832.00  19439.26  8.902  0.772  +0.147  19347.42
  SP  TS_03  256.25    5.96   1428.12  -118.40  +3783.00  19606.47  10.593  0.774  -0.211  19378.32
  SP  TS_04   16.67   88.53   2374.06  -421.90  +3811.00  19455.92  7.435  0.779  -0.256  19379.86
  SP  TS_05   92.63  179.93   2669.47  -126.40  +3883.00  19491.49  10.679  0.744  +0.307  19360.33
  SP  MAX01   90.52  266.37   3694.64  -394.10  +3840.00  19291.35  10.753  0.773  +0.262  19272.76
  SP  MAX02   97.24   89.42   3917.84  -409.30  +3838.00  19288.01  10.693  0.773  -0.233  19261.74
  TOR2DNS   236.63  338.32  346.26  454.03  467.63  540.30  540.30  560.17
  zpe1WHO  56.275
  zpeMSHO  56.275
  zpe2DNS   0.677
  zpeEHR   55.612
  zpeE2DT  56.289
  PF100    2.151E+04 2.744E+04 2.734E+04 1.581E+00 1.432E-01 2.476E+03
  PF100    2.223E-119 2.836E-119 2.631E-119 5.253E-02 1.432E-01 7.170E-119
  PF100    1.521E+06 1.000E+00
  PF100    3.381E-113 4.313E-113 4.001E-113
  PF100    7.350E-01 9.330E-01 5.745E+01 -4.811E+00 1.153E+01
  PF100    7.980E-01 9.970E-01 5.857E+01 -4.860E+00 1.258E+01
  PF100    7.990E-01 9.980E-01 5.857E+01 -4.859E+00 1.269E+01
  PF700    2.121E+08 5.606E+08 6.828E+08 3.011E+01 1.887E+01 4.278E+08
  PF700    5.717E-10 1.511E-09 1.822E-09 1.852E+01 1.887E+01 1.856E-09
  PF700    1.972E+08 1.000E+00
  PF700    1.127E-01 2.979E-01 3.592E-01
  PF700    1.531E+01 1.671E+01 9.992E+01 -5.324E+01 3.870E+01
  PF700    1.573E+01 1.712E+01 1.024E+02 -5.459E+01 3.904E+01
  PF700    1.607E+01 1.746E+01 1.033E+02 -5.487E+01 3.952E+01
  PF2500   2.368E+16 8.220E+16 1.106E+17 3.626E+02 3.164E+02 9.648E+16
  PF2500   2.850E+11 9.896E+11 1.327E+12 3.164E+02 3.164E+02 1.327E+12
  PF2500   4.753E+09 1.000E+00
  PF2500   1.355E+21 4.704E+21 6.308E+21
  PF2500   1.042E+02 1.092E+02 1.629E+02 -2.980E+02 5.740E+01
  PF2500   1.049E+02 1.099E+02 1.656E+02 -3.042E+02 5.746E+01
  PF2500   1.039E+02 1.089E+02 1.658E+02 -3.057E+02 5.598E+01
end_S13
#======================================#
start_S14
  minE   -231.107723
  Vmin      0.00
  Vave   2299.00
  Vmax   4411.46
  V(0,0) 1177.47
  const  2327.00
  SP  MIN01  180.00    0.00      0.00  +144.00  +3847.00  19602.92  8.546  0.744  -0.041  19341.54
  SP  MIN02  169.65  158.13    652.72  +141.00  +3833.00  19559.98  8.393  0.763  -0.362  19320.66
  SP  MIN03   29.56    2.38   1069.28  +102.00  +3838.00  19525.69  8.411  0.763  +0.369  19269.38
  SP  MIN04   41.35  198.06   1576.93  +123.00  +3848.00  19503.46  8.811  0.745  -0.200  19289.18
  SP  TS_01  180.00  180.00    733.05  -190.30  +3865.00  19435.06  8.144  0.745  -0.348  19330.60
  SP  TS_02    0.00    0.00   1177.04   -86.60  +3842.00  19465.90  7.762  0.762  +0.390  19261.81
  SP  TS_03  175.93   88.56   1516.79  -353.20  +3842.00  19379.96  8.577  0.781  -0.210  19280.60
  SP  TS_04    0.00  180.00   2033.21  -142.70  +3860.00  19351.68  7.693  0.736  -0.211  19256.42
  SP  TS_05  102.19    0.30   2576.41  -142.50  +3818.00  19455.68  9.934  0.763  +0.177  19237.73
  SP  TS_06  101.87  184.32   2599.68  -134.60  +3847.00  19423.71  10.012  0.759  -0.345  19259.74
  SP  TS_07   32.45  277.82   2682.64  -375.70  +3849.00  19293.29  8.493  0.770  +0.206  19220.87
  SP  TS_08   27.44   93.83   2871.83  -349.00  +3838.00  19279.71  8.381  0.771  -0.029  19221.92
  SP  MAX01    6.55   90.90   2918.13  -349.70  +3841.00  19232.73  7.848  0.771  +0.067  19219.91
  SP  MAX02  101.36   84.11   4425.05  -416.70  +3841.00  19196.20  9.941  0.770  -0.235  19178.87
  SP  MAX03  253.35   83.95   4442.82  -419.80  +3847.00  19200.74  9.805  0.774  +0.110  19183.20
  TOR2DNS   244.44  380.94  516.33  589.31  650.81  728.88  784.50  825.74
  zpe1WHO  56.048
  zpeMSHO  56.048
  zpe2DNS   0.699
  zpeEHR   55.300
  zpeE2DT  55.999
  PF100    1.716E+04 1.717E+04 1.779E+04 1.172E+00 9.943E-02 1.509E+03
  PF100    5.563E-119 5.565E-119 7.359E-119 3.480E-02 9.943E-02 2.103E-118
  PF100    1.521E+06 1.000E+00
  PF100    8.462E-113 8.464E-113 1.119E-112
  PF100    6.800E-01 8.790E-01 5.645E+01 -4.766E+00 1.055E+01
  PF100    6.810E-01 8.790E-01 5.646E+01 -4.766E+00 1.060E+01
  PF100    6.870E-01 8.860E-01 5.660E+01 -4.774E+00 1.077E+01
  PF700    1.174E+08 2.613E+08 3.017E+08 1.596E+01 9.803E+00 1.853E+08
  PF700    3.727E-10 8.291E-10 9.912E-10 9.655E+00 9.803E+00 1.006E-09
  PF700    1.972E+08 1.000E+00
  PF700    7.348E-02 1.635E-01 1.954E-01
  PF700    1.538E+01 1.677E+01 9.884E+01 -5.242E+01 3.904E+01
  PF700    1.676E+01 1.815E+01 1.024E+02 -5.353E+01 4.110E+01
  PF700    1.684E+01 1.823E+01 1.028E+02 -5.373E+01 4.105E+01
  PF2500   1.415E+16 8.638E+16 8.531E+16 2.914E+02 2.532E+02 7.412E+16
  PF2500   1.783E+11 1.088E+12 1.086E+12 2.532E+02 2.532E+02 1.086E+12
  PF2500   4.753E+09 1.000E+00
  PF2500   8.475E+20 5.174E+21 5.160E+21
  PF2500   1.045E+02 1.094E+02 1.619E+02 -2.955E+02 5.741E+01
  PF2500   1.070E+02 1.120E+02 1.666E+02 -3.044E+02 5.760E+01
  PF2500   1.052E+02 1.101E+02 1.658E+02 -3.044E+02 5.598E+01
end_S14
#======================================#
start_S15
  minE   -231.118206
  Vmin      0.00
  Vave   1338.29
  Vmax   3029.55
  V(0,0)  322.54
  const  1371.51
  SP  MIN01  180.00    0.00      0.00  +119.00  +3157.00  19534.54  7.282  2.971  +0.595  19393.82
  SP  MIN02    3.15    3.89    324.16   +19.00  +3157.00  19470.51  5.943  3.038  -0.708  19390.74
  SP  MIN03    0.00   58.19    325.26  +105.00  +3152.00  19529.11  6.358  3.028  -0.716  19375.47
  SP  TS_01    8.62   33.10    339.09   -40.40  +3149.00  19474.23  6.278  3.030  -0.702  19384.58
  SP  TS_02    0.00    0.00    324.16    -5.50  +3159.00  19462.36  5.927  3.038  -0.709  19392.47
  SP  TS_03  180.00   60.10    421.39  -131.00  +3163.00  19443.56  7.106  2.979  +0.609  19371.97
  SP  TS_04   92.98    1.31   2516.28  -174.20  +3139.00  19319.91  8.921  3.028  +0.035  19235.34
  SP  MAX01   92.18   62.25   3036.87  -184.40  +3146.00  19224.76  8.968  3.024  +0.032  19211.80
  TOR2DNS   134.57  134.58  134.58  253.15  253.15  253.32  278.72  278.72
  zpe1WHO  55.852
  zpeMSHO  55.852
  zpe2DNS   0.385
  zpeEHR   55.450
  zpeE2DT  55.834
  PF100    2.287E+04 2.677E+04 2.568E+04 1.473E+00 2.802E-01 4.886E+03
  PF100    1.983E-118 2.321E-118 2.433E-118 2.124E-01 2.802E-01 3.209E-118
  PF100    1.521E+06 1.000E+00
  PF100    3.016E-112 3.530E-112 3.700E-112
  PF100    7.570E-01 9.560E-01 5.779E+01 -4.823E+00 1.221E+01
  PF100    8.820E-01 1.080E+00 6.153E+01 -5.073E+00 1.688E+01
  PF100    8.090E-01 1.008E+00 5.854E+01 -4.846E+00 1.380E+01
  PF700    2.635E+08 2.163E+09 5.731E+08 3.382E+01 2.572E+01 4.358E+08
  PF700    9.624E-10 7.902E-09 2.120E-09 2.565E+01 2.572E+01 2.125E-09
  PF700    1.972E+08 1.000E+00
  PF700    1.898E-01 1.558E+00 4.180E-01
  PF700    1.505E+01 1.644E+01 9.998E+01 -5.354E+01 3.798E+01
  PF700    1.580E+01 1.720E+01 1.074E+02 -5.800E+01 3.807E+01
  PF700    1.523E+01 1.662E+01 1.018E+02 -5.462E+01 3.781E+01
  PF2500   2.449E+16 3.009E+17 4.916E+16 2.803E+02 2.594E+02 4.549E+16
  PF2500   3.210E+11 3.945E+12 6.466E+11 2.594E+02 2.594E+02 6.466E+11
  PF2500   4.753E+09 1.000E+00
  PF2500   1.526E+21 1.875E+22 3.073E+21
  PF2500   1.040E+02 1.090E+02 1.629E+02 -2.982E+02 5.758E+01
  PF2500   1.048E+02 1.098E+02 1.704E+02 -3.161E+02 5.759E+01
  PF2500   1.025E+02 1.075E+02 1.637E+02 -3.016E+02 5.605E+01
end_S15
#======================================#
start_S16
  minE   -157.131112
  Vmin      0.00
  Vave    972.42
  Vmax   1995.62
  V(0,0) 1392.82
  const   963.91
  SP  MIN01    0.00   60.09      0.00  +166.00  +3205.00  23554.28  5.365  2.954  -1.302  23336.19
  SP  MIN02  120.77   58.24     24.80  +106.00  +3199.00  23514.35  7.691  2.818  -0.408  23342.27
  SP  TS_01  180.00   59.87    741.39  -123.30  +3198.00  23419.80  7.500  2.666  -0.147  23295.66
  SP  TS_02   57.26   59.86    828.52  -136.40  +3199.00  23470.15  6.772  2.944  -0.999  23345.54
  SP  TS_03  237.74    2.02   1118.00  -237.20  +3199.00  23441.71  7.723  2.787  -0.400  23382.56
  SP  TS_04    0.00    0.00   1391.47  -260.50  +3208.00  23478.82  5.064  2.910  -1.276  23424.61
  SP  MAX01  180.00    0.00   1868.17  -245.30  +3199.00  23347.66  7.539  2.631  -0.187  23345.66
  SP  MAX02   61.88    0.54   1997.88  -247.30  +3197.00  23391.19  6.897  2.915  -0.967  23385.86
  TOR2DNS   193.69  193.69  193.69  193.69  193.69  193.69  208.64  208.64
  zpe1WHO  67.302
  zpeMSHO  67.302
  zpe2DNS   0.554
  zpeEHR   66.722
  zpeE2DT  67.275
  PF100    1.135E+04 3.062E+04 2.918E+04 3.601E+00 3.860E-01 3.128E+03
  PF100    9.339E-144 2.518E-143 2.742E-143 2.219E-01 3.860E-01 4.771E-143
  PF100    1.089E+06 1.000E+00
  PF100    1.017E-137 2.743E-137 2.987E-137
  PF100    7.190E-01 9.180E-01 5.536E+01 -4.618E+00 1.113E+01
  PF100    7.190E-01 9.180E-01 5.952E+01 -5.033E+00 1.099E+01
  PF100    7.370E-01 9.360E-01 5.741E+01 -4.805E+00 1.122E+01
  PF700    5.886E+07 1.532E+08 1.923E+08 4.342E+01 2.936E+01 1.300E+08
  PF700    5.724E-14 1.489E-13 1.906E-13 2.916E+01 2.936E+01 1.919E-13
  PF700    1.412E+08 1.000E+00
  PF700    8.082E-06 2.103E-05 2.691E-05
  PF700    1.475E+01 1.614E+01 9.590E+01 -5.099E+01 3.913E+01
  PF700    1.473E+01 1.612E+01 9.996E+01 -5.385E+01 3.913E+01
  PF700    1.482E+01 1.621E+01 9.836E+01 -5.264E+01 3.849E+01
  PF2500   9.926E+15 2.562E+16 2.386E+16 3.148E+02 2.816E+02 2.135E+16
  PF2500   1.298E+10 3.351E+10 3.138E+10 2.816E+02 2.816E+02 3.138E+10
  PF2500   3.404E+09 1.000E+00
  PF2500   4.419E+19 1.140E+20 1.068E+20
  PF2500   1.103E+02 1.152E+02 1.629E+02 -2.920E+02 6.255E+01
  PF2500   1.103E+02 1.152E+02 1.670E+02 -3.022E+02 6.255E+01
  PF2500   1.076E+02 1.126E+02 1.636E+02 -2.964E+02 6.071E+01
end_S16
#======================================#
start_S17
  minE   -157.134201
  Vmin      0.00
  Vave    295.08
  Vmax    560.24
  V(0,0)    0.00
  const   302.89
  SP  MIN01    0.00    0.00      0.00  +120.00  +3129.00  23419.85  2.920  2.920  -0.051  23292.37
  SP  TS_01   60.05    0.00    308.80  -113.20  +3139.00  23329.47  2.925  2.936  -0.033  23269.78
  SP  MAX01   59.48   59.48    560.32   -91.70  +3148.00  23249.44  2.938  2.938  -0.020  23248.87
  TOR2DNS   123.18  123.23  123.23  123.23  123.23  123.28  123.28  123.28
  zpe1WHO  66.961
  zpeMSHO  66.961
  zpe2DNS   0.352
  zpeEHR   66.596
  zpeE2DT  66.948
  PF100    6.566E+03 6.566E+03 9.609E+03 1.632E+00 3.354E-01 1.975E+03
  PF100    3.005E-143 3.005E-143 4.679E-143 2.773E-01 3.354E-01 5.659E-143
  PF100    1.089E+06 1.000E+00
  PF100    3.273E-137 3.273E-137 5.096E-137
  PF100    7.560E-01 9.550E-01 5.464E+01 -4.509E+00 1.188E+01
  PF100    7.560E-01 9.550E-01 5.901E+01 -4.946E+00 1.188E+01
  PF100    8.650E-01 1.064E+00 5.649E+01 -4.585E+00 1.334E+01
  PF700    5.511E+07 5.511E+07 1.016E+08 2.124E+01 1.652E+01 7.905E+07
  PF700    6.849E-14 6.849E-14 1.274E-13 1.649E+01 1.652E+01 1.277E-13
  PF700    1.412E+08 1.000E+00
  PF700    9.670E-06 9.670E-06 1.799E-05
  PF700    1.499E+01 1.638E+01 9.611E+01 -5.090E+01 3.915E+01
  PF700    1.499E+01 1.638E+01 1.005E+02 -5.396E+01 3.915E+01
  PF700    1.435E+01 1.574E+01 9.641E+01 -5.175E+01 3.732E+01
  PF2500   1.052E+16 1.052E+16 8.214E+15 9.531E+01 8.879E+01 7.652E+15
  PF2500   1.474E+10 1.474E+10 1.154E+10 8.879E+01 8.879E+01 1.154E+10
  PF2500   3.404E+09 1.000E+00
  PF2500   5.017E+19 5.017E+19 3.926E+19
  PF2500   1.105E+02 1.155E+02 1.631E+02 -2.923E+02 6.255E+01
  PF2500   1.105E+02 1.155E+02 1.675E+02 -3.032E+02 6.255E+01
  PF2500   1.064E+02 1.114E+02 1.610E+02 -2.911E+02 6.058E+01
end_S17
#======================================#
start_S18
  minE   -269.189840
  Vmin      0.00
  Vave   3004.28
  Vmax   5744.20
  V(0,0) 1814.61
  const  3036.01
  SP  MIN01  180.00  180.00      0.00   +93.00  +3216.00  20676.96  8.786  8.206  -4.498  20515.11
  SP  MIN02  180.00    0.00    533.76  +126.00  +3216.00  20688.51  7.729  4.817  -2.583  20522.20
  SP  MIN03   27.12  178.91   1088.16  +117.00  +3217.00  20660.95  5.543  6.784  -2.078  20512.88
  SP  MIN04  330.66    0.15   1699.39   +80.00  +3216.00  20653.86  8.535  7.114  -5.087  20509.05
  SP  TS_01    0.00  180.00   1172.21   -97.20  +3218.00  20631.18  5.124  6.678  -1.764  20547.63
  SP  TS_02    0.00    0.00   1814.18   -87.90  +3216.00  20625.67  8.485  7.052  -5.215  20544.99
  SP  TS_03   98.95  179.95   2575.53  -176.80  +3213.00  20499.83  7.967  7.635  -3.603  20416.64
  SP  TS_04   99.35  359.28   3191.60  -180.50  +3212.00  20493.59  8.222  6.259  -3.616  20425.97
  SP  TS_05  179.82   92.31   3445.97  -222.70  +3214.00  20447.77  8.571  7.310  -3.873  20353.48
  SP  TS_06   34.93   91.34   4468.50  -216.40  +3213.00  20433.94  6.781  6.714  -3.005  20351.58
  SP  TS_07   35.45  266.51   4510.86  -213.70  +3212.00  20425.80  8.121  7.629  -4.173  20349.23
  SP  MAX01    1.31  267.56   4726.17  -243.30  +3213.00  20395.66  7.178  7.149  -3.470  20370.38
  SP  MAX02  258.68   92.73   5734.43  -239.40  +3208.00  20306.73  9.203  7.983  -4.551  20290.73
  SP  MAX03  100.90   92.55   5754.19  -240.20  +3210.00  20309.47  7.369  6.376  -2.735  20292.22
  TOR2DNS   154.00  261.60  353.98  369.06  461.02  476.36  553.15  567.85
  zpe1WHO  59.118
  zpeMSHO  59.118
  zpe2DNS   0.440
  zpeEHR   58.656
  zpeE2DT  59.096
  PF100    3.924E+04 3.925E+04 3.798E+04 1.347E+00 2.236E-01 6.305E+03
  PF100    2.473E-125 2.474E-125 2.680E-125 1.469E-01 2.236E-01 4.079E-125
  PF100    1.928E+06 1.000E+00
  PF100    4.769E-119 4.771E-119 5.167E-119
  PF100    8.040E-01 1.003E+00 5.981E+01 -4.978E+00 1.316E+01
  PF100    8.050E-01 1.004E+00 5.982E+01 -4.978E+00 1.321E+01
  PF100    8.000E-01 9.990E-01 5.970E+01 -4.971E+00 1.331E+01
  PF700    1.831E+09 3.250E+09 3.457E+09 2.448E+01 1.797E+01 2.537E+09
  PF700    6.390E-10 1.134E-09 1.226E-09 1.784E+01 1.797E+01 1.235E-09
  PF700    2.500E+08 1.000E+00
  PF700    1.597E-01 2.836E-01 3.065E-01
  PF700    1.712E+01 1.851E+01 1.073E+02 -5.657E+01 4.285E+01
  PF700    1.824E+01 1.963E+01 1.100E+02 -5.737E+01 4.510E+01
  PF700    1.846E+01 1.986E+01 1.104E+02 -5.745E+01 4.563E+01
  PF2500   2.199E+18 1.002E+19 1.204E+19 5.251E+02 4.806E+02 1.102E+19
  PF2500   1.494E+13 6.805E+13 8.214E+13 4.806E+02 4.806E+02 8.214E+13
  PF2500   6.025E+09 1.000E+00
  PF2500   9.001E+22 4.100E+23 4.949E+23
  PF2500   1.161E+02 1.211E+02 1.771E+02 -3.217E+02 6.338E+01
  PF2500   1.187E+02 1.237E+02 1.812E+02 -3.292E+02 6.363E+01
  PF2500   1.183E+02 1.232E+02 1.813E+02 -3.301E+02 6.239E+01
end_S18
#======================================#
start_S19
  minE   -346.621591
  Vmin      0.00
  Vave    627.38
  Vmax   1396.20
  V(0,0) 1396.20
  const   622.53
  SP  MIN01   35.65   58.67      0.00   +54.00  +3836.00  29028.46  16.091  0.799  -0.749  28845.02
  SP  MIN02    0.00  180.00    327.46    +3.00  +3871.00  28929.33  15.105  0.773  +0.400  28819.65
  SP  TS_01   88.73    0.00    536.18  -261.20  +3857.00  28869.63  17.198  0.788  -1.087  28842.93
  SP  TS_02    4.74  236.69    560.32  -236.80  +3896.00  28871.85  15.247  0.777  +0.057  28825.71
  SP  TS_03   89.21  180.00    703.86   -29.60  +3844.00  28970.40  16.783  0.784  +0.500  28834.32
  SP  TS_04  133.48   70.09    782.21   -78.00  +3861.00  28937.02  16.028  0.788  -0.564  28832.31
  SP  MAX01   78.08  127.16   1065.11  -304.90  +3885.00  28846.32  16.886  0.778  +0.161  28829.82
  SP  MAX02    0.00    0.00   1404.86  -440.60  +3892.00  28815.24  14.895  0.763  -1.056  28804.00
  TOR2DNS   183.89  183.89  183.89  183.89  236.06  236.06  236.08  236.08
  zpe1WHO  82.997
  zpeMSHO  82.997
  zpe2DNS   0.526
  zpeEHR   82.472
  zpeE2DT  82.998
  PF100    1.218E+05 2.901E+05 2.665E+05 4.073E+00 5.819E-01 3.808E+04
  PF100    5.022E-177 1.196E-176 1.091E-176 2.890E-01 5.819E-01 2.197E-176
  PF100    2.915E+06 1.000E+00
  PF100    1.464E-170 3.484E-170 3.181E-170
  PF100    8.170E-01 1.016E+00 6.301E+01 -5.285E+00 1.291E+01
  PF100    9.330E-01 1.132E+00 6.727E+01 -5.595E+00 1.659E+01
  PF100    8.490E-01 1.048E+00 6.489E+01 -5.441E+00 1.398E+01
  PF700    3.182E+10 4.071E+11 1.270E+11 7.440E+01 5.123E+01 8.747E+10
  PF700    3.896E-16 4.984E-15 1.554E-15 5.098E+01 5.123E+01 1.561E-15
  PF700    3.778E+08 1.000E+00
  PF700    1.472E-07 1.883E-06 5.871E-07
  PF700    2.169E+01 2.308E+01 1.203E+02 -6.111E+01 5.733E+01
  PF700    2.241E+01 2.380E+01 1.278E+02 -6.562E+01 5.749E+01
  PF700    2.168E+01 2.307E+01 1.230E+02 -6.304E+01 5.599E+01
  PF2500   6.164E+22 1.172E+24 1.536E+23 4.441E+02 3.995E+02 1.381E+23
  PF2500   3.423E+15 6.508E+16 8.527E+15 3.995E+02 3.995E+02 8.527E+15
  PF2500   9.108E+09 1.000E+00
  PF2500   3.118E+25 5.927E+26 7.766E+25
  PF2500   1.564E+02 1.613E+02 2.144E+02 -3.746E+02 8.558E+01
  PF2500   1.572E+02 1.621E+02 2.219E+02 -3.927E+02 8.560E+01
  PF2500   1.531E+02 1.581E+02 2.149E+02 -3.792E+02 8.365E+01
end_S19
#======================================#
start_S20
  minE   -346.633849
  Vmin      0.00
  Vave    853.29
  Vmax   1533.71
  V(0,0)  780.05
  const   859.02
  SP  MIN01  180.00   59.26      0.00   +91.00  +3849.00  28747.73  0.765  3.057  +0.047  28534.97
  SP  MIN02    0.00   60.64    158.90  +158.00  +3848.00  28761.44  0.774  3.043  -0.028  28543.03
  SP  TS_01  180.00    0.00    340.84  -127.70  +3852.00  28729.05  0.767  3.053  +0.046  28573.61
  SP  TS_02    0.00    0.00    778.04  -177.30  +3865.00  28699.48  0.770  3.016  -0.029  28569.26
  SP  TS_03   92.79   50.87   1278.44  -316.80  +3834.00  28550.21  0.783  3.055  +0.010  28483.46
  SP  MAX01  274.37   10.15   1539.61  -286.60  +3834.00  28518.04  0.785  3.051  +0.005  28505.32
  TOR2DNS   225.61  225.64  225.64  336.63  336.63  337.46  378.55  378.55
  zpe1WHO  82.194
  zpeMSHO  82.194
  zpe2DNS   0.645
  zpeEHR   81.586
  zpeE2DT  82.231
  PF100    8.373E+04 8.950E+04 8.266E+04 1.414E+00 1.332E-01 7.791E+03
  PF100    1.959E-175 2.094E-175 1.608E-175 5.503E-02 1.332E-01 3.893E-175
  PF100    2.915E+06 1.000E+00
  PF100    5.710E-169 6.104E-169 4.686E-169
  PF100    7.830E-01 9.820E-01 6.193E+01 -5.211E+00 1.300E+01
  PF100    8.130E-01 1.011E+00 6.454E+01 -5.442E+00 1.361E+01
  PF100    8.080E-01 1.007E+00 6.215E+01 -5.208E+00 1.380E+01
  PF700    3.020E+10 4.469E+10 4.496E+10 1.801E+01 1.143E+01 2.854E+10
  PF700    6.583E-16 9.742E-16 9.547E-16 1.133E+01 1.143E+01 9.635E-16
  PF700    3.778E+08 1.000E+00
  PF700    2.487E-07 3.681E-07 3.607E-07
  PF700    2.218E+01 2.357E+01 1.209E+02 -6.104E+01 5.777E+01
  PF700    2.232E+01 2.371E+01 1.240E+02 -6.311E+01 5.781E+01
  PF700    2.227E+01 2.366E+01 1.218E+02 -6.159E+01 5.697E+01
  PF2500   7.984E+22 1.281E+23 8.493E+22 1.186E+02 1.042E+02 7.459E+22
  PF2500   5.212E+15 8.360E+15 5.503E+15 1.042E+02 1.042E+02 5.503E+15
  PF2500   9.108E+09 1.000E+00
  PF2500   4.747E+25 7.614E+25 5.012E+25
  PF2500   1.571E+02 1.620E+02 2.152E+02 -3.759E+02 8.562E+01
  PF2500   1.572E+02 1.622E+02 2.184E+02 -3.837E+02 8.562E+01
  PF2500   1.543E+02 1.593E+02 2.142E+02 -3.762E+02 8.375E+01
end_S20
#======================================#
"""


if __name__ == '__main__': main()


