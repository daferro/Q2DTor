#!/usr/bin/env python

"""
This Python script is used to check the results for Q2DTor
Needed files:
   * 
"""

import os
import shutil

#---------------------------#
# path to Q2DTor executable #
#---------------------------#
#Q2DTor = "/home/david/PyFerro/Q2DTor_v1.1-20181226/Q2DTor.py"
Q2DTor = "/".join(os.path.realpath(__file__).split("/")[:-2])+"/Q2DTor.py"
if not os.path.exists(Q2DTor):
   Q2DTor = "/".join(os.path.realpath(__file__).split("/")[:-2])+"/src/Q2DTor.py"
#---------------------------#


#=================================#
# GENERAL FUNCTION FOR READING    #
#=================================#
def read_lines(filename):
    # read lines
    with open(filename,'r') as ff: lines = ff.readlines()
    # ignore comments
    lines = [line.split("#")[0] for line in lines]
    # strip lines
    lines = [line.strip()       for line in lines]
    # Ignore blank lines
    lines = [line               for line in lines if line != ""]
    return lines
#---------------------------------#
def get_block_of_lines(lines,key1,key2=None):
    record   = False
    selected = []
    for line in lines:
        if key2 is not None and key2 in line: record = False
        if key1 in line: selected,record = [], True
        if record: selected.append(line)
    return selected
#---------------------------------#
def get_startend_blocks(lines):
    keep = False
    blocks = {}
    for line in lines:
        if line.startswith("end_"):
           keep = False
           blocks[system] = slines
        if keep:
           slines.append(line)
        if line.startswith("start_"):
           system = line.split("_")[1]
           keep   = True
           slines = []
    return blocks
#=================================#

#=================================#
# READ REFERENCE FILES            #
#=================================#
def extract_Tlist(blocks):
    Tlist = []
    if "Tlist" in blocks.keys():
       for line in blocks["Tlist"]:
           Tlist += line.split()
       Tlist = [float(T) for T in Tlist]
       Tlist.sort()
    return Tlist
#---------------------------------#
def block_data(lines):
    data = {}
    for line in lines:
        if line.startswith("==>"):
            props = line.split("==>")[1].split()
            continue
        values = line.split()
        for idx,prop in enumerate(props):
            data[prop] = data.get(prop,[])+[float(values[idx])]
    return data
#---------------------------------#
def read_ref_data(filename):
    # read lines
    lines = read_lines(filename)
    # Extract blocks
    blocks = get_startend_blocks(lines)
    lines = None
    # Temperatures
    Tlist = extract_Tlist(blocks)
    blocks.pop("Tlist")
    # Data for each system
    data = {}
    for system,lines in sorted(blocks.items()):
        data[system] = block_data(lines)
        blocks.pop(system)
    return Tlist, data
#---------------------------------#
def write_ref_data(Tlist,data,filename=None):
    string = ""
    # temperatures
    string += "start_Tlist\n"
    for idx in range(0,len(Tlist),5):
        stringT = "  ".join(["%6.1f"%T for T in Tlist[idx:idx+5]])
        string += "    "+stringT+"\n"
    string += "end_Tlist\n"
    string += "\n"
    # systems
    for system,sdata in sorted(data.items()):
        string += "start_%s\n"%system
        for prop,values in sorted(sdata.items()):
            string += "  ==> %s \n"%prop
            for value in values:
                string += "  %12.5E \n"%value
        string += "end_%s\n"%system
        string += "\n"
    if filename is not None:
       with open(filename,'w') as ff: ff.write(string)
    else:
       print "\n"+string
#=================================#

#=================================#
# READ Q2DTOR OUT FILE            #
#=================================#
def get_table(lines,text):
    save = False
    data  = {}
    for line in lines:
        if save:
           if idx > 0 and line.strip() == "": break
           idx += 1
           if idx == 2:
               props = [prop.strip() for prop in line.split("|")]
           if idx >3:
               values = line.split("|")
               for prop,value in zip(props,values):
                   data[prop] = data.get(prop,[]) + [float(value)]
        if text in line:
           save = True
           idx  = 0
           data = {}
    return data
#---------------------------------#
def data_from_q2dtor_out(outfile):

    ldata  = ["rv1WHO*","rvMSHO*","rvE2DT*","Q2DNS*","QEHR*","QTorClas"]
    ldata += ["rv1WHO" ,"rvMSHO" ,"rvE2DT" ,"Q2DNS" ,"QEHR" ,"Fq2DNS"  ]
    ldata += ["Qtr"    ,"Qele"                                         ]
    ldata += ["Q1WHO"  ,"QMSHO"  ,"QE2DT"                              ]
    ldata += ["U1WHO"  ,"H1WHO"  ,"S1WHO"  ,"G1WHO" ,"Cp1WHO"          ]
    ldata += ["UMSHO"  ,"HMSHO"  ,"SMSHO"  ,"GMSHO" ,"CpMSHO"          ]
    ldata += ["UE2DT"  ,"HE2DT"  ,"SE2DT"  ,"GE2DT" ,"CpE2DT"          ]
    Tlist = []
    data  = {key:[] for key in ldata}

    # read file
    lines  = read_lines(outfile)

    # read rovib pfns with regard to ZPE
    key1  = "(a) Rovibrational partition functions using as zero of energy"
    key2  = "(b) Rovibrational partition functions obtained using as zero of energy"
    block = get_block_of_lines(lines,key1,key2)

    subblock = get_block_of_lines(block,"T (K)  |  rv(1WHO)  |","Components of E2DT")
    for line in subblock[2:]:
        T,rv1WHO,rvMSHO,rvE2DT,ratio = [float(val) for val in line.split("|")]
        Tlist.append(T)
        data["rv1WHO*"].append(rv1WHO)
        data["rvMSHO*"].append(rvMSHO)
        data["rvE2DT*"].append(rvE2DT)

    subblock = get_block_of_lines(block,"T (K)  |    2DNS    |  Tors")
    for line in subblock[2:]:
        T,Q2DNS,QTorClas,QEHR = [float(val) for val in line.split("|")]
        data["Q2DNS*"].append(Q2DNS)
        data["QTorClas"].append(QTorClas)
        data["QEHR*"].append(QEHR)

    # read rovib pfns with regard to the PES bottom
    key1  = "(b) Rovibrational partition functions obtained using as zero of energy"
    key2  = "Translational and electronic partition functions"
    block = get_block_of_lines(lines,key1,key2)

    subblock = get_block_of_lines(block,"T (K)  |  rv(1WHO)  |  rv(MSHO)","Components of E2DT")
    for line in subblock[2:]:
        T,rv1WHO,rvMSHO,rvE2DT,ratio = [float(val) for val in line.split("|")]
        data["rv1WHO"].append(rv1WHO)
        data["rvMSHO"].append(rvMSHO)
        data["rvE2DT"].append(rvE2DT)

    subblock = get_block_of_lines(block,"T (K)  |    2DNS    |  Tors")
    for line in subblock[2:]:
        T,Q2DNS,QTorClas,QEHR,Fq = [float(val) for val in line.split("|")]
        data["Q2DNS" ].append(Q2DNS)
        data["QEHR"  ].append(QEHR)
        data["Fq2DNS"].append(Fq)

    # read translational and electronic partition functions
    key1  = "Translational and electronic partition functions"
    key2  = "Total partition functions, from rovib"
    block = get_block_of_lines(lines,key1,key2)
    for line in block[3:]:
        T, Qtr, Qele = [float(val) for val in line.split("|")]
        data["Qtr" ].append(Qtr)
        data["Qele"].append(Qele)
    
    # read total partition functions with regard to PES bottom
    key1  = "Total partition functions, from rovibrational partition functions in (b)"
    key2  = "Tables were stored at "
    block = get_block_of_lines(lines,key1,key2)
    for line in block[3:]:
        T, Q1WHO, QMSHO, QE2DT = [float(val) for val in line.split("|")]
        data["Q1WHO"].append(Q1WHO)
        data["QMSHO"].append(QMSHO)
        data["QE2DT"].append(QE2DT)

    # Get thermodynamics for 1WHO
    key1  = "Thermodynamic functions for 1WHO partition function"
    key2  = "Thermodynamic functions for MSHO partition function"
    block = get_block_of_lines(lines,key1,key2)
    for line in block[3:]:
        T, U, H, S, G, Cp = [float(val) for val in line.split("|")]
        data[ "U1WHO"].append( U)
        data[ "H1WHO"].append( H)
        data[ "S1WHO"].append( S)
        data[ "G1WHO"].append( G)
        data["Cp1WHO"].append(Cp)

    # Get thermodynamics for MSHO
    key1  = "Thermodynamic functions for MSHO partition function"
    key2  = "Thermodynamic functions for E2DT partition function"
    block = get_block_of_lines(lines,key1,key2)
    for line in block[3:]:
        T, U, H, S, G, Cp = [float(val) for val in line.split("|")]
        data[ "UMSHO"].append( U)
        data[ "HMSHO"].append( H)
        data[ "SMSHO"].append( S)
        data[ "GMSHO"].append( G)
        data["CpMSHO"].append(Cp)

    # Get thermodynamics for E2DT
    key1  = "Thermodynamic functions for E2DT partition function"
    key2  = "Internal Energy"
    block = get_block_of_lines(lines,key1,key2)
    for line in block[3:]:
        T, U, H, S, G, Cp = [float(val) for val in line.split("|")]
        data[ "UE2DT"].append( U)
        data[ "HE2DT"].append( H)
        data[ "SE2DT"].append( S)
        data[ "GE2DT"].append( G)
        data["CpE2DT"].append(Cp)

    return Tlist, data
#=================================#



#=================================#
# GEOMETRIES (XYZ FILES)          #
#=================================#
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
#---------------------------------#
def num_sp_in_file(filesplist):
    if not os.path.exists(filesplist): return -1, -1
    with open(filesplist,'r') as ff: lines = ff.readlines()
    nsp = 0
    countNO  = 0
    countYES = 0
    for line in lines:
        line = line.split("#")[0].strip()
        if len(line.split()) != 6: continue
        if "NO"  in line: countNO  += 1
        if "YES" in line: countYES += 1
    return countYES, countNO
#---------------------------------#
def get_num_sp(CASE):
    dnumsp        = {}
    systems  = "S01,S02,S03,S04,S05,S06,S07,S08,S09,S10".split(",")
    systems += "S11,S12,S13,S14,S15,S16,S17,S18,S19,S20".split(",")
    if CASE == "3":
       numsp = [10,14,13,8,11,10,15,7,20,6,6,6,11,15,8,8,3,14,8,6]
       #dnumsp["S21"] = 8
       #dnumsp["S22"] = 7
    else:
       numsp = [10,10,13,8,11,10,13,7,20,7,7,6,12,13,6,8,3,14,8,6]
    for key,value in zip(systems,numsp): dnumsp[key] = value
    return dnumsp
#=================================#


#=================================#
# OPTION (1) - CREATE INPUT FILES #
#=================================#
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
#---------------------------------#
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
#---------------------------------#
def option_createINPs(FOLDER,FREQSCAL,LEVEL):

    # Temperatures
    Tlist = "    100.0  150.0  200.0      "
    Tlist+= "    250.0  300.0  400.0      "
    Tlist+= "    500.0  700.0 1000.0      "
    Tlist+= "   1500.0 2000.0 2500.0      "
    Tlist = [float(T) for T in Tlist.split()]

    # Geometries
    dict_SXX = get_geometries()

    # Create folders
    print "    - Creating test files in: %s"%FOLDER
    if not os.path.exists(FOLDER): os.mkdir(FOLDER)

    # Go system by system
    for SXX in sorted(dict_SXX.keys()):
        print "      * %s"%SXX

        # Create SXX folder
        if os.path.exists(FOLDER+SXX):
           print "        folder '%s' already exists..."%(FOLDER+SXX)
           continue
        if not os.path.exists(FOLDER+SXX): os.mkdir(FOLDER+SXX)

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

        # apply case
        level    = LEVEL
        freqscal = FREQSCAL

        # Get data for SXX
        xyz_string, torsion1, torsion2, tsigma1, tsigma2, symmetry, ttype, charge, multiplicity = dict_SXX[SXX]
        charge  = str(charge)
        multiplicity = str(multiplicity)

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
        thefile = open(FOLDER+SXX+"/"+SXX+".xyz",'w')
        thefile.write(xyz_string)
        thefile.close()

        # Write input file
        thefile = open(FOLDER+SXX+"/"+SXX+".inp",'w')
        thefile.write(inp_string)
        thefile.close()
    print
#=================================#


#=================================#
# OPTION (2) - EXECUTE Q2DTOR     #
#=================================#
def option_execute(FOLDER,CASE):
    # ask user for option
    print
    print "        Option?"
    print "        (1) --init"
    print "        (2) --pes"
    print "        (3) --fourier"
    print "        (4) --findsp"
    print "        (5) --optsp"
    print "        (6) --tor2dns"
    print "        (7) --rovibpf thermo"
    print "        (8) --pdf"
    print "        --------------------------"
    print "        (9) check number of SP"
    print "        --------------------------"
    print "        use .. to exit"
    print "        --------------------------"
    print "        type option and system to only"
    print "        apply to one of them. Example:"
    print "          >> 2 S01                    "
    print "        --------------------------"
    # Geometries
    dict_SXX = get_geometries()
    # PWD and nohup file
    PWD = os.getcwd()+"/"
    nohupfile = PWD+"nohup.out"
    # deal with answer
    answer = raw_input("        >> ").strip()
    if   len(answer.split()) == 1: systems = sorted(dict_SXX.keys())
    elif len(answer.split()) == 2: systems = answer.split()[1:2]
    else                         : systems = answer.split()[1:]
    answer = answer.split()[0]
    # Go system by system
    while True:
       for SXX in systems:
           SXXout    = PWD+FOLDER+SXX+"/"        +SXX+".out"
           SXXsplist = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".splist"
           # generate info line about SP
           nYES,nNO  = num_sp_in_file(SXXsplist)
           nsp = nNO+nYES
           NSP = get_num_sp(CASE)[SXX]
           strnsp = "        %3s: %2i vs %2i (YES,NO = %2i,%2i) "%(SXX,NSP,nsp,nYES,nNO)
           if NSP != nsp or nNO != 0: strnsp += "<== WARNING"
           else                     : strnsp += "<== OK"
           #----------------#
           # check num SP   #
           #----------------#
           if answer in "5679":
               case1 = (answer == "9")
               case2 = (answer == "5"  and NSP != nsp)
               case3 = (answer in "67" and (nNO !=  0 ))
               case4 = (answer in "67" and (NSP != nsp))
               if case1 or case2 or case3 or case4:
                  print strnsp; continue
           #----------------#
           # execute Q2DTor #
           #----------------#
           if answer in "12345678":
              if   answer == "1":
                 if os.path.exists(SXXout): os.remove(SXXout)
                 command = "python %s %s --print --init"%(Q2DTor,SXX)
                 if FOLDER == "ORCA/": command += " orca"
              elif answer == "8":
                 command = "python %s %s --pdf"%(Q2DTor,SXX)
              else:
                 # read SXX.out if exists
                 if not os.path.exists(SXXout): continue
                 with open(SXXout,'r') as out: lines = out.readlines()
                 lines = "\n".join(lines)
                 # check what was executed
                 endopt   = "End of execution of Q2DTor with %s"
                 startopt = "Executing Q2DTor with option: %s"
                 go2 = endopt%"--init"    in lines and  startopt%"--pes"     not in lines
                 go3 = endopt%"--pes"     in lines and  startopt%"--fourier" not in lines
                 go4 = endopt%"--fourier" in lines and  startopt%"--findsp"  not in lines
                 go5 = endopt%"--findsp"  in lines and  startopt%"--optsp"   not in lines
                 go6 = endopt%"--optsp"   in lines and  startopt%"--tor2dns" not in lines
                 go7 = endopt%"--tor2dns" in lines and  startopt%"--rovibpf" not in lines
                 # command for each option
                 command = "nohup python %s %s "%(Q2DTor,SXX)
                 if   answer == "2" and go2: command += " --pes            "
                 elif answer == "3" and go3: command += " --fourier        "
                 elif answer == "4" and go4: command += " --findsp         "
                 elif answer == "5" and go5: command += " --optsp          "
                 elif answer == "6" and go6: command += " --tor2dns        "
                 elif answer == "7" and go7: command += " --rovibpf thermo "
                 else                      : continue
                 command += " 1> %s 2> %s &"%(nohupfile,nohupfile)
              # execute
              if "nohup" in command: print "        %3s: executing Q2DTor (with nohup)"%SXX
              else                 : print "        %3s: executing Q2DTor"%SXX
              #print command
              os.chdir(PWD+FOLDER+SXX)
              os.system(command)
              os.chdir(PWD)
           #------------#
           # do nothing #
           #------------#
           else: continue
       # break infinite loop?
       if answer in ["a","b","c","d","e","f","g","h"]:
          answer = raw_input("        >> ").strip()
          continue
       break
#=================================#


#=================================#
# OPTION (3) - CHECK DATA         #
#=================================#
def max_abs_diff(v1,v2):
    if len(v1) != len(v2): return -1
    mdiff  = max([abs(v2-v1) for v1,v2 in zip(v1,v2)])
    return mdiff
#---------------------------------#
def max_rel_diff(v1,v2):
    if len(v1) != len(v2): return -1
    mdiff  = 100*max([abs(v2-v1)/v1 for v1,v2 in zip(v1,v2)])
    return mdiff
#---------------------------------#
def option_check(FOLDER,root,CASE):
    EPS_REL = 2.5
    EPS_ABS = 0.25
    blank   = "      "
    # Geometries
    dict_SXX = get_geometries()

    # Reference file
    if root:
       print "    Reference file? "
       reffile = raw_input("    >> ").strip()
    else:
       if CASE == "1": reffile = "Q2DTor_RefData_GAUSSIAN.txt"
       if CASE == "2": reffile = "Q2DTor_RefData_ORCA.txt"
    if not os.path.exists(reffile): print "    File '%s' does not exist!"%reffile; exit()
    Tlist_ref,data_ref = read_ref_data(reffile)

    # Go system by system
    PWD = os.getcwd()+"/"
    DATA = {}
    for SXX in sorted(dict_SXX.keys()):
        string = ""
        if not os.path.exists(PWD+FOLDER+SXX): continue
        os.chdir(PWD+FOLDER+SXX)
        outfile = PWD+FOLDER+SXX+"/"+SXX+".out"
        # read data
        Tlist, data = data_from_q2dtor_out(outfile)
        # Compare temperatures
        maxdiff = max_rel_diff(Tlist_ref,Tlist)
        if maxdiff == -1: exit("Number of temperatures differs!")
        if maxdiff > 1.0: print "kjgsadkjgasflh..."
        # common data
        variables = list(set(data.keys()).intersection( data_ref[SXX].keys()))
        # Compare data
        for variable in variables:
            vec1 = data_ref[SXX][variable]
            vec2 = data[variable]
            if variable[0] in  ["U","H","G","S","C"]:
                mdiff = max_abs_diff(vec1,vec2)
                if mdiff == -1: exit("Number of %s differs!"%variable)
                if mdiff > EPS_ABS: string += blank+" %-9s: %.1f\n"%(variable,mdiff)
            else:
                mdiff = max_rel_diff(vec1,vec2)
                if mdiff == -1: exit("Number of %s differs!"%variable)
                if mdiff > EPS_REL: string += blank+" %-9s: %.1f%%\n"%(variable,mdiff)
        # Print data
        fstring  = blank+"=============================\n"
        fstring += blank+"|            %3s            |\n"%SXX
        fstring += blank+"-----------------------------\n"
        if string != "":
           fstring += string
        else:
           fstring += blank+" Rel.Err.(pfns) < %.1f%%\n"%EPS_REL
           fstring += blank+" Abs.Err.(ther) < %.2f units\n"%EPS_ABS
        fstring += blank+"=============================\n"
        print fstring
        
    os.chdir(PWD)
#=================================#


#=================================#
# OPTION (4) - GENERATE DATA FILE #
#=================================#
def option_gendata(FOLDER):
    # Geometries
    dict_SXX = get_geometries()
    # Go system by system
    PWD = os.getcwd()+"/"
    DATA = {}
    for SXX in sorted(dict_SXX.keys()):
        os.chdir(PWD+FOLDER+SXX)
        outfile = PWD+FOLDER+SXX+"/"+SXX+".out"
        print "       reading %s"%outfile
        Tlist, data = data_from_q2dtor_out(outfile)
        DATA[SXX] = data
    os.chdir(PWD)
    print
    # Write file
    filename = raw_input("      Name of file with data? ").strip()
    if filename in "": return
    if os.path.exists(filename):
        print "      File '%s' already exists!"%filename
        answer = raw_input("      Overwrite (y/N)? ").strip().lower()
        if answer not in ["y","yes"]: return
    print "      Writing '%s'..."%filename
    write_ref_data(Tlist,DATA,filename)
    print
#=================================#

#=================================#
# OPTION (5) - BACKUP DATA        #
#=================================#
def option_backup(FOLDER,BKFOLDER):
    PWD = os.getcwd()+"/"
    # Geometries
    dict_SXX = get_geometries()
    # backup folder
    print "    Backup folder: '%s'"%BKFOLDER
    if not os.path.exists(BKFOLDER): os.mkdir(BKFOLDER)
    # Go system by system
    for SXX in sorted(dict_SXX.keys()):
        # files to backup
        outfile = PWD+FOLDER+SXX+"/"+SXX+".out"
        ics     = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".ics"
        pes     = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".pes"
        vfit    = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".vfit"
        dfit    = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".dfit"
        sqrt    = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".sqrt"
        evals   = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".evals"
        splist  = PWD+FOLDER+SXX+"/IOfiles/"+SXX+".splist"
        ff      = [outfile,ics,pes,vfit,dfit,sqrt,evals,splist]
        # gts files
        DIRGTS = PWD+FOLDER+SXX+"/IOfiles/SP/"
        if os.path.exists(DIRGTS):
           ff += [DIRGTS+gts for gts in os.listdir(DIRGTS) if gts.endswith(".gts")]
        # copy them
        print "       copying files of %s"%SXX
        for fi in ff:
            command = "cp %s %s/."%(fi,PWD+BKFOLDER)
            if os.path.exists(fi): os.system(command)
        print
#=================================#



#=================================#
# OPTION (6) - COPY FROM BACKUP   #
#=================================#
def option_copyfrombackup(FOLDER,BKFOLDER):
    # ask user for option
    print
    print "        Option?"
    print "        (a) copy PES from backup"
    print "        (b) copy VFIT from backup"
    print "        (c) copy SPLIST from backup"
    print "        (d) copy GTSs from backup"
    print "        (e) copy DFIT from backup"
    print "        (f) copy EVALS from backup"
    print "        (g) copy SQRT from backup"
    print "        --------------------------"
    print "        use .. to exit"
    print "        --------------------------"
    print "        type option and system to only"
    print "        apply to one of them. Example:"
    print "        >> a S01                    "
    print "        --------------------------"
    print
    # Geometries and PWD
    dict_SXX = get_geometries()
    PWD = os.getcwd()+"/"
    # Go system by system
    while True:
       # deal with answer
       answer = raw_input("        >> ").strip()
       if   len(answer.split()) == 1: systems = sorted(dict_SXX.keys())
       elif len(answer.split()) == 2: systems = answer.split()[1:2]
       else                         : systems = answer.split()[1:]
       answer = answer.split()[0]
       # go system by system
       for SXX in systems:
           # copy individual files
           if   answer in "abcefg":
              # select file
              if   answer == "a": ff = "%s.pes"%SXX
              elif answer == "b": ff = "%s.vfit"%SXX
              elif answer == "c": ff = "%s.splist"%SXX
              elif answer == "e": ff = "%s.dfit"%SXX
              elif answer == "f": ff = "%s.evals"%SXX
              elif answer == "g": ff = "%s.sqrt"%SXX
              # copy file
              file1 = "%s%s"%(BKFOLDER,ff)
              file2 = "%s/%s/IOfiles/%s"%(FOLDER,SXX,ff)
              if os.path.exists(file1): shutil.copyfile(file1,file2)
              else: print "        '%s' does not exist..."%file1
           # copy gts files
           elif answer == "d":
              # list gts files
              gtss = [gts for gts in os.listdir(PWD+BKFOLDER) if gts.endswith(".gts")]
              gtss = [gts for gts in gtss if gts.startswith(SXX)]
              # create SP/ folder
              SPfolder = PWD+"%s/%s/IOfiles/SP/"%(FOLDER,SXX)
              if not os.path.exists(SPfolder): os.mkdir(SPfolder)
              # copy files
              for gts in gtss:
                  gts = PWD+BKFOLDER+gts
                  os.system("cp %s %s/."%(gts,SPfolder))
           continue
       # break infinite loop?
       if answer not in "abcdefgh": return
#=================================#

#===========================================#
#             THE MAIN FUNCTION             #
#===========================================#
def get_case(root):
    ibs = "    "
    # Select case
    print ibs+"Select case:"
    print ibs+" (1) HF/sto-3g (GAUSSIAN)"
    print ibs+" (2) HF/sto-3g (ORCA)"
    if root:
       print ibs+" (3) MPWB1K (GAUSSIAN)"
       nn = 3
    else:
       nn = 2
    print
    CASE = raw_input(ibs+"your choice: ").strip()
    print
    if CASE not in [str(i+1) for i in range(nn)]: exit()
    return CASE
#-------------------------------------------#
def get_action(root):
    ibs = "    "
    # Select action
    print ibs+"Select action:"
    print ibs+" (1) Create input files"
    print ibs+" (2) Check results against refdata file"
    print ibs+" (3) Execute Q2DTor"
    if root:
       print ibs+" (4) Generate reference-data file"
       print ibs+" (5) Backup data"
       print ibs+" (6) Copy data from backup"
       nn = 6
    else:
       nn = 3
    print ibs+" ..  to exit"
    print
    ACTION = raw_input(ibs+"your choice: ").strip()
    print
    if ACTION not in [str(i+1) for i in range(nn)]: exit()
    return ACTION
#-------------------------------------------#
def case_to_data(CASE):
    if   CASE == "1":
        FOLDER   = "GAUSSIAN/"
        BKFOLDER = "backup_GAUSSIAN/"
        FREQSCAL = "1.000"
        LEVEL    = "hf 3-21g"
    elif CASE == "2":
        FOLDER   = "ORCA/"
        BKFOLDER = "backup_ORCA/"
        FREQSCAL = "1.000"
        LEVEL    = "hf 3-21g"
    elif CASE == "3":
        FOLDER   = "MPWB1K/"
        BKFOLDER = "backup_MPWB1K/"
        FREQSCAL = "0.964"
        LEVEL    = "mpwb95/6-31+G(d,p) IOp(3/76=0560004400)"
        LEVEL   += " int=ultrafine"
    return FOLDER, BKFOLDER, FREQSCAL, LEVEL
#-------------------------------------------#
def main(root=False):
    print
    print " << Test creator for Q2DTor >>"
    print
    # check Q2DTor dir
    print "    Q2DTOR PATH: '%s'"%Q2DTor
    if not os.path.exists(Q2DTor):
       print "    Unable to find it..."
       print
       return
    print
    # Get CASE
    CASE = get_case(root)
    FOLDER, BKFOLDER, FREQSCAL, LEVEL = case_to_data(CASE)
    while True:
          # Get ACTION
          ACTION = get_action(root)
          # Act according to action
          if   ACTION ==  "1": option_createINPs(FOLDER,FREQSCAL,LEVEL)
          elif ACTION ==  "2": option_check(FOLDER,root,CASE)
          elif ACTION ==  "3": option_execute(FOLDER,CASE)
          elif ACTION ==  "4": option_gendata(FOLDER)
          elif ACTION ==  "5": option_backup(FOLDER,BKFOLDER)
          elif ACTION ==  "6": option_copyfrombackup(FOLDER,BKFOLDER)
          else               : break
          print
#===========================================#

if __name__ == '__main__': main(False)

