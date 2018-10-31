import math
import numpy


def calculateCOM(atoms):
    comX = 0.0
    comY = 0.0
    comZ = 0.0
    sumM = 0.0
    for atom in atoms:
        comX += atom["x"]*atom["m"]
        comY += atom["y"]*atom["m"]
        comZ += atom["z"]*atom["m"]
        sumM += atom["m"]
    return({"x": comX/sumM, "y": comY/sumM, "z": comZ/sumM, "m": sumM})

def subtractCOM(atoms):
    com = calculateCOM(atoms)
    for atom in atoms:
        atom["x"] = atom["x"] - com["x"]
        atom["y"] = atom["y"] - com["y"]
        atom["z"] = atom["z"] - com["z"]
    return(atoms)

def calculateA(atoms):
    sxx = sxy = sxz = syy = syz = szz = 0.0
    for atom in atoms:
        sxx += atom["m"]*atom["x"]*atom["x"]
        sxy += atom["m"]*atom["x"]*atom["y"]
        sxz += atom["m"]*atom["x"]*atom["z"]
        syy += atom["m"]*atom["y"]*atom["y"]
        syz += atom["m"]*atom["y"]*atom["z"]
        szz += atom["m"]*atom["z"]*atom["z"]
    a = [[syy+szz,-sxy,-sxz],[-sxy,sxx+szz,-syz],[-sxz,-syz,sxx+syy]]
    return(a)
    
def calculateBe(a):
    inertia = sorted(list(numpy.linalg.eigh(a, UPLO="L")[0]))
    h   = 6.62608E-34  # Plancks constant J*s
    Na  = 6.02214E+23  # Avogadro's number mol-1
    c   = 299792458    # speed of light m/s
    pi2 = 9.8696044010 # pi²
    B = []
    for I in inertia:
        if I > 0.001:
            Isi = I*(1E-20)/1000                        # convert g/mol*A² to kg*m²/mol
            Be = h/(8*pi2*(Isi/Na))                     # convert kg*m²/mol to Hz
            B.append(Be/1E6)                            # convert to MHz
    return(B)

def getRotationalConstants(atoms):
    return(calculateBe(calculateA(subtractCOM(atoms))))      
