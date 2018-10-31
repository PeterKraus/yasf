#!/usr/bin/python3
import inputParser as ip
import geoUtils as gu
import rotUtils as ru
import math
import copy
import sys

gridsize = 1
defaultBondAdjust = 0.0025
defaultAngleAdjust = 0.025
iterations = 10
fn = sys.argv[1]

def nameFromMass(mass):
    mass = round(mass)
    if mass == 0:
        return("X ")
    elif mass == 1:
        return("H ")
    elif mass == 12:
        return("C ")
    elif mass == 14:
        return("N ")
    elif mass == 16:
        return("O ")
    elif mass == 19:
        return("F ")
    elif mass == 20:
        return("Ne")
    elif mass == 32:
        return("S ")
    elif mass == 35:
        return("Cl")
    elif mass == 40:
        return("Ar")
    elif mass == 84:
        return("Kr")
    elif mass == 130:
        return("Te")
    elif mass == 132:
        return("Xe")
    

def getRMS(calculated, reference):
    if len(calculated) != 3:
        calculated.append(calculated[0])
        #raise RuntimeError("getRMS: Linear molecules currently not supported!")
    squaredev = []
    for constant in reference:
        if reference[constant]["fit"]:
            if constant == "A":
              dev = (calculated[0] - reference["A"]["value"])
            elif constant == "B":
              dev = (calculated[1] - reference["B"]["value"])
            elif constant == "C":
              dev = (calculated[2] - reference["C"]["value"])
            elif constant == "A+B":
              dev = (calculated[0] + calculated[1] - reference["A+B"]["value"])
            elif constant == "B+C":
              dev = (calculated[1] + calculated[2] - reference["B+C"]["value"])
            elif constant == "(B+C)/2":
              dev = ((calculated[1] + calculated[2])/2 - reference["(B+C)/2"]["value"])
            elif constant == "B-C":
              dev = (calculated[1] - calculated[2] - reference["B-C"]["value"])
            elif constant == "(B-C)/2":
              dev = ((calculated[1] - calculated[2])/2 - reference["(B-C)/2"]["value"])
            elif constant == "A-(B+C)/2":
              dev = (calculated[0] - (calculated[1] + calculated[2]) / 2 - reference["A-(B+C)/2"]["value"])
            squaredev.append(dev**2)
        #print(reference[constant]["value"], calculated, squaredev, sumsq, math.sqrt(sumsq/(nconst+1E-20)))
    return(squaredev)

def getResiduals(cartesians, constants):
    RMSres = []
    for structureIndex in range(0,len(constants)):
        structure = copy.deepcopy(cartesian)
        for atom in structure:
            atom["m"] = atom["masses"][structureIndex]
        calculatedConstants = ru.getRotationalConstants(structure)
        squaredev = getRMS(calculatedConstants, constants[structureIndex])
        for dev in squaredev:
            RMSres.append(dev)
    return(math.sqrt(min(RMSres)), RMSres.index(min(RMSres)),
           math.sqrt(max(RMSres)), RMSres.index(max(RMSres)),
           math.sqrt(sum(RMSres)/len(RMSres)))

# 0) load file
lines = open(fn,"r").readlines()
parsed = ip.parse(lines)

# 1) get average RMS residual from the initial structure
cartesian = gu.zmatToXYZ(parsed["atoms"])
minRMS, mini, maxRMS, maxi, avgRMS = getResiduals(cartesian, parsed["constants"])
inputgeometry = ""
inputgeometry += "{0:d} \n".format(len(cartesian))
inputgeometry += " RMS {0:10.6f} \n".format(maxRMS)
for atom in cartesian:
    inputgeometry += "{0:s} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(nameFromMass(atom["masses"][0]), atom["x"], atom["y"], atom["z"])
open("input.xyz","w").writelines(inputgeometry)
print(" Original RMS:            min res. {0:10.3f} MHz, max res. {1:10.3f} MHz, RMS res. {2:10.3f} MHz".format(minRMS, maxRMS, avgRMS))
print(" Closest structure: {:d}, Furthest structure: {:d}".format(mini, maxi))

# 3) prepare array of changes to Z-Matrix
modifications = []
for par in parsed["fitpars"]:
    if par["fix"] != True:
        modifications.append({"atoms": par["atoms"], "type": par["type"]})
combinations = []
combination = [0]*len(modifications)
for combIndex in range(0,(1+gridsize*2)**len(modifications)):
    combinations.append(copy.deepcopy(combination))
    combination[0] += 1
    if len(combinations) == (1+gridsize*2)**len(modifications):
        break
    for position in range(0,len(modifications)):
        if combination[position] == (1+gridsize*2):
            combination[position] = 0
            combination[position+1] += 1

# 4) loop through changed Z-matrices
converged = False
currentAvgRMS, previousAvgRMS, currentMaxRMS, previousMaxRMS, currentMinRMS, previousMinRMS = 99999, 99999, 99999, 99999, 99999, 99999
oldZMAT = copy.deepcopy(parsed["atoms"])
iteration = [0,0]
factor = 1
while iteration[0] < iterations:
    bondAdjust, angleAdjust = defaultBondAdjust*factor, defaultAngleAdjust*factor
    #print(bondAdjust, angleAdjust)
    while True:
        for combination in combinations:
            newZMAT = copy.deepcopy(oldZMAT)
            for modIndex in range(0,len(modifications)):
                modification = modifications[modIndex]
                if modification["type"] == "bond":
                    offset = (combination[modIndex] - gridsize) * bondAdjust
                else:
                    offset = (combination[modIndex] - gridsize) * angleAdjust
                for atom in modification["atoms"]:
                    newZMAT[atom - 1][modification["type"]] += offset
            cartesian = gu.zmatToXYZ(newZMAT)
            MinRMS, mini, MaxRMS, maxi, AvgRMS = getResiduals(cartesian, parsed["constants"])
            if AvgRMS <= currentAvgRMS:
                currentAvgRMS = AvgRMS
                currentMaxRMS = MaxRMS
                currentMinRMS = MinRMS
                minCart = copy.deepcopy(cartesian)
                minZMAT = copy.deepcopy(newZMAT)
        oldZMAT = minZMAT
        print(" RMS after {0:2d}/{1:2d} iterations: min res. {2:10.3f} MHz, max res. {3:10.3f} MHz, RMS res. {4:10.3f} MHz".format(iteration[0],iteration[1], currentMinRMS, currentMaxRMS, currentAvgRMS))
        print(" Closest structure: {:d}, Furthest structure: {:d}".format(mini, maxi))
        iteration[1] += 1
        #print(previousAvgRMS, currentAvgRMS, previousMaxRMS, currentMaxRMS)
        if round(previousAvgRMS,3) == round(currentAvgRMS,3) and round(previousMaxRMS,3) == round(currentMaxRMS,3) and iteration[1] > 2:
            iteration[0] += 1
            iteration[1] = 0
            bondAdjust = bondAdjust / float(gridsize)
            angleAdjust = angleAdjust / float(gridsize)
            break
        previousAvgRMS = currentAvgRMS
        previousMaxRMS = currentMaxRMS
        
        
       #elif iteration[1] >= 10:
       #    iteration[0] += 1
       #    iteration[1] = 0
       #    break
    geometry = ""
    geometry += "{0:d} \n".format(len(minCart))
    geometry += " min {0:10.3f} max {1:10.3f} MHz, RMS {2:10.3f} MHz\n".format(currentMinRMS, currentMaxRMS, currentAvgRMS)
    for atom in minCart:
        geometry += "{0:s} {1:10.6f} {2:10.6f} {3:10.6f}\n".format(nameFromMass(atom["masses"][0]), atom["x"], atom["y"], atom["z"])
    open(fn[:fn.rindex(".")] + ".yasf.xyz","w").writelines(geometry)
    factor = iteration[0]
    outputzmatrix = ""
    counter = 1
    for atom in minZMAT:
        bond_with, angle_with, torsion_with = 0, 0, 0
        bond, angle, torsion = 0.0, 0.0, 0.0
        mass = atom["masses"][0]
        if counter > 1:
            bond_with, bond = atom["bond_with"], atom["bond"]
        if counter > 2:
            angle_with, angle = atom["angle_with"], atom["angle"]
        if counter > 3:
            torsion_with, torsion = atom["torsion_with"], atom["torsion"]
        outputzmatrix += " {0:3d} {1:3d} {2:3d} {3:3d} {4:12.6f} {5:12.6f} {6:12.6f}   {7:15.8f}\n".format(counter, bond_with, angle_with, torsion_with, bond, angle, torsion, mass)
        counter += 1
    open(fn[:fn.rindex(".")] + ".yasf.zmat","w").writelines(outputzmatrix)


