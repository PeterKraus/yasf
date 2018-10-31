import math

##############################################
# vector functions
##############################################

def vectorDotProduct(a, b):
    return(a["x"]*b["x"] + a["y"]*b["y"] + a["z"]*b["z"])

def matrixDotProduct(M, b):
    return({"x": vectorDotProduct(M["x"], b), "y": vectorDotProduct(M["y"], b), "z": vectorDotProduct(M["z"], b)})

def crossProduct(a, b):
    return({"x": a["y"]*b["z"] - a["z"]*b["y"],
            "y": a["z"]*b["x"] - a["x"]*b["z"],
            "z": a["x"]*b["y"] - a["y"]*b["x"]})

def scalarDivision(a, n):
    return({"x": a["x"]/n, "y": a["y"]/n, "z": a["z"]/n})

def scalarMultiplication(a, n):
    return({"x": a["x"]*n, "y": a["y"]*n, "z": a["z"]*n})

def vectorAddition(a, b):
    return({"x": a["x"] + b["x"], "y": a["y"] + b["y"], "z": a["z"] + b["z"]})
    
def vectorSubtraction(a, b):
    return({"x": a["x"] - b["x"], "y": a["y"] - b["y"], "z": a["z"] - b["z"]})

def vectorNorm(a):
    return(math.sqrt(a["x"]**2 + a["y"]**2 + a["z"]**2))

def rotationMatrix(axis, theta):
    axis = scalarDivision(axis, vectorNorm(axis))
    theta = theta*2*math.pi/360
    a = math.cos(theta/2.0)
    b = -axis["x"] * math.sin(theta/2.0)
    c = -axis["y"] * math.sin(theta/2.0)
    d = -axis["z"] * math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return { "x": {"x": aa+bb-cc-dd, "y": 2*(bc+ad), "z": 2*(bd-ac)},
             "y": {"x": 2*(bc-ad), "y": aa+cc-bb-dd, "z": 2*(cd+ab)},
             "z": {"x": 2*(bd+ac), "y": 2*(cd-ab), "z": aa+dd-bb-cc}}


def zmatToXYZ(zmat):
    xyz = []
    for atom in zmat:
        if "torsion" in atom and len(xyz) >= 3:
            # place D along BC vector
            v_bc = vectorSubtraction(xyz[atom["bond_with"]-1], xyz[atom["angle_with"]-1])
            u_bc = scalarDivision(v_bc,vectorNorm(v_bc))
            d0 = scalarMultiplication(u_bc, atom["bond"])
            # rotate d0 by angle around ab x bc
            ab_cross_bc = crossProduct(vectorSubtraction(xyz[atom["angle_with"]-1],xyz[atom["torsion_with"]-1]),u_bc)
            u_n = scalarDivision(ab_cross_bc, vectorNorm(ab_cross_bc))
            M = rotationMatrix(u_n, 180.0 - atom["angle"] + 1E-20)
            d1 = matrixDotProduct(M, d0)
            # rotate d1 by torsion
            M = rotationMatrix(u_bc, atom["torsion"])
            d2 = matrixDotProduct(M, d1)
            D2 = vectorAddition(xyz[atom["bond_with"]-1],d2)
            xyz.append({"x": D2["x"], "y": D2["y"], "z": D2["z"], "masses": atom["masses"]})
        elif "angle" in atom and len(xyz) == 2:
            # place C along AB vector
            v_ab = vectorSubtraction(xyz[atom["bond_with"]-1], xyz[atom["angle_with"]-1])
            u_ab = scalarDivision(v_ab, vectorNorm(v_ab))
            d0 = scalarMultiplication(u_ab, atom["bond"])
            # rotate d0 by angle around u_y x ab
            y_cross_ab = crossProduct({"x": 0.0, "y": 1.0, "z": 0.0}, u_ab)
            u_n = scalarDivision(y_cross_ab, vectorNorm(y_cross_ab))
            M = rotationMatrix(u_n, 180.0 - atom["angle"] + 1E-20)
            d1 = matrixDotProduct(M, d0)
            D1 = vectorAddition(xyz[atom["bond_with"]-1],d1)
            xyz.append({"x": D1["x"], "y": D1["y"], "z": D1["z"], "masses": atom["masses"]})
        elif "bond" in atom and len(xyz) == 1:
            xyz.append({"x": 0.0, "y": 0.0, "z": atom["bond"], "masses": atom["masses"]})
        elif len(xyz) == 0:
            xyz.append({"x": 0.0, "y": 0.0, "z": 0.0, "masses": atom["masses"]})
        else:
            raise
    return(xyz)
