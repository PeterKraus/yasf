import sys


def parse(lines):
    comment = lines.pop(0).strip()
    try:
        line = lines.pop(0).rstrip()
        n_atoms, debug = [int(each) for each in line[0:10].split()]
    except ValueError:
        raise ValueError("ReadSTF: Problem reading number of atoms: \n{0:s}".format(line))
    # ReadZMAT
    atoms = []
    for atom in range(0,n_atoms):
        line = lines.pop(0).rstrip()
        try:
          index, bond_with, angle_with, torsion_with, bond, angle, torsion, mass = line.strip().split()
          index, bond_with, angle_with, torsion_with = int(index), int(bond_with), int(angle_with), int(torsion_with)
          if index <= max([bond_with, angle_with, torsion_with]):
            raise RuntimeError("ReadZMAT: Index of atom {0:d} is <= than referenced atoms.".format(atom+1))
          elif max([bond_with, angle_with, torsion_with]) > len(atoms):
            raise RuntimeError("ReadZMAT: Referenced atoms not present in Z-matrix.")
        except ValueError and RuntimeError:
            raise ValueError("ReadZMAT: Problem reading atom line: \n{0:s}".format(line))
        atoms.append(dict())
        if torsion_with != 0:
            atoms[-1]["torsion_with"] = torsion_with
            atoms[-1]["torsion"] = float(torsion)
        if angle_with != 0:
            atoms[-1]["angle_with"] = angle_with
            atoms[-1]["angle"] = float(angle)
        if bond_with != 0:
            atoms[-1]["bond_with"] = bond_with
            atoms[-1]["bond"] = float(bond)
        atoms[-1]["masses"] = [float(mass)]
    # ReadPARS
    try:
        line = lines.pop(0).rstrip()
        n_par = int(line[32:34].strip())
    except ValueError:
        raise ValueError("ReadPARS: Error reading total number of parameters: \n{0:s}".format(line))
    fitpars = []
    for para in range(0,n_par):
        line = lines.pop(0).rstrip()
        fix = False
        if "FIX" in line[28:32]:
            fix = True
        try:
            atom, parameter, multi = [int(each) for each in line[32:].split()]
            fitpars.append(dict())
            if parameter == 1:
                fitpars[-1]["type"] = "bond"
            elif parameter == 2:
                fitpars[-1]["type"] = "angle"
            elif parameter == 3:
                fitpars[-1]["type"] = "torsion"
            if fix:
                fitpars[-1]["fix"] = True
            else:
                fitpars[-1]["fix"] = False
            fitpars[-1]["atoms"] = [atom]
            if multi > 0:
                line = lines.pop(0).rstrip()
                others = line[32:51].split()
                if len(others) != multi:
                    raise ValueError("ReadPARS: Inconsistent multiple parameter input")
                for each in others:
                    fitpars[-1]["atoms"].append(int(each))
        except ValueError:
            raise ValueError("ReadPARS: Error reading parameter description: \n{0:s}".format(line))
    # ReadCONST
    try:
        line = lines.pop(0).rstrip()
        n_const = int(line[32:34].strip())
    except ValueError:
        raise ValueError("ReadCONST: Error reading total number of rotational constants: \n{0:s}".format(line))
    constants = []
    for const in range(0, n_const):
        line = lines.pop(0).rstrip()
        fit = True
        if "XXX" in line[28:32]:
            fit = False
        try:
            typ, struct, constant = line[32:51].split()
            typ, struct, constant = int(typ), int(struct), float(constant)
            if typ > 9:
                raise RuntimeError("ReadCONST: Rotational constants other than A, B, C, B+C, A+B, B-C, A-(B+C)/2 are not supported!")
            elif typ == 1:
                typ = "A"
            elif typ == 2:
                typ = "B"
            elif typ == 3:
                typ = "C"
            elif typ == 4:
                typ = "B+C"
            elif typ == 5:
                typ = "A+B"
            elif typ == 6:
                typ = "B-C"
            elif typ == 7:
                typ = "(B+C)/2"
            elif typ == 8:
                typ = "(B-C)/2"
            elif typ == 9:
                typ = "A-(B+C)/2"
            if struct > len(constants):
                constants.append(dict())
            constants[struct-1][typ] = {"fit": fit, "value": constant}
        except ValueError:
            raise ValueError("ReadCONST: Error reading rotational constant description: \n{0:s}".format(line))
    # ReadCHNG
    for struct in range(1,len(constants)):
        for atom in atoms:
            atom["masses"].append(atom["masses"][0])
        try:
            line = lines.pop(0).rstrip()
            n_chng = int(line[32:34].strip())
        except ValueError:
            raise ValueError("ReadCHNG: Error reading number of changes: \n{0:s}".format(line))
        for chng in range(0, n_chng):
            line = lines.pop(0).rstrip()
            atom, typ, newval = line[32:51].split()
            atom, typ, newval = int(atom), int(typ), float(newval)
            if typ != 4:
                raise RuntimeError("ReadCHNG: Changing parameters other than mass not supported!")
            atoms[atom-1]["masses"][struct] = newval
    return({"atoms": atoms, "constants": constants, "fitpars": fitpars})
            
            
    
    
    

