#!/usr/bin/python
"""Clangini 14-01-2022
"""
from __future__ import print_function
import sys, os

def get_gaff_params(data):

    charges = []
    atmtypes = []

    in_atom = False
    in_mol = False
    count = 0
    for line in data:
        if "@<TRIPOS>MOLECULE" in line:
            in_mol = True
            count += 1
            continue

        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            in_mol = False
            continue

        if "@<TRIPOS>BOND" in line:
            in_atom = False
            break

        if in_mol:
            count += 1
            if count == 3:
                natom = int(line.split()[0])
                nbond = int(line.split()[1])
                nsub = int(line.split()[2])
            continue

        if in_atom:
            line_split = line.split()
            atmtypes.append(line_split[5])
            charges.append(float(line_split[8]))
            continue
    if len(atmtypes) != natom:
        sys.exit(f"Length of atom list ({len(atmtypes)}) and declared number of atoms ({natom}) differ!")
    
    header = (natom, nbond, nsub)
    return charges, atmtypes, header

def get_mol2_out(orig, charges, header):
    
    in_atom = False
    in_mol = False 
    in_bond = False

    count = 0
    out = ""
    lps = []

    for line in orig:

        if len(line.split()) == 0:
            continue
        if line.strip().startswith("#"):
            continue

        if "@<TRIPOS>MOLECULE" in line:
            in_mol = True
            count += 1
            out += line
            continue

        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            in_mol = False
            atm_idx = 0
            out += line 
            continue

        if "@<TRIPOS>BOND" in line:
            in_atom = False
            in_bond = True 
            bnd_idx = 0 
            out += line
            continue

        if "@<TRIPOS>ALT_TYPE" in line:
            break

        if in_mol:
            count += 1
            if count == 3:
                # natom = int(line.split()[0])
                # nbond = int(line.split()[1])
                out+= "%d %d 0 0 0\n"%(header[0], header[1])
                continue
            if count == 5:
                out += "USER_CHARGES\n"
                out += "Gaff 2 parametrization\n"
                continue

        if in_atom:
            ls = line.split()
            if ls[1].startswith("LP"): # lone pair -> skip it
                lps.append(ls[0])
                continue
            out += ("%7s %-8s%10.4f%10.4f%10.4f %-7s%5s %-11s%7.4f\n") % (int(ls[0]), ls[1], float(ls[2]), float(ls[3]), float(ls[4]), ls[5], "1", "LIG", float(charges[atm_idx]))
            atm_idx += 1
            continue

        if in_bond:
            if bnd_idx >= header[1]:
                break
            ls = line.split()
            if (ls[1] in lps) or (ls[2] in lps): # skip bonds to lone pairs
                continue 
            bnd_idx += 1

        out += line
    
    return out

def main(argv):
    mol2in = argv[0]
    gaffin = argv[1]
    mol2out = argv[2]

    with open(mol2in, "r") as f:
        orig = f.readlines()
    with open(gaffin, "r") as f:
        gaff = f.readlines()
    
    charges, atmtypes, header = get_gaff_params(gaff)

    # i = 1
    # for q, at in zip(charges, atmtypes):
    #     print(i, q, at)
    #     i += 1

    out = get_mol2_out(orig, charges, header)

    # Now we add the CGenFF alternative atom types
    out += "@<TRIPOS>ALT_TYPE\nGAFF_2_ALT_TYPE_SET\nGAFF_2 "
    for i, at in enumerate(atmtypes):
        out += str(i+1)+" "+at+" "
    out += "\n"

    # mol2out = os.path.splitext(mol2in)[0] + "_par.mol2"
    with open(mol2out, "w") as f:
        f.write(out)

if __name__ == "__main__":
    main(sys.argv[1:])
