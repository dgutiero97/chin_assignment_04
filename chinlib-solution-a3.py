import argparse
import os, sys

from rdkit import Chem


smiles = ""
debug_mode = False


# Perform relaxation step to assign EC labels
def morgan_relaxation(mol):
    pass # TODO


# Derive canonical numbering based on final EC labelling
def morgan_enumeration(mol):
    pass # TODO


# Assign a custom ID to the atoms in a given molecule
def assign_custom_atom_id(mol, canonical):

    if canonical:
        # Use Morgan algorithm to deriva a canonical graph numbering
        morgan_relaxation(mol)
        morgan_enumeration(mol)
    else:
        # Most simple strategy: just copy RDKit internal atom ID
        for atom in mol.GetAtoms():
            atom.SetProp("CID", str(atom.GetIdx()))


# Return SMILES primitive for a given bond
def get_bond_symbol(bond):
    if bond.GetBondTypeAsDouble() == 1:
        return ""
    elif bond.GetBondTypeAsDouble() == 2:
        return "="
    elif bond.GetBondTypeAsDouble() == 3:
        return "#"
    elif bond.GetBondTypeAsDouble() == 1.5:
        return ":"
    else:
        print("Bond of unknown type: EXIT")
        sys.exit(1)


# Return SMILES primitive for a given atom
def get_atom_symbol(atom):
    if atom.GetFormalCharge() > 0:
        return "[" + atom.GetSymbol() + "+" + str(atom.GetFormalCharge()) + "]"
    if atom.GetFormalCharge() < 0:
        return "[" + atom.GetSymbol() + str(atom.GetFormalCharge()) + "]"

    if atom.GetSymbol() in ["B", "b", "C", "c", "N", "n", "O", "o", "S", "s", "P", "p", "F", "Cl", "Br", "I"]:
        return atom.GetSymbol()
    else:
        return "[" + atom.GetSymbol() + "]"


# Traverse molecular graph in depth-first order
def mol_dft(atom):
    global smiles

    # Atom has already been visited
    if atom.HasProp("VISITED"):
        return

    # Tag atom as visited and append element symbol to SMILES string
    atom.SetBoolProp("VISITED", True)
    smiles += get_atom_symbol(atom)

    # Proceed if atom has any neighbors
    if atom.GetDegree():

        # Get sorted list of incident bonds
        # that have not been visited yet
        bonds = [bond for bond in atom.GetBonds() if not bond.HasProp("VISITED")]
        bonds.sort(key=lambda x: int(x.GetOtherAtom(atom).GetProp("CID")))

        # Iterate over all remaining bonds
        for i in range(len(bonds)):

            # Tag bond as visited
            bonds[i].SetBoolProp("VISITED", True)

            # Branch opening if not last bond
            if i < len(bonds) - 1:
                smiles += "("

            # Append bond symbol to SMILES string
            smiles += get_bond_symbol(bonds[i])

            # Recursive call of mol_dft with partner atom
            mol_dft(bonds[i].GetOtherAtom(atom))

            # Branch closing if not last bond
            if i < len(bonds) - 1:
                smiles += ")"

    return


# Generate a SMILES string for a given molecule
def generate_smiles(mol):
    global smiles

    # Clear SMILES string
    smiles = ""

    # Get list of atoms
    atoms = mol.GetAtoms()

    # Start at the first atom
    mol_dft(atoms[0])

    # Make sure that disconnected structures are recognized
    for atom in atoms:
        if not atom.HasProp("VISITED"):
            # Disconnected component identified
            smiles += "."

            # Start SMILES generation
            mol_dft(atom)




# ----------------------------------------------------------
# Main script
# ----------------------------------------------------------


# Command-line argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("i", help="SDF MOL input file")
parser.add_argument("o", help="CSV output file")
parser.add_argument("-d", "--debug", action="store_true", help="Run in debug mode")
parser.add_argument("-o", "--overwrite", action="store_true", help="Overwrite existing output file")

args = parser.parse_args()

if not os.path.isfile(args.i):
    parser.print_help(sys.stderr)
    sys.exit(1)
if not args.i.endswith(".sdf"):
    parser.print_help(sys.stderr)
    sys.exit(1)
if os.path.isfile(args.o) and not args.overwrite:
    print("Output file already exists.")
    print("To overwrite use '--overwrite'.")
    sys.exit(1)

# Set debug mode according to presence of command line flag
debug_mode = args.debug


# Read SD input file
file_i = Chem.SDMolSupplier(args.i)

# Open output file
file_o = open(args.o, "w")
file_o.write("mol_id\tSMILES\n")

# Iterate over molecules
mol_id = 1
for mol in file_i:
    print("-- Processing molecule " + str(mol_id))

    # Generate SMILES
    assign_custom_atom_id(mol, False)
    generate_smiles(mol)

    # Append SMILES to output file
    out_str = "{}\t{}".format(format(mol_id,'05d'), smiles)
    file_o.write(out_str + "\n")

    mol_id += 1

file_o.close()
