import argparse
import os, sys
import argparse
import os, sys
from rdkit.Chem.rdchem import Mol, Atom, Bond
from rdkit import Chem


smiles = ""
debug_mode = False
idx_can = "idx_canonical"
score = "score"
max_score_idx = "maximun_nboors_score_atom_index"
canon_label = "canonical_order_label"

# Relaxation - Initialization
def m_relax_initialize(mol: Mol):
    mol.SetIntProp("i",1)
    mol.SetIntProp("c",1)
    for atom in mol.GetAtoms():
        atom.SetIntProp(score, 1)

#Neighbors addition
def neighbors_addition(atom : Atom):
    neighbors = atom.GetNeighbors()
    sum = 0
    for neighbor in neighbors:
        sum = sum + neighbor.GetProp("i")
    atom.SetIntProp(score, sum)

# Relaxation - Neighbors addition
def m_relax(mol: Mol):
    prev_Mol : Mol
    while(True):
        #Set the current number of tags to 0        
        C = []
        #Iterate all atoms, perform the addition, and append distinct scores
        for atom in mol.GetAtoms():
            neighbors_addition(atom)
            if atom.GetProp(score) not in C:
                C.append(atom.GetProp(score))
        
        #If the number of current tags does not increases, regarding number of previous tags
        if C.count < mol.GetProp("c"):
            #Use previous molecule iteration
            mol = prev_Mol
            #Finalize the algorithm
            break

        #Use current number of tags as previous number of tags
        mol.SetIntProp("c",C.count)
        #Increase the number of the iteration
        mol.SetIntProp("i", mol.GetIntProp("i")+1)
        #Save current molecule as the previous molecule
        prev_Mol = mol.__copy__()

# Relaxation - General handler
def morgan_relax_handler(mol : Mol):
    
    #Initialize
    m_relax_initialize(mol)
    #Relaxation
    m_relax(mol)

#Initialization phase of Morgan Canonical enumeration phase
def morg_enum_initialization(mol : Mol):

    max = 0;
    mol.SetIntProp("c",2)

    atom: Atom
    for atom in mol.GetAtoms():
        atom.SetIntProp(idx_can, 0)
        if max > atom.GetIntProp(score):
            max = atom.GetIntProp(score)
            Mol.SetIntProp(max_score_idx,atom.GetIdx)
    
    mol.GetAtomWithIdx(mol.GetIntProp(max_score_idx)).SetIntProp(canon_label, 1)

#TODO: AQUI
#Iterative phase of Morgan Canonical enumeration phase
def morg_enum(mol : Mol):
    
    #Obtain the atom with max score
    atom: Atom = mol.GetAtomWithIdx(mol.GetIntProp(max_score_idx))

    #Find the best neighboor of all of them
    best_nboor = None
    best_value = None
    for nboor_atom in atom.GetNeighbors():
        #Initialize
        if best_value is None:
            best_value = nboor_atom.GetIntProp(score)
            best_nboor = nboor_atom
        #Find best score
        if nboor_atom.GetIntProp(score) > best_value:
            best_value = nboor_atom.GetIntProp(score)
            best_nboor = nboor_atom
        #Check by other conditions
        elif nboor_atom.GetIntProp(score) == best_value:
            pass #TODO:
        #In case we have full equality
        else:
            pass #TODO: pick at random





# Derive canonical numbering based on final EC labelling
def morgan_enum_handler(mol : Mol):
    #Initialize the labels
    morg_enum_initialization(mol)
    #Perform the algorithm
    morg_enum(mol)


# Assign a custom ID to the atoms in a given molecule
def assign_custom_atom_id(mol, canonical):

    if canonical:
        # Use Morgan algorithm to deriva a canonical graph numbering
        morgan_relax_handler(mol)
        morgan_enum_handler(mol)
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
