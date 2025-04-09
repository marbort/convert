from rdkit import Chem
from rdkit.Chem import Draw
import sys
import pubchempy as pcp
from rdkit.Chem.Draw import rdMolDraw2D



mols=[]
smiles=[]
names=[]
numbers=[]
with open(sys.argv[1]) as f:
    for line in f:
        names.append(line.split(',')[0])
        numbers.append(line.split(',')[1])
        result = pcp.get_compounds(line.split(',')[0], namespace="name")
        print(line.split(',')[0], result)
        smiles.append(result[0].isomeric_smiles)
names_numbers=[f"{name} {number}" for name,number in zip(names,numbers)]    

# Get SMILES
for mol in smiles:
    mols.append(Chem.MolFromSmiles(mol) if mol else None)

print("SMILES:", smiles)

wdt=len(mols)*150
hgt=600
# Set up the drawer
drawer = rdMolDraw2D.MolDraw2DSVG(wdt, hgt, 300, 300)  # (canvas width, height, mol width, height)

# Customize drawing options
options = drawer.drawOptions()
options.legendFontSize = 25  # Increase legend font size

# Draw molecules with names
drawer.DrawMolecules(mols, legends=names_numbers)
drawer.FinishDrawing()

# Save the image
with open("molecules_with_large_text.svg", "w") as f:
    f.write(drawer.GetDrawingText())


