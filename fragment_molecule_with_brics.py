import re
import collections
from rdkit import Chem
from rdkit.Chem.BRICS import BRICSDecompose
from rdkit.Chem import Draw
#from rdkit.Chem.Draw import IPythonConsole

drugbank_input = Chem.SDMolSupplier('drugbank_approved_experimental_investigational_structures.sdf')
drugbank = [m for m in drugbank_input if m]
#drugbank_scaffolds_list = [list(BRICSDecompose(mol,returnMols=True)) for mol in drugbank]
n_mols = len(drugbank)
drugbank_scaffolds_list = []
for i in range(n_mols):
	print('processing ' + str(i) + ' molecule')
	Draw.MolToFile(drugbank[i], './imgs/' + 'drug' + str(i) + '.png', size=(500, 500), fitImage=True)
	try:
		drugbank_scaffolds_list.append(list(BRICSDecompose(drugbank[i],minFragmentSize=3, returnMols=True)))
	except:
		continue

scaffold_smiles = [Chem.MolToSmiles(scaffold, kekuleSmiles=True) for drugbank_scaffolds in drugbank_scaffolds_list for scaffold in drugbank_scaffolds if scaffold != None]


counter=collections.Counter(scaffold_smiles)
print(counter.most_common(5))

with open('drugbank_approved_experimental_investigational_scaffolds.csv', 'a+') as f:
	f.write('fragment,counter\n')
	for smi, num in counter.items():
		f.write(re.sub('\[\d+\*\]', '[*]', smi) + ',' + str(num) + '\n')