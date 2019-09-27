from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdEHTTools


# imatinib as an example
m = Chem.MolFromSmiles('Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1')
mh = Chem.AddHs(m)
Chem.EmbedMolecule(mh)
# should create two output files ('run.out' and 'nul')
passed,res = rdEHTTools.RunMol(mh)
nat = len(mh.GetAtoms())
print(nat)
print(len(res.GetAtomicCharges()))
charges = res.GetAtomicCharges()
print(charges)

#creates a mol object with charges set as double properties on atoms of that mol object
for i in range(nat):
    mh.GetAtomWithIdx(i).SetDoubleProp('EHTcharge',charges[i])
