# Atom Charge Shell Descriptors

This project was created on the RDKit UGM 2019 Hackathon and it comprises of a simple Python package to calculate circular atomic descriptors based on EHT charges of neighboring atoms. So far only a few very basic descriptors are calculated, but it should be easy to add more variants under the common `ShellDescriptorBase` interface.

# Example Usage

```python
from descriptors.average import AverageShell
from rdkit import Chem

#imatinib
smiles = 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'
mol = Chem.MolFromSmiles(smiles)
descriptor = AverageShell(mol, 5, True)
descriptor.calculate()
print(descriptor.atomDescriptors)
```
