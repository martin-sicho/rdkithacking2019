"""
tests

Created by: Martin Sicho, Christian Feldmann
On: 10/3/19, 4:51 PM
"""

import unittest
from rdkit.Chem import AllChem as Chem
from descriptors.average import AverageShell
from descriptors.simple import Simple

class TestStringMethods(unittest.TestCase):

    def test_simple(self):
        """
        Just tests if a simple descriptor calculator is able to
        find shells properly.

        :return:
        """

        smiles = "Cc1ccccc1"
        mol = Chem.MolFromSmiles(smiles)
        descriptor = Simple(mol, 5)
        descriptor.calculate()
        expected = {0: [{1}, {2, 6}, {3, 5}, {4}, set()],
                    1: [{0, 2, 6}, {3, 5}, {4}, set(), set()],
                    2: [{1, 3}, {0, 4, 6}, {5}, set(), set()],
                    3: [{2, 4}, {1, 5}, {6, 0}, set(), set()],
                    4: [{3, 5}, {2, 6}, {1}, {0}, set()],
                    5: [{4, 6}, {1, 3}, {2, 0}, set(), set()],
                    6: [{1, 5}, {0, 2, 4}, {3}, set(), set()],
                    }

        assert descriptor.atom_shell_atoms == expected
        print("Passed check on toluene")

        smiles = "CCO"
        mol = Chem.MolFromSmiles(smiles)
        descriptor = Simple(mol, 3)
        descriptor.calculate()
        expected = {0: [{1}, {2}, set()],
                    1: [{0, 2}, set(), set()],
                    2: [{1}, {0}, set()],
                    }
        assert descriptor.atom_shell_atoms == expected
        print("Passed check on ethanol")

    def test_average(self):
        #imatinib
        smiles = 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'
        mol = Chem.MolFromSmiles(smiles)
        descriptor = AverageShell(mol, 5, True)
        descriptor.calculate()
        print(descriptor.atomDescriptors)

if __name__ == '__main__':
    unittest.main()