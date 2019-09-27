"""
walker

Created by: Martin Sicho
On: 9/27/19, 10:43 AM
"""

from rdkit import Chem
from typing import *


class NeighboroodIterator:

    def __init__(self, mol, shell_count):
        self.mol: Chem.Mol = mol
        self.shell_count: int = shell_count
        self.atom_shell_atoms: dict = dict()

    def _atom_shells(self, atm: int) -> List[Set[int]]:
        """Return the Atoms in the shell
        first list represents the shells (list) containing the atom-indices (int)

        :param atm:
        :return:
        """
        shells: List[Set[int]] = []  # List of Shells (set of atomidx (int))
        current_iter_atomts: List[int] = [atm]  # List of atoms for current iterattion
        prev_atms = set({atm})  # Atoms already in inner shells
        for i in range(self.shell_count):  #type: Chem.Atom
            next_iter_atomts = []
            # Add neighbours for next shell
            for j_atm in current_iter_atomts:
                j_atm_obj = self.mol.GetAtomWithIdx(j_atm)
                next_iter_atomts.extend([k_atm.GetIdx() for k_atm in j_atm_obj.GetNeighbors()])
            # Add neighbours if not in one of the previous shells
            shells.append(set(next_iter_atomts) - prev_atms)
            prev_atms.update(next_iter_atomts)
            current_iter_atomts = next_iter_atomts
        return shells


    def iterate(self):
        for atm in self.mol.GetAtoms():  #type: Chem.Atom
            self.atom_shell_atoms[atm.GetIdx()] = self._atom_shells(atm.GetIdx())

def main():
    smiles = "Cc1ccccc1"
    mol = Chem.MolFromSmiles(smiles)
    iterator = NeighboroodIterator(mol, 5)
    iterator.iterate()
    print(iterator.atom_shell_atoms[0])
    print(len(iterator.atom_shell_atoms[0]))
if __name__ == "__main__":
    main()



