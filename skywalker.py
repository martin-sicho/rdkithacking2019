"""
walker

Created by: Martin Sicho
Changed by: Christian Feldmann
On: 9/27/19, 10:43 AM
"""

from rdkit import Chem
from typing import *


class NeighbourhoodIterator:
    def __init__(self, mol: Chem.Mol, shell_count: int):
        """ Derive shells for every atom in the molecule, where the number of shells equals shell_count
        :param mol:
        :param shell_count:
        """

        self.mol: Chem.Mol = mol
        self.shell_count: int = shell_count
        self.atom_shell_atoms: dict = dict()
        self.iterate()

    def _atom_shells(self, atm: int) -> List[Set[int]]:
        """Return the Atoms in the shell
        first list represents the shells (set) containing the atom-indices (int)

        :param atm:
        :return:
        """
        shells: List[Set[int]] = []  # List of Shells (set of atomidx (int))
        current_iter_atomts: Set[int] = {atm}  # Set of atoms for current iteration, initialized with central-atom
        prev_atms = {atm}  # Atoms already in inner shells, used to avoid occurrence in multiple shells
        for i in range(self.shell_count):  # type: Chem.Atom
            next_iter_atoms: Set[int] = set()
            # Add neighbours of atoms in previous shell (potential candidates for next shell)
            for j_atm in current_iter_atomts:  # type: int
                j_atm_obj: Chem.Atom = self.mol.GetAtomWithIdx(j_atm)
                next_iter_atoms.update([k_atm.GetIdx() for k_atm in j_atm_obj.GetNeighbors()])
            # Add atoms as shell if not in one of the previous shells
            shells.append(next_iter_atoms - prev_atms)
            # Update for next loop
            current_iter_atomts = next_iter_atoms - prev_atms
            prev_atms.update(next_iter_atoms)
        assert len(shells) == self.shell_count
        return shells

    def iterate(self):
        for atm in self.mol.GetAtoms():  #type: Chem.Atom
            self.atom_shell_atoms[atm.GetIdx()] = self._atom_shells(atm.GetIdx())


def check_toluene():
    smiles = "Cc1ccccc1"
    mol = Chem.MolFromSmiles(smiles)
    iterator = NeighbourhoodIterator(mol, 5)
    iterator.iterate()
    expected = {0: [{1}, {2, 6}, {3, 5}, {4}, set()],
                1: [{0, 2, 6}, {3, 5}, {4}, set(), set()],
                2: [{1, 3}, {0, 4, 6}, {5}, set(), set()],
                3: [{2, 4}, {1, 5}, {6, 0}, set(), set()],
                4: [{3, 5}, {2, 6}, {1}, {0}, set()],
                5: [{4, 6}, {1, 3}, {2, 0}, set(), set()],
                6: [{1, 5}, {0, 2, 4}, {3}, set(), set()],
                }

    assert iterator.atom_shell_atoms == expected
    print("Passed check on toluene")


def check_ethanol():
    smiles = "CCO"
    mol = Chem.MolFromSmiles(smiles)
    iterator = NeighbourhoodIterator(mol, 3)
    expected = {0: [{1}, {2}, set()],
                1: [{0, 2}, set(), set()],
                2: [{1}, {0}, set()],
                }
    assert iterator.atom_shell_atoms == expected
    print("Passed check on ethanol")


def run_checks():
    check_toluene()
    check_ethanol()


def main():
    smiles = "Cc1ccccc1"
    mol = Chem.MolFromSmiles(smiles)
    iterator = NeighbourhoodIterator(mol, 5)
    iterator.iterate()
    print(iterator.atom_shell_atoms[1])
    print(len(iterator.atom_shell_atoms[0]))


if __name__ == "__main__":
    main()
    run_checks()
