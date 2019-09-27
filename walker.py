"""
walker

Created by: Martin Sicho
On: 9/27/19, 10:43 AM
"""

from rdkit import Chem

class NeighboroodIterator:

    def __init__(self, mol, shell_count):
        self.mol = mol
        self.shell_count = shell_count

    def iterate(self):
        for atm in self.mol.GetAtoms():
            atm_index = atm.GetIdx()
            neighborhood_map = dict()
            atoms_current = set()
            atoms_next = set()
            atoms_current.add(atm)
            neighborhood_map[atm_index] = set(atoms_current)
            # TODO: collector code for the first atom here

            for i in range(self.shell_count): # this just iterates over the graph, but we could also take into account real distance in 3D and stuff
                atoms_next.clear()
                for atm_current in atoms_current:
                    neighbors = atm_current.GetNeighbors()
                    atoms_next.update(neighbors)
                    # remove hydrogens code here if we need it

                if i > 0:
                    atoms_next.difference_update(neighborhood_map.get(i-1))

                neighborhood_map[i+1] = set(atoms_next)
                # TODO: collection code here again
                atoms_current.clear()
                atoms_current.update(atoms_next)


def main():
    smiles = "Cc1ccccc1"
    mol = Chem.MolFromSmiles(smiles)
    iterator = NeighboroodIterator(mol, 5)
    iterator.iterate()

if __name__ == "__main__":
    main()



