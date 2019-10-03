"""
simple

Created by: Martin Sicho
On: 9/27/19, 3:25 PM
"""

from rdkit import Chem
from shelldesc._base import ShellDescriptorBase

class Simple(ShellDescriptorBase):
    """
    This is an example of how the descriptor interface could look like.

    """

    def __init__(self, mol: Chem.Mol, shell_count: int, include_Hs: bool = False):
        """
        Initialize the instance (say what molecule
        and how many shells to use).

        :param mol:
        :param shell_count:
        """

        super().__init__(mol, shell_count, include_Hs)
        self.atom_shell_atoms: dict = dict()

    def _processShells(self, atm_idx, shells):
        """
        Process shells for a certain atom index.
        This is where the aggregation/averaging/sorting code
        will go to. Right now we just store the shells.

        :param atm_idx:
        :param shells:
        :return:
        """

        self.atom_shell_atoms[atm_idx] = shells
