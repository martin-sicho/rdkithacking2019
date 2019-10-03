"""
This is the shell descriptor base class.

Created by: Martin Sicho, Christian Feldmann, Christoph Bauer
On: 10/3/19, 3:39 PM
"""

from abc import ABC, abstractmethod
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdEHTTools
from typing import *

class ShellDescriptorBase(ABC):

    def __init__(self, mol: Chem.Mol, shell_count: int, include_Hs : bool = False):
        """
        Initialize the descriptor (provide the molecule
        we want to calculate the descritptors for
        and how many shells to consider).

        :param mol: RDKit molecule to use in the calculations
        :param shell_count:
        :param include_Hs: specify whether hydrogen atoms should be considered as shell members
        """

        self.shellCount : int = shell_count
        self.includeHs = include_Hs
        self.mol: Chem.Mol = self._getMolWithEHTcharges(mol)
        self.atomDescriptors: dict = dict()
        super().__init__()

    def _getMolWithEHTcharges(self, mol : Chem.Mol):
        """
        Prepare the molecule and calculate the charges.

        TODO: perhaps this should be divided into more methods for more flexibility...

        :param mol:
        :return:
        """

        mol_ = Chem.RemoveHs(mol)
        mol_ = Chem.AddHs(mol_)
        Chem.EmbedMolecule(mol_)
        # should create two output files ('run.out' and 'nul') TODO: maybe this is not really needed or even desirable in the final version?
        passed, res = rdEHTTools.RunMol(mol_)
        nat = len(mol_.GetAtoms())
        charges = res.GetAtomicCharges()
        if not self.includeHs:
            mol_ = Chem.RemoveHs(mol)
        # creates a mol object with charges set as double properties on atoms of that mol object
        for i in range(nat):
            if i < len(mol_.GetAtoms()):
                mol_.GetAtomWithIdx(i).SetDoubleProp('EHTcharge', charges[i])
            else:
                break
        return mol_

    def _iterate(self):
        """
        Iterates over the molecular graph searching for shells of each atom.

        Written by: Christian Feldmann
        """

        for atom in self.mol.GetAtoms():  #type: Chem.Atom
            atm_idx = atom.GetIdx()
            shells: List[Set[int]] = []  # List of Shells (set of atomidx (int))
            current_iter_atomts: Set[int] = {atm_idx}  # Set of atoms for current iteration, initialized with central-atom
            prev_atms = {atm_idx}  # Atoms already in inner shells, used to avoid occurrence in multiple shells
            for i in range(self.shellCount):  # type: Chem.Atom
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
            assert len(shells) == self.shellCount
            self._processShells(atm_idx, shells)

    @abstractmethod
    def _processShells(self, atm_idx, shells):
        pass

    def calculate(self):
        self._iterate()
        return self.atomDescriptors
