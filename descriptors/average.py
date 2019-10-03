"""
A simple descriptor which averages partial charges in each shell.

Created by: Martin Sicho, Christoph Bauer
On: 9/27/19, 3:33 PM
"""

from descriptors._base import ShellDescriptorBase

class AverageShell(ShellDescriptorBase):
    """
    Calculate simple charge averages in all shells.

    """

    def _processShells(self, atm_idx: int, shells : list):
        """
        This method is called each time we fetch shells for an atom.

        :param atm_idx: index of the atom for which neighboring shells are considered
        :param shells: neighborhoods as atom indices, list of sets of atoms belonging to each shell: [{1stshell}, {2ndshell}, ...]
        :return: the calculated descriptor (lenght determined by `shellCount`)
        """

        descriptor = []
        atom = self.mol.GetAtomWithIdx(atm_idx)
        descriptor.append(atom.GetDoubleProp('EHTcharge'))
        for j in range(self.shellCount):
            sum_shell = 0
            indices = list(shells[j])
            nat_in_this_shell = len(indices)
            for k in range(nat_in_this_shell):
                sum_shell += self.mol.GetAtomWithIdx(indices[k]).GetDoubleProp('EHTcharge')
            shell_average = sum_shell/nat_in_this_shell
            descriptor.append(shell_average)

        self.atomDescriptors[atm_idx] = descriptor
        return descriptor

