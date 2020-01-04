import openbabel as ob
import numpy as np

class ChiralIsomorphism():
    """ Disambiguate chiral centers from graph isomorphism
        Based on OpenBabel"""

    def __call__(self, obmol1, obmol2):
        """ this is the main function inside this class """

        # OpenBabel graph isomorphism
        isomorphs, sympairs_atob, sympairs_btoa = self._graph_isomorphism(obmol1, obmol2)

        # groups of atoms deemed equivalent (i.e. symmetric) by openbabel
        symgroups = self._remove_duplicates(sympairs_btoa)        

        # common atom, and filtered groups of symmetric atoms
        commons, symgroups_filt = self._get_groups_sharing_neighbor(obmol1, symgroups)

        # identify common atoms with an additional two atoms to calculate dihedrals
        # 'hooks' is a list of length-3 tuples, in which the first and second values
        # are the indeces of the two aditional atoms, and the third value is the
        # index of the common atom.
        # 'symgroups3' is the list of groups of equivalent atoms
        # a dihedral is formed by the three atoms in the hooks, plus each of
        # the atoms in the corresponding group of symgroups3
        hooks, symgroups3 = self._get_dihedrals(obmol1, commons, symgroups_filt)

        # iterate over all mappings in 'isomorphs', compute cost function
        # and return the mapping with lowest cost
        mapping, idx = self.find_best_mapping(isomorphs, hooks, symgroups3, obmol1, obmol2)

        return mapping, idx

    def _idx(self, mapping):
        """ openbabel outputs pairs of indexes.
            convert that into a good old array of indeces, for example:
            [(0, 2), (1, 0), (2, 1)] -> [2, 0, 1]
        """

        idx = [None] * len(mapping)
        for i, j in mapping:
            idx[i] = j
        return idx

    def _get_identity_idx(self, n):
        idx = []
        for i in range(n):
            idx.append(i)
        return idx

    def _angle_rmsd(self, array1, array2):
        square_diff = 0.0
        for a, b in zip(array1, array2):
            diff = self.angle_subtraction(a, b)
            square_diff += diff ** 2
        return square_diff ** 0.5


    def find_best_mapping(self, isomorphs, hooks, symgroups3, obmol1, obmol2):

        # get reference dihedrals
        n_atoms = len(isomorphs[0])
        identity_idx = self._get_identity_idx(n_atoms)
        ref_dih_diffs = self._calc_dihedral_diffs(hooks, symgroups3, obmol1, identity_idx)

        # calculate  
        lowest_rmsd = float('inf')
        i = 0
        for mapping in isomorphs:
            idx = self._idx(mapping)
            dih_diffs = self._calc_dihedral_diffs(hooks, symgroups3, obmol2, idx)
            rmsd = self._angle_rmsd(ref_dih_diffs, dih_diffs)
            if rmsd < lowest_rmsd:
                best_mapping = mapping
                best_idx = idx
                lowest_rmsd = rmsd
            i += 1
        return best_mapping, best_idx


    def _calc_dihedral_diffs(self, hooks, symgroups3, obmol, idx): 

        dihedral_diffs = []

        for h, g in zip(hooks, symgroups3):

            a = self._get_coords(obmol, idx[h[0]])
            b = self._get_coords(obmol, idx[h[1]])
            c = self._get_coords(obmol, idx[h[2]])

            # calculate angle diff between all possible pairs in symgroup 'g'
            n = len(g)
            for i in range(n):
                d_i = self._get_coords(obmol, idx[g[i]])
                dih_a = self._dihedral(a, b, c, d_i)
                for j in range(i + 1, n):
                    d_j = self._get_coords(obmol, idx[g[j]])
                    dih_b = self._dihedral(a, b, c, d_j)
                    diff = self.angle_subtraction(dih_a, dih_b)
                    dihedral_diffs.append(diff)

        return dihedral_diffs

    def _get_coords(self, obmol, i):
        atom = obmol.GetAtom(i + 1)
        x = atom.GetX()
        y = atom.GetY()
        z = atom.GetZ()
        return np.array([x, y, z])


    def angle_subtraction(self, a, b):
        """ Positive, if a > b
            Negative, if a < b
            The shortest arc between 'a' and 'b' is considered.
            For example, if 'a' = 355 deg, and 'b' = 10 deg,
            this function returns -15.
            Use RADIANS!
        """

        a = a % (2.0 * np.pi)
        b = b % (2.0 * np.pi)     
        shortest_arc = float('inf')
        for x in [-1, 0, 1]:
            tmp = a + x * 2.0 * np.pi
            arc_length = abs(tmp - b)
            if arc_length < shortest_arc:
                shortest_arc = arc_length
                out = tmp - b
        return out


    def _remove_duplicates(self, list_of_sets):
        """ return list of unique sets """

        y = set()
        for set_ in list_of_sets:
            list_ = list(set_)
            list_.sort()
            t = tuple(list_)
            y.add(t)
        return list(y)

        
    def _graph_isomorphism(self, obmol1, obmol2):
        """ Graph isomorphism from openbabel.

            "sympairs" is similar to rdkit's output,
            in which symmetric atoms get the same index.
        """

        # fundamental openbabel graph isomorphism
        query = ob.CompileMoleculeQuery(obmol1)
        mapper = ob.OBIsomorphismMapper.GetInstance(query)
        isomorphs = ob.vvpairUIntUInt()
        mapper.MapAll(obmol2, isomorphs)

        # sympairs 
        n = obmol1.NumAtoms()
        sympairs_atob = [set() for _ in range(n)]
        sympairs_btoa = [set() for _ in range(n)]
        for isomorph in isomorphs:
            for pair in isomorph:
                i = pair[0]
                j = pair[1]
                sympairs_atob[i].add(j)    
                sympairs_btoa[j].add(i)    

        return isomorphs, sympairs_atob, sympairs_btoa


    def _get_dihedrals(self, obmol, commons, symgroups_filt):
        """ this is the continuation of _get_groups_sharing_neighbor
            where a chain of two extra atoms 'hook' is identified to be
            subsequently used to calculate dihedral angles
        """

        hooks = []      # indeces of the 2 atoms bonded to 'common' that are
                        # not part of the symgroup, plus the index of 'common' 
        symgroups3 = [] # same as 'symgroups_filt', but likely to have repetitions
                        # because there may be more than 1 hook for a symgroup

        # 'c' is the index of the common atom
        # 'g' is a tuple of indeces of atoms bound to the common atom
        for c, g in zip(commons, symgroups_filt):

            # get indeces of neighbors of c
            c_atom = obmol.GetAtom(c + 1)
            neigh_of_c = self._get_neighbors(c_atom)

            for n in neigh_of_c:
                if n not in g:
                    neigh_atom = obmol.GetAtom(n + 1)
                    neigh_of_neigh = self._get_neighbors(neigh_atom) 
                    for k in neigh_of_neigh:
                        if k != c:
                            hooks.append((k, n, c))
                            symgroups3.append(g)

        return hooks, symgroups3



    def _get_groups_sharing_neighbor(self, obmol, symgroups):
        """ atom_idxs is a list of graph-symmetric atoms
            we will find groups of atoms that share a neighbor in common
            returns dictionary where key is the index of common atom
            and values are the indexes of the atoms sharing the common atom 
            groups can be repeated for different keys (e.g. cyclobutane)
            the same key can be used for multiple groups, hence no dicts
        """

        # list of atom indexes for common neighbor of a set of symmetric atoms
        commons = []

        # the corresponding list of lists of symmetric atoms
        symgroups2 = []
        

        for atom_idxs in symgroups:

            groups = {}
            
            # get neigbors for every atom
            neighbor_matrix = {}
            for i in atom_idxs:
                atom = obmol.GetAtom(i + 1)
                neighbors = self._get_neighbors(atom) 
                neighbor_matrix[i] = neighbors

            # form groups
            groups = {}
            for i in atom_idxs:
                for j in atom_idxs:
                    if i == j:
                        continue
                    # iterate over neighbors of 'i'
                    for k in neighbor_matrix[i]:
                        # if 'k' is also neighbor of 'j', group 'i' and 'j'
                        if k in neighbor_matrix[j]:
                            groups.setdefault(k, set()) 
                            groups[k].add(i)
                            groups[k].add(j)

            for common in groups:
                commons.append(common)
                symgroups2.append(tuple(groups[common]))

        return commons, symgroups2


    def _get_neighbors(self, obatom):
        """ trick or treat? """

        neighbors = []
        for bond in ob.OBAtomBondIter(obatom):
            end = bond.GetEndAtom()
            bgn = bond.GetBeginAtom()
            if bgn == obatom:
                neighbors.append(end.GetIndex()) # Index is 0-indexed
            elif end == obatom:
                neighbors.append(bgn.GetIndex())
            else:
                raise RuntimeError
        return neighbors

    def _dihedral(self, A,B,C,D):
        """Calculate dihedral considering A in the beggining"""
        A, B, C, D = [np.array(x) for x in (A,B,C,D)]
        b1 = B - A
        b2 = C - B
        b3 = D - C
        temp = np.linalg.norm(b2) * b1
        y = np.dot(temp, np.cross(b2, b3))
        x = np.dot(np.cross(b1, b2), np.cross(b2, b3))
        return np.arctan2(y, x)



