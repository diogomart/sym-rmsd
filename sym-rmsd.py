#!/usr/bin/env python

import argparse, sys
import pybel
import math
import os
import sys
import openbabel as ob

class OBRMSD():
    def __init__(self, ref_obmol, test_obmol, graph='isomorph'):
        self.graph      = graph
        self.refxyz     = [(a.GetX(), a.GetY(), a.GetZ()) for a in 
                            ob.OBMolAtomIter(ref_obmol)]
        self.ref = ref_obmol

        # run symmetry operations (obabel)
        if graph == 'automorph':
            obmorphs = self._get_automorphisms()
        elif graph == 'isomorph':
            obmorphs = self._get_isomorphisms(test_obmol)
        else:
            raise RuntimeError('graph must be automorph or isomorph')

        self.morphs, self.static = self._process_obmorphs(obmorphs)

        # check if elements match in ref_atoms and test_atoms
        if graph == 'automorph':
            self._check_elements_order(test_atoms)

        # Moving atoms from the reference molecule
        self.dynamicref = [coords for (coords, c) in zip(
                        self.refxyz, self.static) if not c]

    def rmsd(self, obmol):
        """ calculate sym-rmsd """

        xyz = [(a.GetX(), a.GetY(), a.GetZ()) for a in ob.OBMolAtomIter(obmol)]

        # calc squared distance for static atoms
        static_sq_dist = 0
        for i,static in enumerate(self.static):
            if static:
                static_sq_dist += self._sqdist(self.refxyz[i], xyz[self.morphs[0][i]])

        # build moving atoms list
        min_sq_dist = 9999999.9
        for morph in self.morphs:
            this_mol = [xyz[i] for (i, c) in zip(morph, self.static) if not c]
            this_sq_dist = sum(self._sqdist(a, b) for
                                    (a, b) in zip(self.dynamicref, this_mol))
            min_sq_dist = min(min_sq_dist, this_sq_dist)

        return math.sqrt((min_sq_dist + static_sq_dist) / float(len(self.refxyz)))

    def byatom(atoms_list):
        pass

    def _check_elements_order(self, test_atoms):
        """ this may spot incorrect uses of graph=automorph"""

        test_atoms = [a for a in test_atoms]# if a.element != 'H']
        for (i, refatom) in enumerate(self.ref_atoms):
            if refatom.element != test_atoms[i].element:
                raise RuntimeError('Order of elements does not match in test_atoms')

    def _sqdist(self, coordsA, coordsB):
         """Find the squared distance between two 3-tuples"""
         sqrdist = sum( (a-b)**2 for a, b in zip(coordsA, coordsB))
         return sqrdist

    def _get_isomorphisms(self, test_atoms):
    
        ref = self.ref
        mol = test_atoms
    
        # Make all bond orders 1 and reassign atom types
        # this is a magic recipe
        for bond in ob.OBMolBondIter(ref):   bond.SetBondOrder(1)
        for bond in ob.OBMolBondIter(mol):   bond.SetBondOrder(1)
        ob.OBAtomTyper().AssignTypes(ref)
        ob.OBAtomTyper().AssignTypes(mol)
    
        # DEBUG
        if False:
            rt = [a for a in ob.OBMolAtomIter(ref)]
            mt = [a for a in ob.OBMolAtomIter(mol)]
            for (r,m) in zip(rt, mt):
                print r.GetType(), m.GetType(), r.GetType()==m.GetType()
            r = self._atoms2obabel(ref_atoms, ref_type)
            m = self._atoms2obabel(mol_atoms)
            ro = pybel.Outputfile('mol2','ref.mol2')
            mo = pybel.Outputfile('mol2','mol.mol2')
            ro.write(r)
            mo.write(m)
            ro.close()
            mo.close()
    
        # Mapping magic
        query = pybel.ob.CompileMoleculeQuery(ref)
        mapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
        isomorphs = pybel.ob.vvpairUIntUInt()
        mapper.MapAll(mol, isomorphs)
        return isomorphs

    def _get_automorphisms(self):
    
        mol = self._atoms2obabel(self.ref_atoms).OBMol
        automorphs = ob.vvpairUIntUInt()
        ob.FindAutomorphisms(mol, automorphs)
        return automorphs

    def _morph2idx(self, a):
        """ Turns the isomorphism output into index changes"""
        idx = [0 for i in a]
        for i,j in a:
            idx[i] = j
        return idx

    def _process_obmorphs(self, obmorphs):
        """
        1st :: convert automorphisms to index replaces
        2nd :: find out atoms that have no symmetry mates, i.e., static atoms
        """

        # 1st :: conversion
        my_morphs = [self._morph2idx(obmorph) for obmorph in obmorphs]
    
        # 2nd :: find static atoms, no symmetry mates
        counts = [[] for (i,j) in obmorphs[0]] # stores indxs for a given position
        for morph in my_morphs:
            for (i, idx) in enumerate(morph):
                counts[i].append(idx)   # append morphisms
        static = [len(set(count)) == 1 for count in counts] # no symmetry mates!
    
        return my_morphs, static

def readmol(fname):
    extension = fname.split('.')[-1]
    obmol_generator = pybel.readfile(extension, fname)
    obmols = [mol.OBMol for mol in obmol_generator]
    for obmol in obmols:
        to_del = []
        for atom in ob.OBMolAtomIter(obmol):
            #if atom.IsNonPolarHydrogen():
            if atom.IsHydrogen():
                to_del.append(atom)
        for atom in to_del[::-1]:
            obmol.DeleteAtom(atom)
    return obmols


def main():

    # parse user input
    parser = argparse.ArgumentParser()
    parser.add_argument('ref', help='reference molecule for rmsd calculation')
    parser.add_argument('query', help='query molecule for rmsd calculation')
    args = parser.parse_args(sys.argv[1:])

    # read molecules
    refs = readmol(args.ref)
    if len(refs) != 1:
        sys.stderr.write('The reference molecule must contain 1 conformer\n')
        sys.exit(2)
    queries = readmol(args.query)

    # check if molecule is the same - _can_onical smiles
    ext_r = os.path.splitext(args.ref)[1].replace('.', '')
    ext_q = os.path.splitext(args.query)[1].replace('.', '')
    r = pybel.readfile(ext_r, args.ref).next()
    q = pybel.readfile(ext_q, args.query).next()
    r.OBMol.CorrectForPH(7.4)
    q.OBMol.CorrectForPH(7.4)
    can_ref = r.write('can').split()[0]
    can_q = q.write('can').split()[0]
    if can_ref != can_q:
        sys.stderr.write('Different canonical smiles for input molecules\n')
        sys.stderr.write('%s\n' % can_ref)
        sys.stderr.write('%s\n' % can_q)
       # print "-2"
       # sys.exit()

    # initialize rmsd calculator
    obrmsd = OBRMSD(refs[0], queries[0])

    for query in queries:
        rmsd = obrmsd.rmsd(query)
        print '%6.2f' % rmsd

if __name__ == '__main__':
    main()


