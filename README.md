# About

Calculate RMSD between two molecules. Hydrogens are ignored.
Openbabel maps atoms to each other using graph isomorphism, 
which is useful when the order of the atoms differs.
The bond orders are set to 1 to make the atom mapping more forgiving.
There was one case where converting the molecule from SDF to PDB made the mapping
possible, but the underlying reason was not investigated.

# Example
```sh
python sym-rmsd.py xray-lig.pdb docked-lig.pdb
```

# Dependencies

```sh
conda install -c conda-forge openbabel
pip install meeko
```
