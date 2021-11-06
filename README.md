# About

Calculate RMSD between two molecules using openbabel to assign
atoms to each other. This is relevant when the order of the atoms
differs. Hydrogens are ignored. The bond orders are set to 1 to
make the matching more forgiving. There was one case where
converting the molecule from SDF to PDB made the mapping possible.
The underlying reason was not investigated.

# Example
```sh
python sym-rmsd.py xray-lig.pdb docked-lig.pdb
```

# Dependencies

```sh
conda install 
pip install meeko
```
