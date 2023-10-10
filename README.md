# `just-psf`
*Just get me a topology, for god’s sake!*

## Purpose

This package gives a quick and dirty way to get [PDB](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html), [RTF](https://www.charmm-gui.org/?doc=lecture&module=molecules_and_topology&lesson=2) or [PSF](https://www.charmm-gui.org/?doc=lecture&module=pdb&lesson=6) files (used by CHARMM, NAMD and VMD) from a XYZ geometry.
In particular,

1. It uses the [distance matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html) between atoms to infer bonds (using the covalent radii from [10.1039/B801115J](https://dx.doi.org/10.1039/B801115J)).
2. then it extracts the [connected components](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.connected_components.html) of the corresponding molecular graph to create resdidues.
3. If the topologies of two residues are similar  (i.e., their molecular graphs are [isomorphic](https://networkx.org/documentation/stable/reference/algorithms/isomorphism.html)), they are considered as the same residue.
4. If needed, angles and dihedrals are also inferred from bonds obtained in step 1.

Thus, everything depends on step 1: if some bonds are too large for the algorithm to detect them (e.g., metalic bonds), this will lead to incorrect results.

**Note:** for developers, this packages contains parsers that more or less faithfully extract data from said files, see [`just_psf.parers`](just_psf/parsers).

## Install

To install this package, you need a running Python 3 installation (Python >= 3.10 recommended), and

```bash
pip3 install git+https://github.com/pierre-24/just-psf.git
```

Note: as this script install programs, you might need to add their location (such as `$HOME/.local/bin`, if you use `--user`) to your `$PATH`, if any.

## Usage

You can directly obtain a PSF using `just-psf`:

```bash
just-psf tests/tests_files/7H2O.xyz -o 7H2O.psf 
```

In [this example](tests/tests_files/7H2O.xyz), all waters molecules are recognized as one residue, called `RES1`. 
Atom names are assigned asequentially in said residue (so: `O1`, `H2` and `H3`).
You may want to do find/replace, but don't forget that the PSF format is a column-based format, so **keep the alignment** when doing so.

If you prefer, you can also use [`psfgen`](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/) to build your PSF file.
For that, you need a PDB:

```bash
just-pdb tests/tests_files/7H2O.xyz -o 7H2O.pdb
```

... And a topology (also referred to as RTF, RTop, or toppar file):

```bash
just-rtf tests/tests_files/7H2O.xyz -o 7H2O.rtf
```

In this RTF file, you can more easily control things like residue name, atom names, and atom types. 
Note that in this example, it only contains one residue.

When you are happy, use [`psfgen`](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/) available in VMD:

```
# See the manual at https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf
# load topology
topology 7H2O.rtf

# load PDB
segment X { pdb 7H2O.pdb }

# write psf
writepsf 7H2O_psfgen.psf

# if you want, you can generate a more complete PDB:
coordpdb 7H2O.pdb X
writepdb 7H2O_psfgen.pdb
```

You can notice that `7H2O_psfgen.psf` and `7H2O_psfgen.pdb` are pretty similar to their `just-*` counterpart.
One main difference is that `psfgen` changes the order of the atom to match the one found in the topology (which can be an issue if you try to analyse a trajectory *a posteriori*).

## Who?

My name is [Pierre Beaujean](https://pierrebeaujean.net), and I have a Ph.D. in quantum chemistry from the [University of Namur](https://unamur.be) (Belgium).
I'm the main (and only) developer of this project, used in our lab.
I use AIMD in the frame of my post-doctoral research in order to study batteries and solid electrolyte interphrase, and I developed this project to ease my life (because to analyse an AIMD trajectory, a PSF might also be useful!).