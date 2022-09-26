import os
import numpy as np
from scipy.linalg import eigh
from collections import OrderedDict
from ase.atoms import Atoms
from collections import defaultdict
from netCDF4 import Dataset
import sisl


def Lowdin(S):
    """
    Calculate S^(-1/2).
    Which is used in lowind's symmetric orthonormalization.
    psi_prime = S^(-1/2) psi
    """
    eigval, eigvec = eigh(S)
    return eigvec @ np.diag(np.sqrt(1.0/eigval)) @ (eigvec.T.conj())


def symbol_number(symbols):
    """
    symbols can be also atoms. Thus the chemical symbols will be used.
    Fe Fe Fe O -> {Fe1:0 Fe2:1 Fe3:2 O1:3}
    """
    try:
        symbs = symbols.copy().get_chemical_symbols()
    except Exception:
        symbs = symbols
    symdict = {}
    result = OrderedDict()
    for i, sym in enumerate(symbs):
        if sym in symdict:
            symdict[sym] = symdict[sym]+1
        else:
            symdict[sym] = 1
        result[sym+str(symdict[sym])] = i
    return result


class SislWrapper():

    def __init__(self, sisl_geom, spin, ispin=None, shift_fermi=None):
        self.shift_fermi = shift_fermi
        self.spin = spin
        self.ispin = ispin

        self.orbs = []
        self.orb_dict = defaultdict(lambda: [])
        g = sisl_geom
        _atoms = g._atoms
        atomic_numbers = []
        atom_positions = g.xyz
        self.cell = np.array(g.sc.cell)
        for ia, a in enumerate(_atoms):
            atomic_numbers.append(a.Z)
        self.atoms = Atoms(numbers=atomic_numbers,
                           cell=self.cell,
                           positions=atom_positions)
        xred = self.atoms.get_scaled_positions()
        sdict = list(symbol_number(self.atoms).keys())

        if self.spin.is_colinear:
            if ispin is None:
                raise ValueError("For colinear spin, ispin must be given")
        else:
            if ispin is not None:
                raise ValueError(
                    "For non-colinear spin and unpolarized spin, ispin should be None"
                )

        self.positions = []
        if self.spin.is_colinear:
            for ia, a in enumerate(_atoms):
                symnum = sdict[ia]
                orb_names = []
                for x in a.orbitals:
                    name = f"{symnum}|{x.name()}|{ispin}"
                    orb_names.append(name)
                    self.positions.append(xred[ia])
                self.orbs += orb_names
                self.orb_dict[ia] += orb_names
            self.norb = len(self.orbs)
            self.nbasis = self.norb
        elif self.spin.is_spinorbit:
            for ispin in ['up', 'down']:
                for ia, a in enumerate(_atoms):
                    symnum = sdict[ia]
                    orb_names = []
                    for x in a.orbitals:
                        name = f"{symnum}|{x.name()}|{spin}"
                        orb_names.append(name)
                        self.positions.append(xred[ia])
                    self.orbs += orb_names
                    self.orb_dict[ia] += orb_names
            self.norb = len(self.orbs) / 2
            self.nbasis = len(self.orbs)
        else:
            for ia, a in enumerate(_atoms):
                symnum = sdict[ia]
                orb_names = []
                for x in a.orbitals:
                    name = f"{symnum}|{x.name()}|None"
                    orb_names.append(name)
                    self.positions.append(xred[ia])
                self.orbs += orb_names
                self.orb_dict[ia] += orb_names
            self.norb = len(self.orbs)
            self.nbasis = len(self.orbs)
        self.positions = np.array(self.positions, dtype=float)

    def print_orbs(self):
        print(self.orb_dict)

    def solve_all(self, kpts):
        """ Get eigenvalues and eigenvectors for all kpts. Should be
        implemented in children."""
        pass

    def get_fermi_level(self):
        if self.shift_fermi:
            return self.shift_fermi
        else:
            return 0.0


class SislHSWrapper(SislWrapper):
    def __init__(self,
                 sisl_hamiltonian,
                 ispin=None,
                 shift_fermi=None,
                 format='dense',
                 nbands=10):

        self.format = format
        self.nbands = nbands
        self.ham = sisl_hamiltonian
        super().__init__(self.ham.geometry, self.ham.spin, ispin, shift_fermi)

    def solve(self, k):
        if self.format == 'sparse':
            # Calculates a subset of eigenvalues (nbands) using the sparse
            # algorithms.
            if self.ispin is None:
                evals, evecs = self.ham.eigsh(
                    k=k, n=self.nbands, eigvals_only=False, which='SA', mode='normal', sigma=10)
            else:
                evals, evecs = self.ham.eigsh(
                    k=k, n=self.nbands, spin=self.ispin, eigvals_only=False,  which='SA', mode='normal', sigma=1)
        else:
            if self.ispin is None:
                evals, evecs = self.ham.eigh(k=k, eigvals_only=False)
            else:
                evals, evecs = self.ham.eigh(
                    k=k, spin=self.ispin, eigvals_only=False)
        if self.shift_fermi:
            evals += self.shift_fermi
        return evals, evecs

    def Hk(self, k, format='dense'):
        if self.ispin is not None:
            return self.ham.Hk(k, format=format, spin=self.ispin)
        else:
            return self.ham.Hk(k, format=format)

    def solve_all(self, kpts, orth=False):
        evals = []
        evecs = []
        for ik, k in enumerate(kpts):
            if orth and self.ham.orthogonal:
                if self.format == 'sparse':
                    raise NotImplementedError('Currently format="sparse" is not'
                                              ' compatible with orth="True", because there is no'
                                              ' sparse matrix square root implemented.')
                S = self.ham.Sk(k, format='dense')
                Smh = Lowdin(S)
                H = self.Hk(k, format='dense')
                Horth = Smh.T.conj() @ H @ Smh
                evalue, evec = eigh(Horth)
            else:
                evalue, evec = self.solve(k)
            evals.append(evalue)
            evecs.append(evec)
        return np.array(evals, dtype=float), np.array(evecs,
                                                      dtype=complex,
                                                      order='C')


class SislWFSXWrapper(SislWrapper):
    """ Wrapper for retrieving eigenvalues and eigenvectors from siesta WFSX file

    Parameters
    ----------
    geom : sisl.Geometry
        the geometry containing information about atomic positions and orbitals
    wfsx_sile: sisl.io.siesta.wfsxSileSiesta
        file reader for WFSX file
    spin : sisl.physics.Spin 
        spin object carrying information about spin configurations and spin component.
    ispin : None or int
        index of spin channel to be considered. Only takes effect for collinear spin calculations (UP: 0, DOWN: 1). 
        (default: None)
    shift_fermi: None or float
        energy shift to be applied to all energies. If `None` no shift is applied. (default: None)
    """

    def __init__(self, geom, wfsx_sile, spin, ispin=None, shift_fermi=None):
        super().__init__(geom, spin=spin, ispin=ispin, shift_fermi=shift_fermi)
        self.geom = geom
        self.wfsx_sile = wfsx_sile
        self.read_all()

    def read_all(self):
        """ Read all wavefunctions, eigenvalues, and k-points from WFSX file."""
        evals = []
        evecs = []
        wfsx_kpts = []
        wfsx_weights = []

        def change_gauge(k, evec):
            """ Change the eigenvector gauge """
            phase = np.dot(
                self.geom.xyz[self.geom.o2a(np.arange(self.geom.no)), :],
                np.dot(k, self.geom.rcell))
            if self.spin.has_noncolinear:
                # for NC/SOC we have a 2x2 spin-box per orbital
                phase = np.repeat(phase, 2)
            # r -> R gauge tranformation.
            return evec * np.exp(1j * phase).reshape(1, -1)

        # Read wavefunctions and eigenvalues
        for wfc in self.wfsx_sile.yield_eigenstate():
            wfsx_kpts.append(wfc.info['k'])
            wfsx_weights.append(wfc.info['weight'])
            evals.append(wfc.c)
            # To get the same eigvecs as eigh returns we need to transpose the
            # array and change the gauge from 'r' to 'R'
            evecs.append(change_gauge(wfc.info['k'], wfc.state).T)

        wfsx_kpts = np.asarray(wfsx_kpts)
        wfsx_weights = np.asarray(wfsx_weights)

        # If any k-point occurs more than once in the WaveFuncKPoints block,
        # we discard the duplicates
        is_duplicate = self._is_duplicate(wfsx_kpts)
        self.wfsx_kpts = wfsx_kpts[~is_duplicate]
        self.wfsx_weights = wfsx_weights[~is_duplicate]
        # Normalize weights
        self.wfsx_weights /= np.sum(self.wfsx_weights)
        self.evals = np.array(evals, dtype=float)[~is_duplicate]
        if self.shift_fermi is not None:
            self.evals += self.shift_fermi
        self.evecs = np.array(evecs, dtype=np.complex64,
                              order='C')[~is_duplicate]

    def write_to_netcdf(self, fname):
        root = Dataset(fname, 'w', data_model='NETCDF4')
        nkpt, three = self.wfsx_kpts.shape
        print(f"nkpt: {nkpt}")
        print(f"three: {three}")
        nkpt, nband = self.evals.shape
        print(f"nkpt: {nkpt}")
        print(f"nband: {nband}")

        print(f"evecs.shape: {self.evecs.shape}")
        nkpt, nbasis, nband = self.evecs.shape
        natom = len(self.atoms)

        id_nkpt = root.createDimension("nkpt", size=nkpt)
        id_three = root.createDimension("three", three)
        id_nband = root.createDimension("nband", nband)
        id_norb = root.createDimension("norb", self.norb)
        id_natom = root.createDimension("natom", natom)
        id_nbasis = root.createDimension("nbasis", nbasis)
        # only soc.
        id_spintype = root.createDimension("spintype", 1)
        root.createVariable("cell", datatype=float,
                            dimensions=("three", "three"))
        root.createVariable("xred", datatype=float,
                            dimensions=("natom", "three"))
        root.createVariable("atomic_numbers", datatype=int,
                            dimensions=("natom"))
        root.createVariable("spintype", datatype=int)
        root.createVariable("kpts", datatype=float,
                            dimensions=("nkpt", "three"))
        root.createVariable("kweights", datatype=float,
                            dimensions=("nkpt"))

        root.createVariable("eigenvalues", datatype=float,
                            dimensions=("nkpt", "nband"))
        root.createVariable("eigenvectors_real", datatype=float,
                            dimensions=("nkpt", "nbasis", "nband"))
        root.createVariable("eigenvectors_imag", datatype=float,
                            dimensions=("nkpt", "nbasis", "nband"))
        #root.createVariable("orb_names",  datatype=str, dimensions=("nbasis"))

        root.variables["cell"][:] = self.atoms.get_cell().array
        root.variables["xred"][:] = self.atoms.get_scaled_positions()
        root.variables["atomic_numbers"][:] = self.atoms.get_atomic_numbers()
        root.variables["eigenvalues"][:] = self.evals
        root.variables["kpts"][:] = self.wfsx_kpts
        root.variables["kweights"][:] =  self.wfsx_weights
        # print(self.orbs)
        # print(len(self.orbs))
        #root.variables["orb_names"][:] = self.orbs
        root.variables["eigenvectors_real"][:] = self.evecs.real
        root.variables["eigenvectors_imag"][:] = self.evecs.imag
        root.close()

    def _is_duplicate(self, array):
        # TODO: Move into utils
        # Find all matches (i,j): array[i] == array[j]
        matches = np.all(np.isclose(array[None, :], array[:, None]), axis=-1)
        # Remove double counting of matches: (i,j) and (j,i)
        # Also, remove matches of elements with themselves: (i,i)
        matches = np.triu(matches, 1)

        # Finally determine elements which are duplicates
        return np.any(matches, axis=0)

    def find_all(self, kpts):
        """ Find the correct eigenvectors and eigenvalues and sort them
        to match the order in kpts.

        Parameters
        ----------
        kpts : list of float (3,) or (nk, 3)
            list of k points

        Returns
        -------
        evals : list of float (nk, n)
            list of eigenvalues for every requested k point
        evecs :
            list of eiegnvector for every requested k point
        """
        kpts = np.asarray(kpts)
        sort_idx = np.where(
            np.all(np.isclose(self.wfsx_kpts[None, :], kpts[:, None]),
                   axis=-1))[1]
        if len(sort_idx) < len(kpts):
            # k-point not found
            raise ValueError(
                f"{self.__class__.__name__}._read_all unable to find at least one "
                "required k point in '{self.wfsx_sile.file}'. Please, ensure that "
                "all k points listed below are included:\n" +
                "\n".join([str(k) for k in kpts]))
        if not np.all(np.isclose(self.wfsx_kpts[sort_idx], kpts)):
            raise ValueError(
                f"{self.__class__.__name__}._read_all was unable to match k points "
                "in {self.wfsx_sile.file} against required kpts:" +
                "\n".join([str(k) for k in kpts]))

        return self.evals[sort_idx], self.evecs[sort_idx]

    # Define alias for find_all to unify in wrapper's interface
    solve_all = find_all

    def set_kweights(self, kpts, kweights):
        
        kpts = np.asarray(kpts)
        sort_idx = np.where(
            np.all(np.isclose(self.wfsx_kpts[None, :], kpts[:, None]),
                   axis=-1))[1]
        print(self.wfsx_weights)
        print(self.wfsx_weights[sort_idx])
        self.wfsx_weights[sort_idx] = kweights


def write_to_netcdf(folder=None, spin="soc", fdf_fname="siesta.fdf",
                    wfxfname="siesta.fullBZ.WFSX", kwfname=None,
                    ncfname="wf.nc"):
    if spin != "soc":
        raise NotImplementedError("currently only spin=soc is implemented.")
    from sisl.physics.spin import Spin
    fdf = sisl.get_sile(os.path.join(folder, fdf_fname))
    geom = fdf.read_geometry()
    wfsx_sile = sisl.io.siesta.wfsxSileSiesta(
        wfxfname, parent=geom)
    try:
        efermi = fdf.read_fermi_level().data[0]
    except:
        efermi = fdf.read_fermi_level()
    model = SislWFSXWrapper(geom, wfsx_sile, spin=Spin(spin),
                            ispin=None, shift_fermi=None)

    if kwfname is not None:
        tmp = np.genfromtxt(kwfname)
        kpts = tmp[:,:3]
        kweights = tmp[:,3]
        model.set_kweights(kpts, kweights)

    model.write_to_netcdf(ncfname)


write_to_netcdf(folder='../siesta', fdf_fname='in.fdf',
                wfxfname="../siesta/graphene.selected.WFSX",
                kwfname='kpoints',
                ncfname='wfk.nc')
