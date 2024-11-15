from ase.calculators.abacus import Abacus
from ase.calculators.abacus import AbacusProfile
from ase.io import read, write

def abacus_init(fname, profile, **kwargs):
    '''instantiate a instance of Abacus calculator (and bind) with 
    specific structures and ABACUS parameters
    
    Parameters
    ----------
    fname : str
        The filename of the structure, can be any format that ASE supports.
        For a complete list of supported formats, see:
        https://wiki.fysik.dtu.dk/ase/ase/io/io.html#module-ase.io
    profile : AbacusProfile
        The command to run ABACUS
    **kwargs : dict
        The ABACUS parameters, e.g. pp, basis, etc.
    
    Returns
    -------
    a : Atoms
        The Atoms object with the attached ABACUS calculator
    '''
    a = read(fname)
    pp = kwargs.pop('pp', {})
    # if there is no data in pp, pops an warning
    if not pp:
        print('WARNING: No pseudopotential data is provided')
    a.calc = Abacus(pp=pp, profile=profile, **kwargs)
    return a

def relax(a, method='BFGS', fmax=0.05, max_steps=100, **kwargs):
    '''relax the structure using the given method'''
    if method == 'BFGS':
        from ase.optimize import BFGS
        dyn = BFGS(a, trajectory=None, **kwargs)
        dyn.run(fmax=fmax, steps=max_steps)
    else:
        raise NotImplementedError(f'{method} is not supported by ASE')
    return a

if __name__ == '__main__':
    pp={'Si': 'Si.ccECP.upf', 'H': 'H.ccECP.upf'}
    env = AbacusProfile('mpirun -np 16 /usr/local/bin/abacus')
    recipe = {'suffix': 'abacus',
              'pseudo_dir': '.',
              'orbital_dir': '.',
              'nspin': 2,
              'nupdown': 2,
              'symmetry': 0,
              'dft_functional': 'PBE',
              'ecutwfc': 100,
              'scf_thr': 1e-8,
              'scf_nmax': 50,
              'basis_type': 'pw',
              'gamma_only': True,
              'smearing_method': 'fixed',
              'mixing_type': 'broyden',
              'mixing_beta': 0.4,
              'mixing_beta_mag': 0.4,
              'mixing_gg0': 0,
              'mixing_gg0_mag': 0,
              'ks_solver': 'dav_subspace',
              'out_chg': -1}
    ini = 'SiH2-pos-ini.cif'
    fin = 'SiH2-pos-fin.cif'
    a = abacus_init(ini, profile=env, pp=pp, **recipe)
    a = relax(a, fmax=0.01, max_steps=100)
    write(fin, a)
    print(f'Final structure is saved in {fin}')