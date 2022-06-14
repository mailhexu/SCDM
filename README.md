# SCDM

## Tutorial

### Installation

​	  First, these dependencies should be already installed. 

- CMake
- MPI
- libnetcdf/libnetcdf-fortran
- BLAS
- LAPACK
- Fortran compiler 
- json-fortran (https://github.com/jacobwilliams/json-fortran

​	On debain/ubuntu linux, the dependencies except json-fortran can be installed with

```shell
sudo apt install openmpi libnetcdff-dev libopenblas-dev cmake gfortran
```

To install json-fortran, please follow this page. 

https://github.com/jacobwilliams/json-fortran

​	Then we can install SCDM with the following command:

```bash
mkdir build
cd build
cmake ..
make
make install
```

After installation, there will be a sdown command. 

#### Usage:

To use sdown to build the wannier functions from the siesta output, the following steps 

- Run siesta, and generate the wavefunctions. To save the eigenvalues and the wavefunctions, following options should be turned on :

  ```
  COOP.write True
  WriteEigenvalues True
  ```

  After the Siesta scf calculation, there should be a ".WFSX" file.

- Convert the wavefunction file into netcdf format. 

  There is a write_wfk_nc.py  script, which can convert the WFSX file into netcdf format, named wf.nc. The last line of the script should be modified to adapt to the folder and filenames. 

  ```python
  write_to_netcdf(folder='./', fdf_fname='siesta.fdf',
                    wfxfname="siesta.fullBZ.WFSX", ncfname='wf.nc')
  ```

  Note that to use this script, the python libaries including ase, netcdf4, sisl, scipy should be installed.

- Prepare the input files. sdown reads a json file named input.json as input. Here is an example:

  ````
  {                             
      "method": 1,
      "disentangle_func_type": "gauss",
      "mu": -4,
      "sigma": 5,
      "kmesh": [
          6,
          6,
          6
      ],
      "nwann": 5,
      "project_to_anchor": true,
      "anchor_kpt": [
          0,
          0,
          0
      ],
      "anchor_ibands": [
          46,
          47,
          48,
          49,
          50
      ]
  }
  ````

  - method: the method for building wannier functions.

    1. SCDM-k method
    2. projected wannier function method. (not finished yet.)

  - disentangle_func_type: the type of function for disentanglement.  unity, gauss or Fermi.  The Gauss function is centered at mu and has a width of sigma. The Fermi function with Fermi energy of mu and smearing of sigma.  

  - mu

  - sigma

  - kmesh: should be the same as in the siesta input. 

  - nwann: number of wannier functions.

  - project_to_anchor: whether use the projection to the anchor-wavefunctions for disentanglement. 

  - anchor_kpt: the anchor k-point. 

  - anchor_ibands: the index of the band used as anchor points.  It is automatically decided if not specified. The size should be nwann. 

    

  With all these files prepared, we can now run sdown to build the wannier functions.  By running the command "sdown", a file named wannier.nc will be generated, which contains the wannier functions and the Hamiltonian. 

  

- Plot the Wannier band structure. 

  IN the script directory there is an "plot_wannier_band.py" script. Run it to plot the band structures. 

  The python library "wannierbuilder" needs to be installed for this script. 

  ```
  pip install wannierbuilder
  ```

  

### Example:
An example of running sdown can be found in the example directory. 
In the SrMnO3_SOC/siesta directory, there are the siesta input file and pseudopotentials. The write_wfk_nc.py file can then be used to convert the wavefunction file into netcdf format. In the build_wannier directory, the input.json file (also a python script for generating it), and the plot_wannier_band.py file can be found.

  

​	





