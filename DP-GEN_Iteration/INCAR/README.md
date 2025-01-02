In the INCAR file, the charge of H3O+ or OH- is controlled by the NELECT keyword, which sets the number of electrons. All subsystems in the datasets are neutral, as a homogeneous background charge is assumed.

As shown in [VASP Wiki](https://www.vasp.at/wiki/index.php/NELECT), if the number of electrons is not compatible with the number derived from the valence and the number of atoms a homogeneous background charge is assumed. If the number of ions specified in the POSCAR file is 0 and NELECT=n, then the energy of a homogeneous electron gas is calculated.
