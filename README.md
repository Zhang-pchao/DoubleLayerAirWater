# Paper Title

Double-layer distribution of hydronium and hydroxide ions in the air-water interface [ACS Phys. Chem Au 2024, 4, 4, 336–346](https://pubs.acs.org/doi/10.1021/acsphyschemau.3c00076) [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/64a1f26aba3e99daef69917a)

# Dataset and model

The dataset for training the DP model is uploaded to [AIS Square](https://www.aissquare.com/datasets/detail?pageType=datasets&name=SCAN_H2O_H3O_OH&id=243) and [Zenodo](https://zenodo.org/records/14306810).

The compressed DP model is uploaded to [AIS Square](https://www.aissquare.com/models/detail?pageType=models&name=SCAN_H2O_H3O_OH&id=242) and [Zenodo](https://zenodo.org/records/14306810).

In the [INCAR file](https://github.com/Zhang-pchao/DoubleLayerAirWater/tree/main/DP-GEN_Iteration/INCAR), the charge of H₃O⁺ or OH⁻ is controlled by the NELECT keyword, which sets the number of electrons. All subsystems in the DP datasets are neutral because of a homogeneous background charge. As shown in [VASP Wiki](https://www.vasp.at/wiki/index.php/NELECT), if the number of electrons is not compatible with the number derived from the valence and the number of atoms a homogeneous background charge is assumed. If the number of ions specified in the POSCAR file is 0 and NELECT=n, then the energy of a homogeneous electron gas is calculated.

[Voronoi CVs](https://github.com/Zhang-pchao/OilWaterInterface/tree/main/Ion_Diffusion_Coefficient) can be used to calculate the diffusion coefficient for H₃O⁺ or OH⁻ ions.
