## Reproduce the Fig. 5 of Myers et al. (2021) on Nature Climate Change

[Myers et al. (2021)](https://doi.org/10.1038/s41558-021-01039-0): Observational constraints on low cloud feedback reduce uncertainty of climate sensitivity

#### 
#### Data and code
* The forcing data (`.csv`) under `data` directory and python scripts under `scripts` direcotry (except `reproduce_myers2021.py`) are from [Sherwood et al. (2020)](https://doi.org/10.1029/2019RG000678) and public available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3945276.svg)](https://doi.org/10.5281/zenodo.3945276)

#### How to reproduce
```bash
./run_ecs.sh
python -u ./scripts/reproduce_myers2021.py
```
