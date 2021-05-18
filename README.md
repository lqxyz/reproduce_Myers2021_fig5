## Reproduce the Fig. 5 of Myers et al. (2021) on $N$ature Climate Change

[Myers et al. (2021)](https://doi.org/10.1038/s41558-021-01039-0): Observational constraints on low cloud feedback reduce uncertainty of climate sensitivity

#### 
#### Data and code
* The forcing data (`.csv`) under [`data`](https://github.com/lqxyz/reproduce_Myers2021_fig5/tree/main/data) directory and python scripts under [`scripts`](https://github.com/lqxyz/reproduce_Myers2021_fig5/tree/main/scripts) direcotry (except `reproduce_myers2021_fig5.py`) are from [Sherwood et al. (2020)](https://doi.org/10.1029/2019RG000678) and public available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3945276.svg)](https://doi.org/10.5281/zenodo.3945276)

#### How to reproduce
* Run the following command:
```bash
./run_ecs.sh
```
* Or see the reproduced figure in [`reproduce_myers2021_fig5.ipynb`](https://github.com/lqxyz/reproduce_Myers2021_fig5/blob/main/reproduce_myers2021_fig5.ipynb). If Github can not load this notebook, you can try this [link](https://nbviewer.jupyter.org/github/lqxyz/reproduce_Myers2021_fig5/blob/main/reproduce_myers2021_fig5.ipynb) on nbviewer.

#### Input cloud feedback parameters 
Individual cloud feedbacks from Table 1 of Sherwood et al. (2020):

| Individual cloud feedbacks | Value (<img src="https://render.githubusercontent.com/render/math?math=Wm^{-2}K^{-1}">)|
| ------------- |:-------------:|
| High-cloud altitude | N(+0.20, 0.10) |
| Tropical marine low cloud | N(+0.25, 0.16) |
| Tropical anvil cloud area | N(-0.20, 0.20) |
| Land cloud amount | N(+0.08, 0.08) |
| Middle-latitude marine low-cloud amount | N(+0.12, 0.12) |
| High-latitude low-cloud optical depth | N(+0.00, 0.10) |
| Total cloud feedback | N(+0.45, 0.33) |

The `near-global marine low cloud feedback` in Sherwood as "the sum of tropical marine low cloud amount, midlatitude marine low cloud amount and high-latitude low cloud optical depth feedbacks", so the mean and stand deviation are: <img src="https://render.githubusercontent.com/render/math?math=\mu=0.25+0.12+0.00=0.37 Wm^{-2}K^{-1}"> and <img src="https://render.githubusercontent.com/render/math?math=\sigma=\sqrt{0.16^2 + 0.12^2 + 0.10^2}=0.22 Wm^{-2}K^{-1}">, respectively.

| | Sherwood et al. (2020) | Myers et al. (2021) |
| ------------- |:-------------:|:-------------:|
| near-global marine low cloud feedback| N(0.37, 0.22) | N(0.19, 0.07)  |
| Total cloud feedback | N(+0.45, 0.33) | N(+0.27, 0.25) |



