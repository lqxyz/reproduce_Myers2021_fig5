## Reproduce the Fig. 5 of Myers et al. (2021) on Nature Climate Change

[Myers et al. (2021)](https://doi.org/10.1038/s41558-021-01039-0): Observational constraints on low cloud feedback reduce uncertainty of climate sensitivity

#### 
#### Data and code
* The forcing data (`.csv`) under [`data`](https://github.com/lqxyz/reproduce_Myers2021_fig5/tree/main/data) directory and python scripts under [`scripts`](https://github.com/lqxyz/reproduce_Myers2021_fig5/tree/main/scripts) direcotry are from [Sherwood et al. (2020)](https://doi.org/10.1029/2019RG000678) and public available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3945276.svg)](https://doi.org/10.5281/zenodo.3945276)

#### Input cloud feedback parameters 
* Individual cloud feedbacks from Table 1 of Sherwood et al. (2020) are shown as follows and N(mu, sigma) is for Gaussian distribution:

| Individual cloud feedbacks | Value (<img src="https://render.githubusercontent.com/render/math?math=Wm^{-2}K^{-1}">)|
| ------------- |:-------------:|
| High-cloud altitude | N(0.20, 0.10) |
| Tropical marine low cloud | N(0.25, 0.16) |
| Tropical anvil cloud area | N(-0.20, 0.20) |
| Land cloud amount | N(0.08, 0.08) |
| Middle-latitude marine low-cloud amount | N(0.12, 0.12) |
| High-latitude low-cloud optical depth | N(0.00, 0.10) |
| Total cloud feedback | N(0.45, 0.33) |

* The `near-global marine low cloud feedback` in Sherwood et al. (2020) is defined as "the sum of tropical marine low cloud amount, midlatitude marine low cloud amount and high-latitude low cloud optical depth feedbacks", so the mean and standard deviation (sigma) are: `mu=0.25 (tropical marine low cld) + 0.12 (mid-lat marine low cld) + 0.00 (high-lat low cld optical depth) = 0.37 (near-global marine low cloud)` and <img src="https://render.githubusercontent.com/render/math?math=\sigma=\sqrt{0.16^2 %2B 0.12^2 %2B 0.10^2}=0.22 Wm^{-2}K^{-1}">, respectively.
* In Myers et al. (2021), the 90% confidence rather than sigma is used, and 90% confidence should be divided by `3.29/2=1.645` to get the sigma, e.g., `0.37/1.645=0.22`, `0.12/1.645=0.07` and `0.33*1.645=0.54`. 
* The best estimate of `near-global marine low cloud feedback` from Myers et al. (2021) is `N(0.19, 0.07)`, so if we replace the "sum of their tropical marine low cloud amount, midlatitude marine low cloud amount and high-latitude low cloud optical depth feedbacks" with 0.19, and keep all the other individual cloud feedback terms unchanged (high-cloud altitude, tropical anvil cloud area, and land cloud amount), the mean of the total cloud feedback should be: `0.19 (near-global marine low cloud) + 0.2 (altitude) - 0.2 (anvil) + 0.08 (land) = 0.27 (total)`, and the sigma should be: <img src="https://render.githubusercontent.com/render/math?math=\sigma=\sqrt{0.07^2 %2B 0.10^2 %2B 0.20^2 %2B 0.08^2}=0.25 Wm^{-2}K^{-1}">
* The values in the last row (i.e., mean and standard deviation of the total cloud feedback) of follwing table are used in the script ([`ecs-baseline-ec2.py`](https://github.com/lqxyz/reproduce_Myers2021_fig5/blob/main/scripts/ecs-baseline-ec2.py#L63-L70)) to generate the data for reproducing the figure.

| | Sherwood et al. (2020) | Myers et al. (2021) |
| ------------- |:-------------:|:-------------:|
| Near-global marine low cloud feedback| N(0.37, 0.22) | N(0.19, 0.07)  |
| Total cloud feedback | N(0.45, 0.33) | N(0.27, 0.25) |

#### How to reproduce
* Run the following command:
```bash
$ ./run_ecs.sh
```
* Also, you can see the reproduced figure in [`reproduce_myers2021_fig5.ipynb`](https://github.com/lqxyz/reproduce_Myers2021_fig5/blob/main/reproduce_myers2021_fig5.ipynb). If Github can not load this notebook, try this [link](https://nbviewer.jupyter.org/github/lqxyz/reproduce_Myers2021_fig5/blob/main/reproduce_myers2021_fig5.ipynb) on nbviewer.

