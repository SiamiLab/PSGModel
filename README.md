# PSGModel

MATLAB code and videos for the paper

> A. Darabi, H. Montazeri Hedesh, M. Siami, and M. Sznaier,
> **"Predator-Swarm-Guide Dynamics: A Hybrid Approach to Crowd Modeling and Guidance in Mass Shooting Scenarios,"**
> *2024 American Control Conference (ACC)*, Toronto, ON, Canada, 2024, pp. 4402-4408.
> DOI: [10.23919/ACC60939.2024.10645036](https://doi.org/10.23919/ACC60939.2024.10645036)

The Predator-Swarm-Guide (PSG) model describes how a crowd moves during an active-shooter event.
People repel and attract each other, move away from the shooter (the predator), and follow a
trusted guide toward a safe zone, inside a zoned building with walls. The paper uses the model to
choose the crowd and guidance parameters that reduce casualties, and studies where the crowd
settles with a continuum-limit version of the model.

## Contents

| File | What it does |
|---|---|
| `Paper1-LCSS Matlab codes/Evacuation_OL.m` | Evacuation in the baseline case (`alpha1 = 1`, `lambda = 0.1`). Writes `Evacuation_RK4_OL_alpha1_lambda0.1.mp4`. |
| `Paper1-LCSS Matlab codes/Evacuation_CL.m` | Evacuation in the optimized case (`alpha1 = -1.625`, `lambda = 0.7135`). Writes `Evacuation_RK4_CL_alpha_1.625_lambda0.7135.mp4`. |
| `Paper1-LCSS Matlab codes/Equilibrium_A.m` to `Equilibrium_F.m` | Equilibrium of the crowd around the shooter for six parameter settings, with the inner and outer ring radii predicted by the continuum-limit analysis. |
| `Paper1-LCSS Matlab codes/x0.mat` | Initial positions used by the two evacuation scripts. |
| `Paper1-LCSS Videos/` | The two evacuation videos (baseline and optimized). |

Run the scripts from inside `Paper1-LCSS Matlab codes/` so that they find `x0.mat`.

## Citation

If you use this code, please cite the paper. GitHub's "Cite this repository" button (from
[`CITATION.cff`](CITATION.cff)) gives the same reference in other formats.

```bibtex
@inproceedings{darabi2024predator,
  author    = {Darabi, Atefe and Montazeri Hedesh, Hamidreza and Siami, Milad and Sznaier, Mario},
  title     = {Predator-Swarm-Guide Dynamics: A Hybrid Approach to Crowd Modeling and Guidance in Mass Shooting Scenarios},
  booktitle = {2024 American Control Conference (ACC)},
  address   = {Toronto, ON, Canada},
  pages     = {4402--4408},
  year      = {2024},
  doi       = {10.23919/ACC60939.2024.10645036}
}
```
