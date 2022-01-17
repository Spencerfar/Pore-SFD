# Single-file diffusion in pores
Code for the paper "Non-Fickian single-file pore transport" https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.L032102 (pre-print https://arxiv.org/abs/2107.03498).

<p align="center"> 
<img src="diagram.png" width="714" height="300">
</p>

## Running simulations
Simulation code in Sim/ is compiled with "make".

To run simulations for density profiles and flux (figure 2), use main.cpp with command line arguments for L, Ron, dx, D0, koff, run ID. Example: ./main 8 0.1 7 270000 0.1 0.

To run simulations for time series (figure 3 and 4), use time_series.cpp with command line arguments for L, Ron, dx, D0, koff, run ID, scale_series, and scale_numtime. scale_series sets when time for which the time series is recorded, at time=2^scale_series. scale_numtime is the number of time steps the series is recorded for.

## Plotting
Plotting code in PlottingCode/ creates the plots from data from the simulations.

## Citation
The original paper should be cited as
```
@ARTICLE{Farrell2021-sd,
  title     = "{Non-Fickian} single-file pore transport",
  author    = "Farrell, Spencer and Rutenberg, Andrew D",
  journal   = "Phys. Rev. E",
  publisher = "American Physical Society",
  volume    =  104,
  number    =  3,
  pages     = "L032102",
  month     =  sep,
  year      =  2021
}

```
