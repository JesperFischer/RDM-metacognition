# RDM-metacognition

This repository contains code, data, and models for a preregistered study on metacognition in perceptual decision-making using a Random Dot Motion (RDM) task. The study investigates metacognitive uncertainty in correct and error trials across the post-decisional response window.

## Repository structure

```
Preregistration/       Preregistration documents (RMarkdown), figures, and supporting scripts
Experimental scripts/  PsychoPy experiment code for the RDM task
Pilot/                 Pilot data and experimental analysis
Simulations/           Simulation scripts for computational models
  Heurestic/           Heuristic model simulations
  MCMC/                MCMC-based model simulations
Sequential sampling/   Sequential sampling model materials
```

## Requirements

- **Experiment:** Python with PsychoPy
- **Analysis:** R (see `RDM-metacognition.Rproj`)
- **Models:** Stan (called via R)
