In this update to the pre-registration, we report three changes.

First, we identified a bug in the experimental code.
Specifically, the "DotLife" variable was staircased at a coherence level of 0.9 instead of the intended 0.3.
This error occurred due to an overwritten variable and has now been corrected.
We believe this bug is partly responsible for the issue described in the second update.

Second, we now exclude trials with the most extreme coherence values for participants whose performance drops substantially on these trials.
We suspect that this pattern is partly driven by the bug described above.
Staircasing "DotLife" at a high coherence level (0.9) made these trials easier than intended, which in turn caused the staircase procedure to reduce "DotLife"" to compensate for the increased ease.
This combination of low "DotLife" and high coherence, may have induced perceptual artifacts similar to the reverse-phi illusion, potentially leading to responses that look like guessing, but at coherence levels where accuracy should be close to perfect.
As this is a very different perceptual process we decided to exclude these trials.


Third, we were unable to achieve convergence of our hierarchical Bayesian model for the initial 15 participants using our original sampler settings and contingency steps.
To proceed with evaluating our stopping criteria, we removed the lapse-rate parameter from the model.
The main issue appeared to be a strong trade-off between the lapse-rate and the noise parameters governing choice ($sigma_e$ and $sigma_k$), which prevented stable estimation.
Removing the lapse-rate resolved the convergence issues. Therefore, we add removal of the lapse-rate as an explicit contingency step for achieving convergence in our future testing.

