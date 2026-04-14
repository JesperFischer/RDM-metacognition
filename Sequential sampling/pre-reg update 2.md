In this update to the pre-registration, we report a single modification to the model. After successfully fitting data from participants 15 through 25, we encountered convergence issues when extending the hierarchical Bayesian model to 30 participants using our existing sampler settings and updated contingency procedures. Specifically, we observed a couple of divergent transitions.

Upon investigation, we identified numerical instability in the latent mean of the confidence ratings as a likely source of the problem. In particular, values approaching the boundaries (0 or 1) induced regions of high curvature in the posterior distribution, which in turn led to the observed divergences.

To mitigate this issue, we introduced a constraint on the latent confidence mean, ensuring it remains within a bounded interval away from the extremes. Concretely, we modified the Stan model to clamp the latent confidence mean to the interval [0.001, 0.999] immediately after its computation, this addition to the Stan model was implemented as follows just after the confidence mean computation:


if (theta_conf[n] > 0.999) {
  theta_conf[n] = 0.999;
} else if (theta_conf[n] < 0.001) {
  theta_conf[n] = 0.001;
}

This type of adjustment of numerical stabilization technique is sometimes used in Bayesian modeling, to prevent underflow or overflow and to avoid pathological posterior geometries near the boundaries of constrained parameters. Importantly, this modification does not alter the substantive interpretation of the model, but ensures stable and reliable inference.