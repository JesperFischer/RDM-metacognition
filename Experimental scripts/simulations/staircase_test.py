from psychopy import data 
import numpy as np
from scipy.stats import norm
import sys
import os
import random
import json
import pandas as pd


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from RDM_exp import EXP_functions as exp

## Simulated vals 
lapse_vals = [0,0.03,0.1]

threshold_prior_1 = norm.pdf(np.linspace(1, 20, 50), loc=10, scale=5)
threshold_prior = threshold_prior_1 / threshold_prior_1.sum()

slope_prior_1 = norm.pdf(np.linspace(1, 10, 10), loc=5, scale=2)
slope_prior = slope_prior_1 / slope_prior_1.sum()

n = 40

sim = {}
## Staircase
for idx in range(8):
    slope = random.choice(np.linspace(1, 10, 10))
    thresh = random.choice(np.linspace(1, 20, 50))
    lapse = random.choice(lapse_vals)

    SC = data.QuestPlusHandler(nTrials=n,intensityVals=list(range(1, 20, 2)),thresholdVals=np.linspace(1, 20, 50),slopeVals=np.linspace(1, 10, 10),lowerAsymptoteVals=0.5,
    lapseRateVals=0.03,responseVals=[1, 0],prior={"threshold": threshold_prior, "slope": slope_prior},psychometricFunc="weibull",startIntensity=11,
    stimScale="linear",stimSelectionMethod="minEntropy",paramEstimationMethod="mean")


    df = pd.DataFrame(columns = ['C', 'out', 'alpha', 'beta', 'lapse'])

    for i in range(n):
        print(i)
        stim = SC.next()

        print(stim)

        P = exp.weibull_cdf(stim, alpha = thresh, beta = slope, gamma =0.5, lapse = lapse)
        out = np.random.binomial(n=1, p=P)

        SC.addResponse(out)

        df.loc[len(df)] = [stim, out, thresh, slope, lapse]

    thresh_est = SC.paramEstimate["threshold"]
    slope_est = SC.paramEstimate["slope"]

    df["slope_est"] = slope_est
    df["thresh_est"] = thresh_est

    sim[idx] = df 

sim_serializable = {k: v.to_dict(orient="records") for k, v in sim.items()}

with open("Experimental scripts/simulations/sim_data.json", "w") as f:
    json.dump(sim_serializable, f, indent=4)
