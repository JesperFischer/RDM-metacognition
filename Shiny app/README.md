# Heuristic Model Simulator - Shiny App

## Overview

This Shiny application allows you to interactively simulate data from the heuristic model of perceptual decision making and metacognition, and visualize how model parameters affect:

1. **Psychometric functions** (relationship between stimulus strength and choice probability)
2. **Confidence ratings** (mean confidence by accuracy level)

## Installation

### Requirements

Make sure you have R installed (version 4.0 or later recommended).

### Install Dependencies

Run the following in R:

```r
# Install required packages
packages <- c(
  "shiny",
  "tidyverse", 
  "plotly",
  "ordbetareg",
  "truncnorm",
  "shinydashboard",
  "brms"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}
```

## Usage

### Running the App

1. Open R or RStudio
2. Set working directory to the Shiny app folder:
   ```r
   setwd("c:/Users/au645332/Documents/RDM-metacognition/Shiny app")
   ```
3. Run the app:
   ```r
   shiny::runApp()
   ```
4. The app will open in your default browser or in RStudio's viewer pane

### Using the Simulator

**Left Panel - Model Parameters:**
- **Evidence Noise (σ_e)**: Perceptual noise in the evidence [on log scale]
  - Lower values (e.g., -3) = less noise = better discrimination
  - Higher values (e.g., 0) = more noise = worse discrimination

- **Criterion Variability (σ_k)**: Stochasticity in the decision criterion [log scale]
  - Controls how noisy/variable the decision process is
  - Creates choice errors independent of evidence quality

- **Metacognitive Noise (σ_m)**: Noise in introspection [log scale]
  - Controls access to internal evidence for confidence
  - Lower = better metacognitive sensitivity

- **Stimulus Bias (β)**: Systematic bias in evidence
  - Shifts psychometric function left/right

- **Confidence Precision**: How extreme the confidence ratings are
  - Lower values = more uniform, less extreme
  - Higher values = more extreme (closer to 0 or 1)

- **Number of Trials**: Total simulation length

- **Number of Stimulus Levels**: How many different stimulus strengths to use

**Center/Right Panels - Visualizations:**

1. **Psychometric Function**: 
   - X-axis: Stimulus strength
   - Y-axis: Probability of choosing "1"
   - Points colored by accuracy (correct = green, error = red)
   - Shows slope (discrimination) and threshold

2. **Confidence Ratings**:
   - X-axis: Stimulus strength
   - Y-axis: Mean confidence rating (0-1)
   - Colored by accuracy
   - Shows error monitoring and calibration

3. **Summary Statistics Table**:
   - Displays current parameter settings
   - Shows aggregate performance metrics

## Tips for Exploring

### Understanding Parameters

- **Effect of σ_e on psychometric**: Flatten the curve (increase discrimination noise)
- **Effect of σ_k on confidence**: Create "confident errors" - wrong choices with high confidence
- **Effect of σ_m on metacognition**: Weak or strong relationship between confidence and accuracy
- **Effect of β on both**: Shifts both psychometric and confidence curves

### Interesting Manipulations

1. **Compare signal detection models**: Set σ_k to very low (-3) vs high (0)
   - Low σ_k = deterministic choices based mostly on evidence
   - High σ_k = random criterion = more choice errors with good evidence

2. **Test error monitoring**: Set σ_m to very low (-3) 
   - Watch how well the model recognizes errors (low confidence on errors)

3. **Simulate individual differences**: Try different parameter combinations
   - Some people might have high σ_e (perceptually noisy)
   - Others might have high σ_m (poor introspection)

## File Structure

```
Shiny app/
├── app.R                      # Main Shiny application
├── simulation_functions.R     # Simulation and utility functions
├── model_description.html     # Model explanation for About tab
└── README.md                  # This file
```

## Model Details

The model implements the following process:

1. **Evidence Generation**: `e ~ N(D·X - β, σ_e)`
2. **Variable Criterion**: `c ~ N(0, σ_k)`  
3. **Binary Choice**: `choice = 1 if e > c else -1`
4. **Confidence Readout**: `r ~ N(e, σ_m)`
5. **Confidence Rating**: `P(correct|r, choice)` via logistic function

For full mathematical derivation, see the pre-registration documents in the project.

## Troubleshooting

### App won't start
- Check that all packages installed correctly
- Try: `library(shiny)` in R console to verify installation

### Plots not showing
- Click "Simulate Data" button to generate initial data
- Check that you're seeing white/blank panels (not error messages)

### Slow performance
- Reduce "Number of Trials" to speed up simulation
- Close other applications to free up RAM

## Contact

For questions about the model or app, refer to the project documentation or the pre-registration files.
