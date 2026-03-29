# Quick Start Guide

## Get the App Running in 5 Minutes

### Step 1: Install Required Packages (one-time only)

Open R or RStudio and run this:

```r
# Copy and paste this entire block into R console
packages <- c("shiny", "tidyverse", "plotly", "ordbetareg", "truncnorm", "shinydashboard", "brms")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}
```

It will ask to install from source - you can say "no" and use precompiled binaries.

### Step 2: Run the App

```r
# Set working directory
setwd("c:/Users/au645332/Documents/RDM-metacognition/Shiny app")

# Run the app
shiny::runApp()
```

The app should open automatically in your browser (or check RStudio viewer).

### Step 3: Use the App

1. **Left panel**: Adjust the parameter sliders
2. **Click "Simulate Data"** button
3. **View plots** in the center/right panels

That's it! You should see:
- A psychometric function (how choice probability changes with stimulus)
- Confidence ratings by accuracy (error monitoring)
- Summary statistics table

## First Thing to Try

Start with the default parameters and click "Simulate Data" to see how the model behaves.
Then:
- Increase evidence noise (σ_e) → flatter psychometric
- Increase choice noise (σ_k) → confident errors
- Increase metacognitive noise (σ_m) → worse error detection

## Need Help?

- See README.md for full documentation
- See model_description.html for model details (click "About Model" tab in app)
- Check the pre-registration documents for mathematical details

Enjoy exploring!
