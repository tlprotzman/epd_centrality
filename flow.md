# How Data Should Move
## Ingest
We start with hundreds of pico files.  There are preprocessed with PicoDstAnalyzer and simulationDataPreprocessor into singular root files.  Let's call these detector_data.root and sim_data.root.

## Stage 1 Analysis
From these bulk files containing all the data we need, we will then generate histograms using the different methods we are comparing.  These will all live in the root file epd_tpc_relations.root.

## Stage 2 Analysis
Next, we need to run the analysis actually comparing the methods.  For this, we need to find the quantiles for X and Y and compare equal quantiles projections onto the X axis.  To do this, we will open epd_tpc_relations.root and for each histogram complete the analysis.  I should see if I can find a way to do this without creating two new histograms.  That would save on memory, especially if I start running this on larger data sets.  

Ultimately what should be saved is a plot comparing the projections for each method so that they can be overlaid like in figure 11 of the paper, and the variance of each quantile range should be recorded so that a quantitative comparison between methods can be made.  This can be saved in method_comparison.root