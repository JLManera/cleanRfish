# cleanRfish 🐟📊  
**Trajectory Cleaning for Animal Tracking Data**  
*A PhD Student's Tool for Handling Tracking Software Errors*

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## ⚠️ Project Status Note
**Maintainer:** Jack L Manera, PhD Candidate in the School of Biological Sciences, Monash University 

**Maintenance Level:** 🐣 Experimental/Basic (Personal Research Project)  

> **Important Note:** This package is actively used in my PhD research but receives irregular maintenance. While functional for my needs, users should:
> - Expect delayed responses to issues
> - Verify results against raw data

## 🔍 Problem Addressed
Automated video-based tracking systems (e.g., EthoVision, idtracker.ai) are powerful tools, but they routinely produce systematic, non-random errors that disrupt otherwise continuous animal trajectories. Unlike GPS-based movement ecology—where errors are typically small, independent, and easily removed by smoothing—high-resolution video tracking produces clusters of contiguous errors that cannot be treated as stochastic noise.

These problems arise because video tracking systems may:

- momentarily assign the focal animal’s identity to another object,

- track a reflection, ripple, shadow, or background artefact,

- lose the focal animal during occlusions or interactions, and then fail to correctly reacquire the individual after it reappears.

The result is that tracking “errors” may look biologically credible, forming long, smooth, and coherent false trajectories. Standard filtering or smoothing cannot distinguish these artefactual segments from genuine movement because the errors are temporally autocorrelated and often span many frames.

In practice, this produces trajectories that are:

- fragmented into many disjoint segments,

- interrupted by identity swaps or missed detections, and

- unusable for behavioural or kinematic analysis unless they can be reliably reconstructed.

cleanRfish addresses this problem by providing a principled method to identify, rank, and reconnect fragmented trajectory segments using uncertainty-aware likelihood models, rather than relying on conventional smoothing. By modelling how positional uncertainty grows over time, predicting feasible future positions, and penalising implausible kinematic transitions, cleanRfish reconstructs trajectories that would otherwise remain broken or misleadingly corrupted.

The goal is to clean the data and differentiate real behaviour from tracking artefacts, ensuring that downstream biological inferences are based on authentic movement rather than software errors. 

That said, it’s worth noting that the need for a package like cleanRfish will hopefully diminish as tracking software continues to improve. Many controlled laboratory systems already achieve excellent accuracy with minimal errors. My intention here is simply to offer a practical tool for situations where tracking is pushed beyond those ideal conditions—such as ecological or multi-individual assays in visually complex environments—where errors are still a routine challenge. I share the hope that future advances in tracking technology will make tools like this increasingly obsolete.

## 📺 Package Overview
*cleanRfish* provides a comprehensive pipeline for removing errors from animal tracking trajectories. 

*cleanRfish* reconstructs animal tracking trajectories using bidirectional, uncertainty weighted spline propagation. The algorithm first fragments the track by detecting abnormal changes in the animal’s velocity and heading direction. These flags partition the track into segments that represent periods of both contiguous tracking of an object (either the focal animal or a non-target object) and stochastic errors. Large temporal gaps in the tracking data (exceeding the user specified search window) are then identified and used to group segments into isolated tracking groups, ensuring the algorithm does not link across biologically implausible time intervals. Within each group, a reliable ‘ground truth’ segment is selected, and from this anchor point the algorithm propagates both forwards and backwards in time using multivariate GAM splines to predict trajectory continuations, ranking candidate segments by log likelihood scores derived from Ornstein–Uhlenbeck uncertainty modelling. This uncertainty model is constructed through bootstrap sampling of observation pairs at varying time separations, quantifying how positional uncertainty grows with temporal gap, and when taken together with the spline predicted trajectory is used to select the next most probable segment. This procedure is iterated until the full path is reconstructed and segments deemed unlikely are excluded, subject to biological plausibility constraints (velocity limits informed by the data). Reconstructed trajectories can optionally undergo temporal smoothing via Savitzky–Golay or moving average filtering after linear interpolation of remaining gaps. For multi individual datasets, the complete pipeline is run independently for each individual and the results are combined into a single data frame. Visualisation tools include diagnostic plots showing spline fits with uncertainty envelopes and candidate rankings, spatial 2D trajectory plots, and video overlay generation using ffmpeg to burn colour coded tracks onto the source video, while an interactive Shiny gadget allows chronological browsing of diagnostic plots for reconstruction quality assessment.

## 🛠️ Installation
```r
# Install from GitHub
remotes::install_github("JLManera/cleanRfish")
```

## 🐟 Basic Usage
### **For Fish Tracking Data**
```r
library(cleanRfish)

# Load tracking data (expected format: time, x1, y1, x2, y2, ...)
df_multi <- read.csv("trajectories.csv")

# Process all individuals with uncertainty-weighted reconstruction
results <- clean_track(
  df_raw = df_multi,
  window = 20,                        # 20-second search window
  min_candidates = 8,                 # Minimum candidates segments to consider
  speed_threshold_quantile = 0.999,   # Speed feasibility threshold (99.9%)
  smooth_method = "savitzky_golay",   # Smoothing filter
  diagnostic_plots = FALSE            # Set TRUE for visualization
)

# Extract reconstructed tracks
all_individuals <- results$tracks

# View structure
head(all_individuals)
# Columns: time, x_original, y_original, x_reconstructed, y_reconstructed, 
#          x_smooth, y_smooth, individual_number
```

### Visualise Results

```r
library(ggplot2)

# Plot composite trajectory for Individual 1
ggplot() +
  # Original fragmented data
  geom_point(
    data = all_individuals |> filter(individual_number == 1),
    aes(x = time, y = y_original + x_original),
    alpha = 0.4, colour = "black", size = 0.5
  ) +
  # Reconstructed continuous trajectory
  geom_path(
    data = all_individuals |> 
      filter(individual_number == 1, !is.na(x_reconstructed)),
    aes(x = time, y = x_reconstructed + y_reconstructed),
    colour = "purple", linewidth = 0.6
  ) +
  labs(
    title = "Trajectory Reconstruction - Individual 1",
    x = "Time (seconds)",
    y = "Composite Position (X + Y)"
  ) +
  theme_classic()
```

### Interactive Diagnostic Plots

```r
# Generate diagnostic plots during processing
results <- clean_track(
  df_raw = df_multi,
  diagnostic_plots = TRUE
)

# Launch interactive viewer
check_spline_fits(results, individual_number = 1)
# Use arrow keys (← →) to navigate through reconstruction steps
```

### Create Video Overlay

```r
# Overlay all individuals simultaneously on original video
overlay_video(
  tracks_df = all_individuals,
  video_path = "path/to/video.mp4",
  output_name = "reconstructed_overlay.mp4",
  use_smoothed = TRUE
)
```

## Core Functions

The following functions are available for users:

### Main Functions

| Function | Description |
|----------|-------------|
| `clean_track()` | Primary trajectory cleaning function - reconstructs fragmented tracks using temporal priority ground truth selection |
| `symmetrical_clean_track()` | Alternative reconstruction using symmetrical (largest segment) ground truth selection |
| `smooth_track()` | Apply smoothing to tracking data without reconstruction |
| `check_spline_fits()` | Interactive Shiny viewer for diagnostic plots |
| `overlay_video()` | Create video with reconstructed tracks overlaid |

## Parameter Tuning

| Parameter | Default | Description |
|-----------|---------|-------------|
| `window` | 20 | Search window size (seconds). Should match typical occlusion duration |
| `min_candidates` | 5 | Minimum segments to consider as the next potential candidate segment of the track. Window auto-expands if needed |
| `speed_threshold_quantile` | 0.999 | Speed feasibility quantile. Higher = more permissive |
| `min_movement` | 50 | Minimum total movement for ground truth selection |
| `smooth_method` | "savitzky_golay" | Smoothing filter: "savitzky_golay" or "moving_average" |

## Data Format

### Input Format

CSV file with columns:
- `time`: Timestamp (numeric, in seconds)
- `x1, y1, x2, y2, ..., xN, yN`: Coordinates for N individuals

```
time,x1,y1,x2,y2,x3,y3
0.00,245.3,189.7,512.1,334.2,678.9,201.4
0.04,246.1,190.2,513.4,335.1,679.3,202.1
...
```

### Output Format

Data frame with columns per individual:
- `time`: Timestamp
- `x_original, y_original`: Original raw coordinates
- `x_reconstructed, y_reconstructed`: Gap-filled coordinates
- `x_smooth, y_smooth`: Smoothed final coordinates
- `individual_number`: Individual identifier (1 to N)


## 📈 Validation
### **Tested Systems:**

- **Species**: Guppies (*Poecilia reticulata*), Eastern mosquitofish (*Gambusia holbrooki*), but is general applicabile to 2D tracking of any animal
- **Tracking software**: idtrackerai, EthoVision XT
- **Data scale**: Up to 500K points per individual, however more is possible it just takes a long time. 

## Requirements

- R >= 4.1
- Required packages: dplyr, mgcv, signal, zoo, ggplot2
- Optional: plotly, shiny, miniUI (for interactive features)
- ffmpeg (for video overlay functionality)

## Citation

If you use cleanRfish in your research, please cite:

```
Manera, J.L. (2025). cleanRfish: Uncertainty-Weighted Trajectory 
Reconstruction for Animal Tracking. R package version 2.0.0. 
https://github.com/JLManera/cleanRfish
```

## 🌱 Contributing
While this is primarily a research tool, I welcome:
- **Bug reports** with reproducible examples
- **Documentation improvements**

### **Please note:** As a full-time PhD student, I cannot guarantee:
- **Timely responses to issues**
- **Implementation of feature requests**
- **Compatibility with other tracking formats**

## 📝 License
GNU General Public License v3.0 - See LICENSE file

## 📨 Contact
For scientific use inquiries:  
**jack.manera@monash.edu**   
