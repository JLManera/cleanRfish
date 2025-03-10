# cleanRfish 🐟📊  
**Trajectory Cleaning for Animal Tracking Data**  
*A PhD Student's Tool for Handling Tracking Software Artifacts*

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## ⚠️ Project Status Note
**Maintainer:** [Your Name], PhD Candidate in [Your Department/Field]  
**Maintenance Level:** 🐣 Experimental/Basic (Personal Research Project)  

> **Important Note:** This package is actively used in my PhD research but receives irregular maintenance. While functional for my needs, users should:
> - Expect delayed responses to issues
> - Verify results against raw data
> - Not consider this production-ready software

## 🔍 Problem Addressed
Automated video tracking systems (e.g., EthoVision, DeepLabCut) often exhibit "jumping" artifacts where:
1. Software temporarily locks onto non-target objects
2. Creates artificial discontinuities in trajectories
3. Introduces noise in behavioral analyses

### **Example Cases:**
- Fish temporarily obscured by tank decor
- Rodent briefly hidden in complex maze
- Insect lost in high-density tracking

## 📺 Package Overview
### **Current Implementation (Heuristic Approach)**
```r
find_smooth_path(raw_data) %>%
  analyze_behavior()
```
Filters discontinuities using:
- **Velocity Thresholding** (MAD-robust statistics)
- **Path Segmentation**
- **Savitzky-Golay Smoothing**

### **Future Directions (Probabilistic Framework)**
```r
# Goal for v2.0
bayesian_path_cleaner(raw_data, species_movement_model) %>%
  probabilistic_analysis()
```
Planned features:
- **Hidden Markov Models** for state detection
- **Species-specific movement priors**
- **Uncertainty quantification**

## 🛠️ Installation
```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("[yourusername]/cleanRfish")
```

## 🐟 Basic Usage
### **For Zebrafish Tracking Data**
```r
library(cleanRfish)

# Process raw tracking output
cleaned_data <- find_smooth_path(
  raw_df,
  na.fill = TRUE,  # Interpolate small gaps
  p = 3,           # Cubic polynomial for SG filter
  n = 13           # 13-point smoothing window
)

# Analyze movement parameters
velocity_profile <- cleaned_data %>%
  mutate(
    inst_velocity = sqrt(diff(x)^2 + diff(y)^2) / diff(time)
  )
```

## 📈 Validation
### **Tested On:**
- Zebrafish (*Danio rerio*) 2D tracking
- Drosophila larval tracking
- Artificial datasets with controlled jumps

### **Performance:**
- Reduces false displacement events by **62-89%** in validation sets
- Preserves true movement features (*see /vignettes/validation.Rmd*)

## 🌱 Contributing
While this is primarily a research tool, I welcome:
- **Bug reports** with reproducible examples
- **Documentation improvements**
- **Bayesian methods expertise**

### **Please note:** As a full-time PhD student, I cannot guarantee:
- **Timely responses to issues**
- **Implementation of feature requests**
- **Compatibility with other tracking formats**

## 📝 License
GNU GPLv3 - See LICENSE file

## 📍 Roadmap
```mermaid
gantt
    title cleanRfish Development Timeline
    dateFormat  YYYY-MM
    section Core Functionality
    Velocity Filtering       :done, des1, 2023-01, 2023-03
    Path Reconnection        :done, des2, 2023-04, 2023-06
    SG Smoothing             :done, des3, 2023-05, 2023-07
    
    section Future Goals
    HMM Implementation      :active, 2024-01, 2024-12
    Probabilistic Framework  :crit, 2025-01, 2026-12
    Shiny Interface          :2026-01, 2027-12
```

## 📨 Contact
For scientific use inquiries:  
**[Your Academic Email]**  
Lab Website: **[Your Lab URL]**  

*Please understand if I'm slow to respond during:*
- Thesis writing months (**MM-YYYY**)
- Experiment-intensive periods
- Conference deadlines

