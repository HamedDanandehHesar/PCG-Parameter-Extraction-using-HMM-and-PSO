# PCG Segmentation and Morphological Modeling using Springer HSMM and Dual‑Phase Analysis

MATLAB implementation for **heart sound (PCG) segmentation and morphological modeling** using the **Springer Logistic Regression–HSMM heart sound segmentation algorithm** combined with a **dual‑phase cardiac representation and Gaussian–cosine kernel modeling**.

The algorithm first detects **S1 and S2 heart sounds** using a state‑of‑the‑art HSMM segmentation method, then constructs a **cardiac phase representation**, estimates the **mean PCG morphology**, and finally models the signal using a **parametric kernel decomposition**.

This framework is useful for **heart sound analysis, cardiac cycle modeling, and synthetic PCG generation**.

---

# Overview

Heart sound signals (phonocardiograms – PCG) contain important information about cardiac mechanical activity.  
Key events in a cardiac cycle include:

• **S1** – closure of the mitral and tricuspid valves  
• **S2** – closure of the aortic and pulmonary valves  

Accurate detection of S1 and S2 enables:

• cardiac cycle segmentation  
• heart sound classification  
• murmur analysis  
• model‑based PCG signal reconstruction  

This project combines three important components:

1. **HSMM‑based heart sound segmentation** (Springer algorithm)
2. **Dual‑phase cardiac representation using S1 and S2**
3. **Parametric PCG morphology modeling using Gaussian‑cosine kernels**

---

# Methodology

The complete algorithm follows the pipeline below.

---

# 1. Loading PCG and ECG Data

The program loads a MATLAB `.mat` file containing:

• `PCG` → phonocardiogram signal  
• `ECG` → electrocardiogram signal (optional reference)  
• `fs` → sampling frequency  

Example:

```
data = load(file);
PCG = data.PCG;
ECG = data.ECG;
fs  = data.fs;
```

The ECG is only used for **visual comparison** with the PCG segmentation.

---

# 2. Heart Sound Segmentation (Springer Algorithm)

The algorithm uses the **Springer Logistic Regression–Hidden Semi‑Markov Model (HSMM)** for heart sound segmentation.

Reference:

D. Springer et al.  
“Logistic Regression‑HSMM‑based Heart Sound Segmentation”  
IEEE Transactions on Biomedical Engineering, 2015.

The segmentation divides the PCG signal into four physiological states:

1 — S1  
2 — systole  
3 — S2  
4 — diastole  

Implementation:

```
assigned_states = runSpringerSegmentationAlgorithm(...)
```

From the state sequence, S1 and S2 regions are extracted:

```
S1_segments = (assigned_states==1)
S2_segments = (assigned_states==3)
```

---

# 3. Detection of S1 and S2 Peaks

The segmentation output identifies **time regions** of S1 and S2.

To determine the exact peak location:

• a local search window is defined around each segment  
• the **maximum PCG amplitude** inside that window is selected  

This produces:

```
S1_peaks
S2_peaks
```

These peaks mark the main mechanical events of each cardiac cycle.

---

# 4. Dual‑Phase Cardiac Representation

To model the PCG signal in a periodic domain, a **dual‑phase representation** is constructed using S1 and S2.

The mapping is defined as:

S1 → phase = 0  
S2 → phase = π/2  
next S1 → phase = 2π

The phase is wrapped to:

```
[-π , π]
```

Two implementations are provided:

### Linear Phase Mapping

```
calculate_dual_phase_S1_S2_custom
```

Phase evolves linearly between S1 and S2.

### DTW‑Based Phase Mapping

```
calculate_dual_phase_S1_S2_DTW
```

Dynamic Time Warping is used to **nonlinearly align heart sound morphology**, producing a more physiologically consistent phase trajectory.

---

# 5. Mean PCG Morphology Extraction

The PCG signal is then averaged in the **phase domain**.

Procedure:

1. divide the phase interval into bins  
2. collect PCG samples belonging to each phase bin  
3. compute a robust average

A **robust Gaussian weighted mean** based on the **Median Absolute Deviation (MAD)** is used to reduce the influence of outliers.

Function used:

```
pcgsd_extractor_simple_ver_MAD_1
```

Outputs:

• `pcg_mean` → mean PCG waveform  
• `pcgsd` → phase‑dependent variance  

This produces the **average PCG morphology over one cardiac cycle**.

---

# 6. Kernel‑Based PCG Morphology Modeling

The mean PCG waveform is approximated using a **sum of Gaussian‑cosine kernels**.

Each kernel has the form

PCG(θ) = a · exp(-(θ − θᵢ)² / (2b²)) · cos(fθ − φ)

where

a = amplitude  
b = Gaussian width  
θᵢ = center phase  
f = modulation frequency  
φ = phase shift  

This representation captures both:

• localized energy of heart sounds  
• oscillatory structure of PCG components

---

# 7. Parameter Estimation using Particle Swarm Optimization

Kernel parameters are estimated sequentially using **Particle Swarm Optimization (PSO)**.

Optimization variables:

• amplitude (a)  
• width (b)  
• frequency (f)  
• phase shift (φ)

The cost function minimizes the difference between:

```
pcg_mean  –  modeled_kernel_sum
```

Typical configuration:

```
SwarmSize = 200
MaxIter   = 80
```

The algorithm iteratively adds kernels until the PCG morphology is well reconstructed.

---

# 8. Synthetic PCG Signal Reconstruction

After estimating the kernel parameters, a **synthetic PCG signal** is generated:

```
PCG_synth = Σ kernels(Phase)
```

This signal approximates the original PCG waveform using the parametric model.

The result demonstrates how well the extracted kernels represent the heart sound morphology.

---

# Outputs

The script produces several outputs and visualizations.

### Segmentation Plot
Displays:

• PCG signal  
• detected S1 segments  
• ECG reference

### Mean PCG Morphology
Average PCG waveform in the phase domain.

### Kernel Reconstruction
Comparison of:

• extracted mean PCG  
• kernel‑based reconstruction

### Synthetic PCG
Comparison of:

• original PCG  
• reconstructed PCG model

---

# MATLAB Requirements

The following MATLAB toolboxes may be required:

Signal Processing Toolbox  
Global Optimization Toolbox  
Statistics Toolbox  

External dependencies:

Springer Heart Sound Segmentation Toolbox

Required files include:

```
default_Springer_HSMM_options.m
runSpringerSegmentationAlgorithm.m
HMM_init_parameters_Springer.mat
```

These are part of the Springer segmentation framework.

---

# Applications

This framework can be used in:

• heart sound segmentation  
• PCG morphology modeling  
• heart sound synthesis  
• cardiac cycle analysis  
• murmur detection research  
• multimodal ECG–PCG studies

---

# Reference

Springer D., Tarassenko L., Clifford G.  
Logistic Regression‑HSMM‑based Heart Sound Segmentation  
IEEE Transactions on Biomedical Engineering, 2015.

