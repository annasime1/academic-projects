# Comparative Analysis of Chrominance-Based Techniques for Remote Photoplethysmography (rPPG)

## Overview
This project evaluates **chrominance-based** and **blind source separation (BSS)** techniques for **remote photoplethysmography (rPPG)**, a non-contact method to estimate heart rate from facial video recordings.  
The goal was to compare robustness, accuracy, and signal quality of different approaches under varying conditions of **motion and lighting**.

---

## Dataset
- Data: pre-extracted RGB signals from video recordings of **20 subjects**, ~600 seconds each.  
- Sampling rate: **115 Hz**.  
- Ground truth: reference **Blood Volume Pulse (BVP)** recorded simultaneously.  

### Preprocessing
- **Jump detection and correction** for signal discontinuities.  
- **Spline interpolation** for uniform resampling.  
- **Low-pass filtering** (4th-order Butterworth, 10 Hz cutoff).  

---

## Methods
- **Chrominance-based**: RoverG, XoverY, XoverY Fixed, XsminαYs.  
- **Blind Source Separation (BSS)**: Independent Component Analysis (ICA), Principal Component Analysis (PCA).  
- **Signal Processing**: bandpass filtering (40–240 BPM), spectral analysis with FFT, Bland–Altman agreement evaluation, SNR computation.  

---

## Results
- **Best overall performance**: **XsminαYs** (highest SNR, best agreement with reference BVP).  
- **PCA**: surprisingly robust, close to XsminαYs in many metrics.  
- **RoverG**: high correlation but lower robustness to noise.  
- **XoverY Basic & ICA**: weakest performance, highly sensitive to motion and illumination artifacts.  
- Methods like XsminαYs and PCA achieved **>96% agreement** with reference within ±5 BPM for most subjects.  

---

## Discussion
- Adaptive and statistically structured approaches (XsminαYs, PCA) were the most reliable across subjects and heart rate ranges.  
- Static methods (XoverY, ICA) degraded significantly under high motion or illumination variability.  
- Some failures were linked to **subject-specific variability** or **reference signal artifacts**.  
- Spectral analysis revealed limitations in handling harmonics, particularly at low heart rates.  

---

## Future Work
- Introduce **post-processing corrections** for harmonic misclassifications.  
- Expand dataset with controlled acquisition conditions (lighting, motion, skin tone).  
- Compare against modern **deep learning-based rPPG** methods.  

---

## Further Details
For the full methodology, analysis, and results, please see the [project report](docs/report.pdf).

