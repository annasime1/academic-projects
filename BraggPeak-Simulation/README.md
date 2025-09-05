# Simulation of Bragg Peak in Ion Beam Therapy

## Overview
This project focuses on the **simulation of the Bragg Peak** using **Geant4** to model proton beam interactions in **hadrontherapy**.  
The study reproduces the **CATANA proton therapy facility** setup (Catania, Italy) used for **ocular tumor treatment**, with the goal of analyzing how physical and biological doses vary under different simulation parameters.  

---

## Simulation Setup
- **Beam**: Proton, 60 MeV, ~10,000 particles.  
- **Phantom**: Water, 40 × 40 × 40 cm³.  
- **Detector**: 4 × 4 × 4 cm³ volume.  
- **Geant4 Monte Carlo simulation** reproducing CATANA facility components:  
  - Scattering system (beam enlargement).  
  - Collimators.  
  - Range shifters (energy modulation).  
  - Modulator wheel (Spread-Out Bragg Peak, SOBP).  
  - Monitor chambers + MOPI detector (beam symmetry check).  
  - Tumor-shaped patient collimator.  

---

## Parameters Studied
- **Number of particles**: dose amplitude ~linear with particle count.  
- **Beam energy**: linearly shifts Bragg Peak depth.  
- **Energy spread (σE)**:  
  - 0 MeV = ideal sharp peak.  
  - ~0.4 MeV = realistic optimal spread.  
  - >1 MeV = degradation of peak shape.  
- **Detector position**:  
  - Depth shifts → limited impact until far from peak.  
  - Transversal shifts → reduce dose, increase randomness.  

---

## Biological Considerations
- **Relative Biological Effectiveness (RBE)** used to compute biological dose:  
  `Biological Dose = RBE × Physical Dose`.  
- **Cell survival studies**:  
  - U87 Glioma cells → higher sensitivity, more DNA damage at the peak.  
  - AG01522 Healthy fibroblasts → higher survival fraction.  
- **Transition zone (~15–25 mm)**: balance between tumor control and risk of healthy tissue damage.  

---

## Key Findings
- Increasing proton energy shifts peak depth without changing shape.  
- Optimal energy spread is ~0.4 MeV.  
- Personalized therapy achievable by tuning **particle number** and **beam modulation**.  
- Healthy tissue near tumor edges is at risk, highlighting importance of precise beam control.  

---

## Future Work
- Refine SOBP modeling with advanced range modulators.  
- Simulate realistic tumor geometries and patient-specific cases.  
- Compare simulation results with experimental/clinical data for validation.  

---

## Further Details
For additional results and references, please refer to the **presentation slides** included in this repository.

