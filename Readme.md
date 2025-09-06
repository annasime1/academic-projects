# Academic Projects Portfolio

This repository collects a selection of academic projects developed during my studies in **Biomedical Engineering**.  
Each folder contains a self-contained project with its own README, code or notebooks, figures, and, when available, reports or presentations.

---

## Project Index

- **[Blood-samples-classification](Blood-samples-classification/)**  
  Deep learning for microscopic blood cell images (8 classes). Includes baselines and transfer learning models (EfficientNet, VGG) with class-imbalance handling and data augmentation.

- **[LungCancer-Classification](LungCancer-Classification/)**  
  CT-based lung nodule classification (binary and multi-class). The project compares custom CNNs with transfer learning (EfficientNetB2), and includes model explainability using Grad-CAM and LIME.

- **[Mars-terrain-segmentation](Mars-terrain-segmentation/)**  
  Semantic segmentation of Martian terrain using U-Net architectures with ASPP and attention mechanisms. Experiments include different depths, augmentations, and loss functions.

- **[Pulse-rate-extraction](Pulse-rate-extraction/)**  
  Remote photoplethysmography (rPPG) from facial video signals. Comparative analysis of chrominance-based methods and blind source separation (PCA, ICA) for heart rate estimation.

- **[Attention-ASD](Attention-ASD/)**  
  Attention classification for children with Autism Spectrum Disorder using gaze features collected during sessions with a social assistive robot. Implemented with an MLP and balancing strategies.

- **[BraggPeak-Simulation](BraggPeak-Simulation/)**  
  Geant4 Monte Carlo simulation of the Bragg Peak in proton therapy, reproducing the CATANA setup for ocular tumors. Analysis includes energy, spread, number of particles, detector shifts, and biological effectiveness.

- **[Laughter-Effects](Laughter-Effects/)**  
  Study of cardio-respiratory effects during laughter and after a program of laughter therapy. Measurements include ventilation, lung volumes, heart rate, and heart rate variability.

- **[Telemonitoring-Hypertension](Telemonitoring-Hypertension/)**  
  Design of a telemedicine system for hypertensive patients. Includes multi-user workflows, clinical data acquisition, questionnaires, report handling, risk scoring, and a patient forum.

- **[diabetic-readmission-ML](diabetic-readmission-ML/)**  
  Machine learning pipeline for predicting hospital readmissions in diabetic patients. Includes data preprocessing, baseline models, imbalance handling, and feature importance analysis.

---

## Technologies Used

Across the projects, the following tools and technologies are employed:

- Python (PyTorch, TensorFlow, scikit-learn)  
- Jupyter Notebooks  
- MATLAB  
- SQL and data modeling  
- C++/Geant4  
- Classical signal processing and statistical methods  

---

## How to Navigate

1. Each project folder contains a dedicated README file describing the objectives, methodology, and results.  
2. Code and notebooks are organized within each project.  
3. Where datasets cannot be shared, instructions or scripts for downloading public datasets are provided.  

---

## Repository Conventions

- Reports and presentations are stored in a `docs/` folder inside each project.  
- Figures and plots are saved in a `results/` folder.  
- Small data samples or synthetic datasets are stored in a `data/` folder.  
- Each README follows a consistent structure: *Overview → Data → Method → Results → References*.  

---

## Contact

For questions or collaborations, feel free to open an issue or reach out through my GitHub profile.
