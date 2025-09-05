# Blood Cell Classification with Deep Learning

## Overview
This project applies **deep learning techniques** to classify microscopic blood cell images into **eight categories**:  
**Basophil, Eosinophil, Erythroblast, Immature Granulocytes, Lymphocyte, Monocyte, Neutrophil, Platelet.**

The goal is to develop robust models for medical diagnostics, starting from a **custom CNN** and extending to **transfer learning** with **VGG16** and **EfficientNet** variants.

---

## Dataset
- ~13,000 blood cell images, resolution **96×96 RGB**.  
- Eight cell classes, with **imbalanced distribution**.  
- Split: **80% training, 10% validation, 10% testing**.  

### Challenges
- **Class imbalance**: some cell types are underrepresented.  
- **Limited variability**: clean and uniform images risk overfitting.  
- **Generalization**: models must adapt to unseen data (e.g., electron microscopy, pathological cells).  

---

## Methodology
- **Preprocessing**:  
  - Removed ~2,000 outliers.  
  - Stratified split into train/val/test.  
  - Balanced classes using **upsampling** and advanced augmentations (**RandAugment**, flips, rotations, zooming).  
- **Models**:  
  - Custom CNN (4 convolutional blocks + dense layers).  
  - Transfer learning with **VGG16** and **EfficientNet (B0, B2, V2B2, V2M)**.  
- **Training**:  
  - Gradual fine-tuning (starting from top layers).  
  - Optimizers: **Adam**, **Lion** (achieved best results).  
  - Early stopping to avoid overfitting.  

---

## Results
- Custom CNN trained fast but **failed to generalize** (local accuracy high, but low Codabench performance).  
- VGG16 + augmentation achieved decent results, but was outperformed by EfficientNet.  
- **EfficientNetV2-M with Lion optimizer** reached the **best performance: 85% accuracy** on the hidden test set.  

---

## Discussion
- Transfer learning with EfficientNet proved most effective, balancing depth and efficiency.  
- Fine-tuning deeper layers significantly improved classification.  
- Some augmentation strategies (e.g., RandAugment + CutMix) hindered convergence.  
- More powerful architectures (e.g., VGG19) and hyperparameter exploration could further boost performance.  

---

## Further Details
For the complete methodology, experiments, and results, please refer to the **full project report** here: ([full project report] (https://github.com/annasime1/academic-projects/blob/main/Blood-samples-classification/docs/report.pdf)).



