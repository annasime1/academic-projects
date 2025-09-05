# Lung Nodule Classification in CT Images Using Deep Learning

## Overview
This project explores the use of **deep learning for lung cancer diagnosis** through CT scans.  
The objective is to classify lung nodules as **benign vs. malignant (binary classification)** and to predict **malignancy scores (1–5, multi-class classification)**.  

Two types of inputs were compared:  
- **Full-slice CT images**, providing anatomical context.  
- **Zoomed nodule images**, focusing on the lesion.  

The project leverages **custom CNNs** and **transfer learning with EfficientNetB2**, enhanced with explainability techniques such as **Grad-CAM** and **LIME**.

---

## Dataset
- CT scans from ~2,400 patients.  
- Each sample includes both a **full-slice CT image** and a **zoomed nodule crop**.  
- Labels: **malignancy score 1–5**, used for both binary and multi-class classification.  
- Split: **70% training, 15% validation, 15% test**.  

### Challenges
- **Data scarcity** (limited dataset size).  
- **Class imbalance** (e.g., benign nodules overrepresented).  
- **Generalization**: risk of overfitting due to limited and imbalanced data.  

---

## Methodology
- **Preprocessing**: normalization, sigmoidal contrast enhancement, center cropping, resizing, and class balancing via augmentation.  
- **Models**:  
  - Custom CNN (baseline).  
  - Pre-trained MobileNet, ResNet50, EfficientNetB2.  
- **Training**: gradual fine-tuning, Adam optimizer with weight decay, early stopping, and dynamic augmentation.  
- **Explainability**: Grad-CAM, LIME, and image retrieval analysis to interpret model predictions.  

---

## Results
- **Binary classification (nodules, EfficientNetB2 fine-tuned)**: **F1-score ~81%**, best performance overall.  
- **Binary classification (full-slice)**: F1-score ~70%.  
- **Multi-class classification (nodules)**: F1-score ~45%.  
- **Multi-class classification (full-slice)**: F1-score ~39%.  

Models trained on **zoomed nodules consistently outperformed** those trained on full-slice images.  

---

## Discussion
- Zoomed nodule images reduced background noise and improved classification accuracy.  
- Binary classification proved more reliable than multi-class.  
- Explainability confirmed that nodule-based models focused on the correct regions, while full-slice models often misinterpreted background structures.  
- A hybrid approach combining **full-slice + nodule models** could enhance robustness in clinical applications.  

---

## Future Development
- Implement cross-validation and hyperparameter tuning.  
- Integrate multimodal data (clinical + imaging).  
- Explore ensemble methods combining full-slice and nodule models.  

---

## Further Details
For the complete methodology, experiments, and results, please refer to the **[full project report](https://github.com/annasime1/academic-projects/blob/main/LungCancer-Classification/docs/report.pdf)**.

