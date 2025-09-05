# Attention Classification in Children with ASD using Deep Learning

## Overview
This project focuses on developing a **deep learning algorithm** to classify attention in children with **Autism Spectrum Disorder (ASD)** using gaze tracking data collected during interaction with a **social assistive robot**.  
The aim is to support **therapists** in evaluating attentional patterns during therapy sessions, providing an objective and automated tool for behavioral assessment.  

---

## Dataset
- Data acquired from sessions involving **children with ASD, neurotypical subjects, and therapists**.  
- Each sample includes **gaze features** (azimuth, elevation, eye positions, confidence score, and use of robot) extracted from video recordings.  
- Preprocessing steps included:  
  - **Outlier removal** (frames with wrong robot usage or extreme gaze values).  
  - **Feature selection** (confidence, eyes, azimuth, elevation, use_robot).  
  - **Normalization**: Min–Max scaling vs. **Z-score normalization** (the latter yielded better performance).  
- Dataset split into **training, validation, and test sets**.  

---

## Methodology
- **Model**: Multi-Layer Perceptron (MLP).  
- **Training strategies**:  
  - Online **data augmentation** to reduce overfitting.  
  - **Dropout layers** (rate = 0.3) for regularization.  
  - **Activation functions**: LeakyReLU (to address dying neuron problem).  
  - **Focal Loss** to handle class imbalance.  
  - Upsampling and downsampling to balance classes.  
- Special constraint introduced to improve classification of the **"robot" vs. "nowhere"** gaze directions, which were often confused.  

---

## Results
- Best configuration (Focal Loss + data augmentation in the first layer) achieved an **F1-score of ~47%** on the test set.  
- The model struggled with **underrepresented classes** and with distinguishing **"nowhere"** from **"robot"** gaze directions.  
- Feature importance analysis showed that **confidence and gaze features** significantly influenced performance.  

---

## Limitations
- Small dataset with limited variability.  
- Noise and inaccuracies in gaze extraction pipeline.  
- Overlap between "nowhere" and "robot" classes caused frequent misclassifications.  

---

## Future Work
- Acquire additional data to improve model generalization.  
- Enhance the **gaze extraction algorithm** for higher accuracy.  
- Explore a **two-stage classification approach**:  
  1. Binary classifier to distinguish **target vs. nowhere**.  
  2. Target-specific classifier for finer gaze categorization.  

---

## Further Details
For additional insights and methodology, please refer to the **video** included in this repository.

