# Predictive Modeling of Diabetic Patient Readmission

## Overview
This project investigates **hospital readmission prediction for diabetic patients** using statistical learning and machine learning techniques.  
The objective is to classify patients into three categories:
- **<30** → readmission within 30 days  
- **>30** → readmission after 30 days  
- **NO** → no readmission  

The study addresses challenges such as data quality, missingness, class imbalance, and model interpretability, aiming to provide clinically useful insights.

---

## Dataset
- **Source:** 10 years of clinical data (1999–2008) from 130 US hospitals.  
- **Size:** 101,766 encounters → 98,913 after preprocessing.  
- **Features:** 50 attributes including demographics, diagnoses, medications, lab tests, and outcomes.  
- **Target variable:** patient readmission status (`<30`, `>30`, `NO`).  

Key preprocessing steps:
- Removal of inconsistent or biased records (e.g., hospice discharges).  
- Grouping high-cardinality categorical variables into clinically meaningful categories.  
- Handling missing values with imputation strategies (mode, KNN, “Unknown” category).  
- Outlier detection and skewness correction (log transformations).  
- One-hot encoding of categorical features and z-score normalization of numerical features.  

---

## Models and Training
We experimented with several classification algorithms:
- **Logistic Regression** (baseline, interpretable)  
- **Random Forest** (robust, balanced performance)  
- **XGBoost** (flexible, strong minority class recall)  
- **Two-stage Gating Model** (focused on early readmissions `<30`)  

Techniques applied:
- **Class weighting** and **SMOTE** to address imbalance.  
- **5-fold cross-validation** and **bootstrap resampling** for robust evaluation.  
- **Feature importance and SHAP** for interpretability.  

---

## Results
Main findings:
- **Balanced Logistic Regression** → interpretable, stable, good recall on minority classes.  
- **Balanced Random Forest** → best overall F1-score (54.2%), well-balanced across classes.  
- **XGBoost** → highest recall for `<30` (47.7%), but lower global performance.  
- **Gating Model** → strongest sensitivity for `<30` (59.0%), suited for screening workflows.  

Interpretability analysis highlighted **prior admissions, discharge disposition, and payer type** as key predictors of readmission risk.

---

## Report
The full report with detailed methodology, experiments, and discussion is available here: [Read the full report](./docs/report.pdf)

---

## Future Directions
- Explore temporal models (e.g., LSTMs).  
- Validate on external datasets for generalizability.  
- Apply causal inference for better understanding of clinical factors.  
- Integrate into explainable clinical decision support systems.  

