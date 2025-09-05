# Semantic Segmentation of Martian Terrain

## Overview
This project applies **deep learning techniques** to the semantic segmentation of Martian terrain images.  
The task consists of classifying each pixel into five categories: **Background, Soil, Bedrock, Sand, Big Rock**.  

The solution is based on the **U-Net architecture**, enhanced with **Atrous Spatial Pyramid Pooling (ASPP)** and **Attention Blocks** to improve multi-scale feature extraction and segmentation accuracy.

---

## Methodology
- **Dataset**: 64×128 grayscale images with pixel-wise masks.  
- **Preprocessing**: outlier removal, augmentation (flips, contrast/brightness changes, cropping, resizing).  
- **Architecture**: U-Net variations (2–4 blocks) with residual connections, ASPP, and attention mechanisms.  
- **Training**: AdamW optimizer, learning rate schedulers, fine-tuning of decoder and bottleneck.  
- **Loss functions**: categorical cross-entropy, weighted loss, focal loss.  

---

## Results
- Deeper architectures captured more complex terrain features.  
- ASPP and attention blocks significantly improved segmentation quality.  
- Fine-tuning the bottleneck and decoder further boosted performance.  
- Some augmentations (e.g., brightness/contrast) unexpectedly reduced accuracy.  

---

## Further Details
For the complete methodology, experiments, and results, please refer to the **full project report** here: [full project report](https://github.com/annasime1/academic-projects/blob/main/Mars-terrain-segmentation/docs/report.pdf).
