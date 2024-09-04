Bacteriocin Classification Project

This repository contains R scripts used for the classification of potential bacteriocins from wastewater metagenomes using machine learning techniques. These scripts handle feature extraction, feature selection, and the training and evaluation of models. Below is a brief description of each script:

Feature Extraction

all_acc.R – Extracts amino acid composition features.
all_amphipseudo.R – Extracts amphiphilicity pseudo features.
all_dipep.R – Extracts dipeptide composition features.
all_pseudo.R – Extracts pseudo-amino acid features.
all_ss.R – Extracts secondary structure features.
all_distribution.R – Extracts global distribution features.
all_transition.R – Extracts global transition features.
all_composition.R – Extracts global composition features.

Feature Selection

MDG_reduced_features.R – Performs feature extraction and reduction using the Mean Decrease Gini (MDG) method from Random Forest. 
2nd_pearson_correlation.R – Performs Pearson correlation analysis for feature reduction.
rfe_mdg.R – Performs Recursive Feature Elimination (RFE) using Mean Decrease Gini (MDG) for feature selection.
rfe_mdg_ML_RF.R – Combines RFE and MDG for Random Forest models.
rfe_mdg_ML_SVM.R – Combines RFE and MDG for Support Vector Machine (SVM) models.

Machine Learning Models

2.0_ADTree.R – Implements an Alternating Decision Tree (ADTree) model.
ADT_RF.R – Trains a Random Forest (RF) classifier.
ADT_SVM.R – Trains a Support Vector Machine (SVM) classifier
