# Day-35 immune cell composition (ICC35) is associated with overall survival post allogeneic stem cell transplantation
This is a GitHub repository containing main pieces of code for the ICC35 score publication. Below is the description of what each folder contains:

### Delete Data_Preprocessing_Integration_Clustering
Contains the scripts for initial preprocessing, integration & clustering of the flow cytometry data:
<ul>
  <li>Step1.Removing_TopBottom_0.01%Percentiles.R</li>
  <li>Step2.2.cyCondor_SubsettingForVisualisation.R</li>
  <li>Step2.cyCondor_ReadingFiles_Integrating_Clustering.R</li>
  <li>Step3.Assessing_Integration_Efficiency.R</li>
  <li>Step4.cyCondor_ReadingFiles_Integrating_Clustering_CD3+Lymphocytes.R</li>
  <li>Step5.Merging_Cluster_Annotations_CD3-_&_CD3+.R</li>
  <li>Step6.Characterizing_Clusters.R</li>
</ul>
For step 5 - there are also two excel tables:
<ul>
  <li>Step5.Splitting_Clusters.xlsx (which clusters require further specification at higher resolution)</li>
  <li>Step5.Cluster_Identities.xlsx (immune phenotype identities for each cluster)</li>
</ul>

### Survival_Analysis
Contains the scripts for covnerting the immune population counts to fractions and ratios, running the main survival analyses<br></br>
Subfolder <b>Cause-specific proportional hazard models</b> contains scripts for testing every immune parameter in different cause-specific proportional hazard model setups. These python scripts were run on the computational cluster via <i>slurm sbatch</i> scripts.
<ul>
  <li>CauseSpecificHazardModel_CompetingRisks_Relapse_cGVHD_and_aGVHD.py (competing risks are relapse and either aGVHD or cGVHD)</li>
  <li>CauseSpecificHazardModel_CompetingRisks_Relapse_cGVHD_or_aGVHD.py (competing risks are relapse, cGVHD, and aGVHD)</li>
  <li>CauseSpecificHazardModel_GRFS (cause-specific proportional hazard modelling for GVHD/relapse-free survival - competing events are relapse and aGVHD/cGVHD).py</li>
  <li>CauseSpecificHazardModel_NRM.py (cause-specific proportional hazard modelling for non-relapse mortality)</li>
</ul>
Subfolder <b>Univariate & Multivariate Cox proportional hazard (CPH) models</b> contains:
<ol>
  <li>Scripts for running the univariate/ multivariate Cox proportional hazard analyses (for clinical covariates and/or ICC35 score)</li>
  <li>Technical checks on ICC35 score (association with other covariates, selecting penalizer values for multivariate models)</li>
</ol>
The scripts are:
<ul>
  <li>Running_Multivariate_Analysis.ipynb (main Jupyter notebook for running the multivariate CPH models)</li>
  <li>Running_Univariate_Multivariate_Analysis_ClinVars_only.ipynb (Jupyter notebook for running the univariate and multivariate CPH models only for clinical covariates)</li>
  <li>Selecting_Penalizer_Value_for_Multivariate_Analysis.ipynb (Jupyter notebook for selecting the optimal penalizer value for multivariate CPH models with ICC35 score)</li>
  <li>Testing_Association_of_ClinVars_with_ICC35score.ipynb (Jupyter notebook for testing whether there are associations between clinical covariates and ICC35 score)</li>
</ul>

### Visualisations
Contains the scripts for visualisation of the flow cytometry data and survival analysis findings
