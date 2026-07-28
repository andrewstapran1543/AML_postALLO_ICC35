# Day-35 immune cell composition (ICC35) is associated with overall survival post allogeneic stem cell transplantation
This is a GitHub repository containing main pieces of code for the ICC35 score publication. Below is the description for what each file does:

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
