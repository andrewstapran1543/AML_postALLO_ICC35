# importing the packages:
import os
import math
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from tqdm.notebook import tqdm

from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.impute import SimpleImputer

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats.multitest import multipletests
from lifelines.statistics import logrank_test, multivariate_logrank_test

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.table import Table
import matplotlib.transforms as mtransforms

import scipy.stats
from scipy.stats import norm, zscore, fisher_exact, chi2_contingency, mannwhitneyu




# reading the main parameters:
print('reading the main parameters:')
patient_subset = str(sys.argv[1])
print(f'Patient subset: {patient_subset}')
severe_aGVHD_cats = [int(element) for element in str(sys.argv[2]).split(',')]
print(f'severe_aGVHD_cats: {sys.argv[2]}')
cutoff = int(sys.argv[3])
print(f'Censoring cutoff: {sys.argv[3]}')
print('---------------')




# reading the main dataframe:
print('reading the main dataframe:')
# path to the full dataset excel with clinical metadata and flow cytometry parameters for the patients:
full_dataset_clinmetadata_flowcytometrydata_path = ''
d35_input_full = pd.read_excel(full_dataset_clinmetadata_flowcytometrydata_path)
d35_input_full_copy = d35_input_full[(d35_input_full['OC_d35_PassCrit_FCAnalysis'] == 'Yes') | (d35_input_full['VC_d35_PassCrit_FCAnalysis'] == 'Yes')]

PERCENT_FEATURES = [c for c in d35_input_full_copy.columns if ' of ' in c]
RATIO_FEATURES   = [c for c in d35_input_full_copy.columns if '-to-' in c]
d35_input_full_copy[PERCENT_FEATURES] = np.arcsin(np.sqrt(np.clip(d35_input_full_copy[PERCENT_FEATURES] / 100, 0, 1)))

d35_input_full_copy['Sex_Mismatch'] = np.where(d35_input_full_copy['Sex_Recipient'] == d35_input_full_copy['Sex_Donor'], 0, 1)
d35_input_full_copy['Sex_Recipient'] = np.where(d35_input_full_copy['Sex_Recipient'] == 'Male', 0, 1) # Male = 0, Female = 1
d35_input_full_copy['CMV_Status_Mismatch'] = np.where(d35_input_full_copy['CMV_Status_Recipient'] == d35_input_full_copy['CMV_Status_Donor'], 0, 1)
d35_input_full_copy['CMV_Status_Recipient'] = np.where(d35_input_full_copy['CMV_Status_Recipient'] == 'neg', 0, 1) # neg = 0, pos = 1
d35_input_full_copy['MHC_Mismatch_between_Donor&Recipient'] = np.where(d35_input_full_copy['MHC_Mismatch_between_Donor&Recipient'] == 'nein', 0, 1) # no = 0, yes = 1

OUTCOME_COLS = {
    "RFS": {"time": 'Time_to_Relapse_from_TPL', "event": 'Relapse_Present_1_NotObserved_0'},
    "OS":  {"time": 'Time_to_Death_from_TPL',  "event": 'Death_Present_1_NotObserved_0'},
    "aGVHD": {"time": "Time_to_aGVHD", "event": "aGVHD_grading_by_Katja"},
    "cGVHD": {"time": "Time_to_cGVHD", "event": "cGVHD_grading_by_Katja"},
}

d35_input_full_copy[OUTCOME_COLS['RFS']['event']] = np.where(d35_input_full_copy[OUTCOME_COLS['RFS']['time']] >= cutoff, 0, d35_input_full_copy[OUTCOME_COLS['RFS']['event']])
d35_input_full_copy[OUTCOME_COLS['RFS']['time']] = np.where(d35_input_full_copy[OUTCOME_COLS['RFS']['time']] >= cutoff, cutoff, d35_input_full_copy[OUTCOME_COLS['RFS']['time']])

d35_input_full_copy[OUTCOME_COLS['OS']['event']] = np.where(d35_input_full_copy[OUTCOME_COLS['OS']['time']] >= cutoff, 0, d35_input_full_copy[OUTCOME_COLS['OS']['event']])
d35_input_full_copy[OUTCOME_COLS['OS']['time']] = np.where(d35_input_full_copy[OUTCOME_COLS['OS']['time']] >= cutoff, cutoff, d35_input_full_copy[OUTCOME_COLS['OS']['time']])

d35_input_full_copy[OUTCOME_COLS['aGVHD']['event']] = np.where(d35_input_full_copy[OUTCOME_COLS['aGVHD']['time']] >= cutoff, 0, d35_input_full_copy[OUTCOME_COLS['aGVHD']['event']])
d35_input_full_copy[OUTCOME_COLS['aGVHD']['time']] = np.where(d35_input_full_copy[OUTCOME_COLS['aGVHD']['time']] >= cutoff, cutoff, d35_input_full_copy[OUTCOME_COLS['aGVHD']['time']])

d35_input_full_copy[OUTCOME_COLS['cGVHD']['event']] = np.where(d35_input_full_copy[OUTCOME_COLS['cGVHD']['time']] >= cutoff, 0, d35_input_full_copy[OUTCOME_COLS['cGVHD']['event']])
d35_input_full_copy[OUTCOME_COLS['cGVHD']['time']] = np.where(d35_input_full_copy[OUTCOME_COLS['cGVHD']['time']] >= cutoff, cutoff, d35_input_full_copy[OUTCOME_COLS['cGVHD']['time']])


d35_input_full_copy['Relapse_Present_1_NotObserved_0'] = d35_input_full_copy['Relapse_Present_1_NotObserved_0'].astype("Int64")
d35_input_full_copy['aGVHD_grading_by_Katja'] = d35_input_full_copy['aGVHD_grading_by_Katja'].astype("Int64")
d35_input_full_copy['cGVHD_grading_by_Katja'] = d35_input_full_copy['cGVHD_grading_by_Katja'].astype("Int64")
d35_input_full_copy['FC_Sampling_Days_After_Allo_d35'] = round(d35_input_full_copy['FC_Sampling_Days_After_Allo_d35'], 0)
d35_input_full_copy['FC_Sampling_Days_After_Allo_d35'] = d35_input_full_copy['FC_Sampling_Days_After_Allo_d35'].astype("Int64")
d35_input_full_copy['MAX_FOLLOW_UP'] = d35_input_full_copy[['Time_to_aGVHD', 'Time_to_cGVHD',
                                                            'Time_to_Relapse_from_TPL','Time_to_Death_from_TPL']].max(axis = 1)


d35_input_full_copy['Time_to_Death_from_TPL'] -= d35_input_full_copy['FC_Sampling_Days_After_Allo_d35']
d35_input_full_copy['Time_to_Relapse_from_TPL'] -= d35_input_full_copy['FC_Sampling_Days_After_Allo_d35']
d35_input_full_copy['Time_to_aGVHD'] -= d35_input_full_copy['FC_Sampling_Days_After_Allo_d35']
# only 65 patients (40 in VC and 25 in OC)
d35_input_full_copy['Time_to_cGVHD'] -= d35_input_full_copy['FC_Sampling_Days_After_Allo_d35']
d35_input_full_copy['MAX_FOLLOW_UP'] -= d35_input_full_copy['FC_Sampling_Days_After_Allo_d35']
print('---------------')




# setting the main functions:
print('setting the main functions:')
def candidate_time(t, e):
    INF = 10**9
    # e: 0/1 event indicator
    if (e == 0) or pd.isna(t):
        return INF
    if t < 0:
        return INF
    return float(t)

def _p_two_sided_from_z(z: float) -> float:
    # Normal two-sided p without scipy: p = erfc(|z|/sqrt(2))
    return math.erfc(abs(z) / math.sqrt(2.0))

def meta_fixed_from_betas(feature, betas: List[float], ses: List[float]):
    """Inverse-variance fixed-effect meta-analysis over k cohorts."""
    if len(betas) != len(ses) or len(betas) < 2:
        raise ValueError("Provide betas & ses of the same length (k>=2).")
    weights = [1.0 / (s**2) for s in ses]
    w_sum = sum(weights)
    beta_hat = sum(w * b for w, b in zip(weights, betas)) / w_sum
    se_hat = math.sqrt(1.0 / w_sum)
    z_meta = beta_hat / se_hat
    p_meta = _p_two_sided_from_z(z_meta)

    # Heterogeneity
    Q = sum(w * (b - beta_hat) ** 2 for w, b in zip(weights, betas))
    k = len(betas)
    I2 = 0.0
    if Q > 0 and k > 1:
        I2 = max(0.0, (Q - (k - 1)) / Q) * 100.0

    # Back-transform to HR scale
    try:
        hr_pooled = math.exp(beta_hat)
    except:
        hr_pooled = "UnableToCalculate"

    try:
        ci_low = math.exp(beta_hat - 1.96 * se_hat)
    except:
        ci_low = "UnableToCalculate"
    
    try:
        ci_high = math.exp(beta_hat + 1.96 * se_hat)
    except:
        ci_high = "UnableToCalculate"


    return {
        "Marker": feature,
        "beta_pooled": beta_hat,
        "se_pooled": se_hat,
        "z_meta": z_meta,
        "p_meta": p_meta,
        "HR_pooled": hr_pooled,
        "HR_CI95": (ci_low, ci_high),
        "Q": Q,
        "I2_percent": I2,
    }

def compet_risk_converter_median_splitter_NRM(input_df, treat_as_continuous):
    INF = 10**9
    df = input_df.copy()

    df["aGVHD_pre35"] = ((df['severe_aGVHD_simplified'] == 1) & (df[OUTCOME_COLS['aGVHD']['time']] < 0)).astype(int)
    
    # candidate post-d35 times (only if event truly happened)
    df["t_NRM_post"] = [candidate_time(t, e) for t, e in zip(df[OUTCOME_COLS['OS']['time']], df['NonRelapse_Mortality_Present_1_Unobserved_0'])]
    df["t_Relapse_post"] = [candidate_time(t, e) for t, e in zip(df[OUTCOME_COLS['RFS']['time']], df[OUTCOME_COLS['RFS']['event']])]

    # first event time after d35 (capped by follow-up)
    df["t_first"] = df[["t_NRM_post", "t_Relapse_post"]].min(axis=1)
    df["t_first"] = np.minimum(df["t_first"], df["MAX_FOLLOW_UP"])

    # event indicators for cause-specific Cox
    df["NRM_first"] = ((df["t_NRM_post"] == df["t_first"]) & (df["t_first"] <= df["MAX_FOLLOW_UP"])).astype(int)
    df["relapse_first"]    = ((df["t_Relapse_post"] == df["t_first"]) & (df["t_first"] <= df["MAX_FOLLOW_UP"])).astype(int)

    d35_cohort_OC = df[df['Cohort'] == 'OC']
    # print('OC', len(d35_cohort_OC))
    d35_cohort_VC = df[df['Cohort'] == 'VC']
    # print('VC', len(d35_cohort_VC))

    if treat_as_continuous:
        for feature in PERCENT_FEATURES + RATIO_FEATURES:
            d35_cohort_OC[feature] = zscore(d35_cohort_OC[feature], nan_policy = 'omit')
            d35_cohort_VC[feature] = zscore(d35_cohort_VC[feature], nan_policy = 'omit')
    else:
        split = 2
        for feature in PERCENT_FEATURES + RATIO_FEATURES:
            d35_cohort_OC[feature] = pd.qcut(d35_cohort_OC[feature], q=split, labels=list(range(0,split)))
            d35_cohort_VC[feature] = pd.qcut(d35_cohort_VC[feature], q=split, labels=list(range(0,split)))
    
    return d35_cohort_OC, d35_cohort_VC

def HR_risk_plotter(input_df, to_visualize, input_col, data_type, analysis_type, save_path):
    used_pvalue = None
    df = None  

    for p_value_cutoff in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:
        df = input_df[input_df["max_pvalue_in_2cohorts"] < p_value_cutoff]
        if df.empty:
            continue   # try next threshold
        used_pvalue = p_value_cutoff
        break

    if used_pvalue is None:
        print('No parameters with matching direction & p-values below p-value thresholds of 0.05, 0.1, 0.2, 0.3, 0.4, 0.5')
    else:
        df['HR_OC'] = np.exp(df['beta_coef_OC'])
        df['HR_VC'] = np.exp(df['beta_coef_VC'])

        df['beta_coef_OC_lower'] = df['beta_coef_OC'] - 1.96*df['SE_OC']
        df['HR_OC_lower'] = np.exp(df['beta_coef_OC_lower'])
        df['beta_coef_OC_upper'] = df['beta_coef_OC'] + 1.96*df['SE_OC']
        df['HR_OC_upper'] = np.exp(df['beta_coef_OC_upper'])
        df['beta_coef_VC_lower'] = df['beta_coef_VC'] - 1.96*df['SE_VC']
        df['HR_VC_lower'] = np.exp(df['beta_coef_VC_lower'])
        df['beta_coef_VC_upper'] = df['beta_coef_VC'] + 1.96*df['SE_VC']
        df['HR_VC_upper'] = np.exp(df['beta_coef_VC_upper'])

        n = len(df)
        y = np.arange(n)
        if n <= 5:
            fig, ax = plt.subplots(figsize=(10, 1*n))
        else:
            fig, ax = plt.subplots(figsize=(10, 0.6*n))
        offset = 0.15

        if to_visualize == 'beta': 
            # --- Cohort 1 (blue) ---
            oc = ax.errorbar(
                df["beta_coef_OC"], y - offset,
                xerr=[df["beta_coef_OC"] - df["beta_coef_OC_lower"], df["beta_coef_OC_upper"] - df["beta_coef_OC"]],
                fmt="s", markersize=6,
                color="tab:blue", ecolor="tab:blue",
                capsize=3,
                label=f"Cohort 1\n(2006-2018)\n{list(df['total_patients_OC'].unique())[0]} patients"         # <--- label here
            )

            # --- Cohort 2 (red) ---
            vc = ax.errorbar(
                df["beta_coef_VC"], y + offset,
                xerr=[df["beta_coef_VC"] - df["beta_coef_VC_lower"], df["beta_coef_VC_upper"] - df["beta_coef_VC"]],
                fmt="s", markersize=6,
                color="tab:red", ecolor="tab:red",
                capsize=3,
                label=f"Cohort 2\n(2019 - 2023)\n{list(df['total_patients_VC'].unique())[0]} patients"         # <--- and here
            )
            ax.axvline(0, color="gray", linestyle="--", linewidth=1)  # no label

        elif to_visualize == 'HR':
            # --- Cohort 1 (blue) ---
            oc = ax.errorbar(
                df["HR_OC"], y - offset,
                xerr=[df["HR_OC"] - df["HR_OC_lower"], df["HR_OC_upper"] - df["HR_OC"]],
                fmt="s", markersize=6,
                color="tab:blue", ecolor="tab:blue",
                capsize=3,
                label=f"Cohort 1\n(2006-2018)\n{list(df['total_patients_OC'].unique())[0]} patients"         # <--- label here
            )

            # --- Cohort 2 (red) ---
            vc = ax.errorbar(
                df["HR_VC"], y + offset,
                xerr=[df["HR_VC"] - df["HR_VC_lower"], df["HR_VC_upper"] - df["HR_VC"]],
                fmt="s", markersize=6,
                color="tab:red", ecolor="tab:red",
                capsize=3,
                label=f"Cohort 2\n(2019 - 2023)\n{list(df['total_patients_VC'].unique())[0]} patients"         # <--- and here
            )
            ax.axvline(1, color="gray", linestyle="--", linewidth=1)  # no label

        ax.set_yticks(y)
        ax.set_yticklabels(df["Marker"])
        ax.set_ylim(-0.5, n - 0.5)
        ax.invert_yaxis()

        ax.legend(bbox_to_anchor=(1,1), framealpha = 1.0)                  # <--- no args, uses the labels above
        if to_visualize == 'beta': 
            ax.set_xlabel("log(HR) (95% CI)")
        elif to_visualize == 'HR': 
            ax.set_xlabel("HR (95% CI)")

        if input_col == 'NRM_first':
            variable = 'NRM'
        if input_col == 'relapse_first':
            variable = 'RFS'
        
        title = f"First event predicting: {variable}\nVariables are: {data_type}\nCompetitive risk model setup is: {analysis_type}\np-value cutoff used: {used_pvalue}"
        plt.title(title)
        plt.tight_layout()
        save_path += f'_{used_pvalue}.png'
        plt.savefig(save_path, bbox_inches = 'tight')
        plt.show()

def runner_coxph(OC_intput, VC_input, input_col, treat_as_continuous):
    cohort_list = []
    marker_list = []
    beta_list = []
    se_list = []
    z_list = []
    p_list = []
    how_many_events = []
    total_patients = []

    for cohort in ['OC', 'VC']:
        if cohort == 'OC':
            dataset = OC_intput.copy()
        else:
            dataset = VC_input.copy()

        for marker in PERCENT_FEATURES + RATIO_FEATURES:
            df_model = dataset[[marker, 't_first', input_col, "aGVHD_pre35"]]
            formula_cph = "marker + aGVHD_pre35"
            
            df_model = df_model.rename(columns={marker: 'marker'})
            if df_model['marker'].isna().any():
                df_model = df_model.dropna(subset = ['marker'])
                # print(cohort, marker, len(dataset), len(df_model))
            else:
                df_model = df_model
        
            try:
                cph = CoxPHFitter()
                cph.fit(df_model, duration_col="t_first", event_col=input_col,formula=formula_cph)
                if treat_as_continuous:
                    s = cph.summary.loc['marker']
                else:
                    s = cph.summary.loc['marker[T.1]']

                total_patients.append(len(df_model))
                how_many_events.append(len(df_model[df_model[input_col] == 1]))
                cohort_list.append(cohort)
                marker_list.append(marker)
                beta_list.append(float(s["coef"]))
                se_list.append(float(s["se(coef)"]))
                z_list.append(float(s["z"]))
                p_list.append(float(s["p"]))

            except:
                print(f'{cohort}, {marker}: Unable to run the cox model\nTotal: {len(df_model)}, With event: {len(df_model[df_model[input_col] == 1])}' )

    output = pd.DataFrame([cohort_list,marker_list,beta_list,se_list,z_list,p_list,how_many_events,total_patients]).T
    output.columns = ['Cohort', 'Marker', 'beta_coef', 'SE', 'z', 'p', 'how_many_events', 'total_patients']

    output_OC = output[output['Cohort'] == 'OC']
    pvals = output_OC["p"].values
    output_OC["q"] = multipletests(pvals, method="fdr_bh")[1]
    output_OC.drop(['Cohort'], axis = 1, inplace = True)
    output_OC.columns = [f'{element}_OC' if element not in ['Marker'] else element for element in output_OC.columns]


    output_VC = output[output['Cohort'] == 'VC']
    pvals = output_VC["p"].values
    output_VC["q"] = multipletests(pvals, method="fdr_bh")[1]
    output_VC.drop(['Cohort'], axis = 1, inplace = True)
    output_VC.columns = [f'{element}_VC' if element not in ['Marker'] else element for element in output_VC.columns]

    common_markers = list(set(list(output_OC['Marker'])) & set(list(output_VC['Marker'])))
    total_markers = PERCENT_FEATURES + RATIO_FEATURES
    if len(common_markers) < len(total_markers):
        difference = len(total_markers) - len(common_markers)
        print(f"{difference} out of {len(total_markers)} missing in the final output")
    output_OC = output_OC[output_OC['Marker'].isin(common_markers)]
    output_VC = output_VC[output_VC['Marker'].isin(common_markers)]
    merged_output = pd.merge(output_OC, output_VC, on='Marker', how='outer')

    merged_output['MATCHING_DIRECTION'] = np.where(np.sign(merged_output['beta_coef_OC']) == np.sign(merged_output['beta_coef_VC']), 'MATCH', 'NON-MATCH')
    merged_output['EFFECT'] = np.where(
        merged_output['MATCHING_DIRECTION'] == 'MATCH', np.where(
            np.sign(merged_output['beta_coef_OC']) < 0, 'protective', 'non-protective'
        ), 
        merged_output['MATCHING_DIRECTION']
    )

    merged_output['max_pvalue_in_2cohorts'] = merged_output[['p_OC', 'p_VC']].max(axis = 1)
    
    merged_output['beta_coef_OC_ABS'] = merged_output['beta_coef_OC'].abs()
    merged_output['beta_coef_VC_ABS'] = merged_output['beta_coef_VC'].abs()
    merged_output['min_beta_in_2cohorts_ABS'] = merged_output[['beta_coef_OC_ABS', 'beta_coef_VC_ABS']].min(axis = 1)
    merged_output['min_beta_in_2cohorts'] = np.where(merged_output['EFFECT'] == 'protective', merged_output['min_beta_in_2cohorts_ABS'] * (-1), merged_output['min_beta_in_2cohorts_ABS'])
    
    merged_output = merged_output.sort_values(by = ['max_pvalue_in_2cohorts','min_beta_in_2cohorts_ABS'], ascending = [True, False])
    merged_output.drop(['beta_coef_OC_ABS', 'beta_coef_VC_ABS', 'min_beta_in_2cohorts_ABS'], axis = 1, inplace = True)

    meta_output = []
    for feature in merged_output['Marker']:
        betas = [
            float(merged_output.loc[merged_output['Marker'] == feature, 'beta_coef_OC']),
            float(merged_output.loc[merged_output['Marker'] == feature, 'beta_coef_VC'])
        ]
        ses = [
            float(merged_output.loc[merged_output['Marker'] == feature, 'SE_OC']),
            float(merged_output.loc[merged_output['Marker'] == feature, 'SE_VC'])
        ]
        meta_output.append(meta_fixed_from_betas(feature, betas, ses))
    
    meta_output = pd.DataFrame(meta_output)
    merged_output = pd.merge(merged_output, meta_output, on = 'Marker')

    # specify the base_dir for recording the results (folder name should specify # of days used to horizon the clinical follow-up)
    base_dir = ''
        
    if len(severe_aGVHD_cats) == 1:
        aGVHD_part = 'aGVHD' + str(severe_aGVHD_cats[0])
    elif len(severe_aGVHD_cats) > 1:
        aGVHD_part = 'aGVHD' + ''.join([str(element) for element in severe_aGVHD_cats])
    
    if treat_as_continuous:
        add_on = 'CONTINUOUS'
    else:
        add_on = 'SPLIT_BY_MEDIAN'

    analysis_type = f'NRM_Relapse_{aGVHD_part}'
    save_path_table = f'{base_dir}{patient_subset}/{analysis_type}/{input_col}_{add_on}.xlsx'
    merged_output.to_excel(save_path_table, index = False)

    merged_output = pd.read_excel(save_path_table)
    save_path_plot = f'{base_dir}{patient_subset}/{analysis_type}/{input_col}_{add_on}'
    HR_risk_plotter(merged_output[merged_output['MATCHING_DIRECTION'] == 'MATCH'], 'beta', input_col, add_on, analysis_type, save_path_plot)
print('---------------')




# running the main code:
print('running the main code:')
if patient_subset == 'all_patients':
    d35_input_scores = d35_input_full_copy.copy()
elif patient_subset == 'CMVRecPos':
    d35_input_scores = d35_input_full_copy[d35_input_full_copy['CMV_Status_Recipient'] == 1]
elif patient_subset == 'CMVRecNeg':
    d35_input_scores = d35_input_full_copy[d35_input_full_copy['CMV_Status_Recipient'] == 0]


d35_input_scores['severe_aGVHD_simplified'] = np.where(d35_input_scores['aGVHD_grading_by_Katja'].isin(severe_aGVHD_cats), 1, 0)
for treat_as_continuous in [True, False]:
    print('treat_as_continuous', treat_as_continuous)
    d35_cohort_OC, d35_cohort_VC = compet_risk_converter_median_splitter_NRM(d35_input_scores, treat_as_continuous)
    
    for input_col in ['NRM_first', 'relapse_first']:
        runner_coxph(d35_cohort_OC, d35_cohort_VC, input_col, treat_as_continuous)
print('Complete')
print('---------------')
        