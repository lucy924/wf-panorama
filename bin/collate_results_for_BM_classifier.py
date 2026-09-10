#!/usr/bin/env python3

import argparse
import pandas as pd
from shared_functions import get_BM_TYPE_FULL, get_full_SCORING_TYPE, variant_prep, filter_to_variant_type, preclin_stage_panel_result_header, BIOMARKER_ID, BIOMARKER_NAME, SCORING_TYPE, BIOMARKER_TYPE, RESULT_OPTIONS, RESULT, NOTES

ULTRA_LOW = "very lowly expressed/not expressed"
LOW = "lowly expressed"
MODERATE = "moderately expressed/neutral"
HIGH = "highly expressed"
ULTRA_HIGH = "very highly expressed"

######### Set thresholds for expression #########
# TODO: WARNING: These are arbitary values set during code development. Requires review.
# REVIEW
methylation_threshold = 0.8  # if region is >=80% methylated, then consider it methylated
very_low_expression_threshold = 0.8  # if region is >=80% methylated, then consider it no/low expressed
low_expression_threshold = 0.65  # if region is 65%-80% methylated, then consider it low expressed
# if region is 40%-65% methylated, then consider it moderately expressed
high_expression_threshold = 0.40  # if region is 10%-40% methylated, then consider it highly expressed
very_high_expression_threshold = 0.10  # if region is <=10% methylated, then consider it very highly expressed
######### ######### ######### ######### #########


def get_header(metadata_df, identifier):
        # Get header from panel metadata file
        metadata_df_header = list(metadata_df.columns)
        mol_group_header = None
        for col in metadata_df_header:
            if col.startswith(identifier):
                mol_group_header = col
                return mol_group_header
        if mol_group_header is None:
            raise ValueError("Molecular Group header not found in panel metadata file")
        
        
def get_HR_result(marker_components, results_dict):
    """HR = Hormone Receptor"""
    
    # marker componenents for HR+ are [ESR1, ESR2, PGR]
    marker_score = 0
    for component in marker_components:
        if component in results_dict.keys():
            # check if result is an expression-based marker
            if results_dict[component][BIOMARKER_TYPE] not in ['mod - promoter region', 'expression']:
                print(f"Component: {component} is not an expression-based marker, skipping...")
                continue
            
            # check expression level
            exp_value = results_dict[component][RESULT]
            print(f"Component: {component}, Expression level: {exp_value}")
            
            if exp_value <= very_high_expression_threshold:
                exp_level = ULTRA_HIGH
                marker_score += 2
            elif exp_value <= high_expression_threshold:
                exp_level = HIGH
                marker_score += 1
            elif exp_value <= low_expression_threshold:
                exp_level = MODERATE
                marker_score += 0
            elif exp_value <= very_low_expression_threshold:
                exp_level = LOW
                marker_score -= 1
            elif exp_value > very_low_expression_threshold:
                exp_level = ULTRA_LOW
                marker_score -= 2
            else:
                raise ValueError(f"Component: {component} expression level is unknown, something went wrong.")
            
            print(f"Component: {component} is {exp_level}")
    
    if marker_score >= 2:
        composite_result = "Positive (High)"
    elif marker_score == 1: 
        composite_result = "Positive"
    elif marker_score == 0:
        composite_result = "Neutral"
    elif marker_score == -1:
        composite_result = "Negative"
    elif marker_score <= -2:
        composite_result = "Negative (Low)"
    else:
        raise ValueError(f"Marker score: {marker_score}. Something went wrong.")
        
    return composite_result


def get_gene_exp_level(results_dict, component):
    
    # marker components for HER2-low and HER2-ultralow is ERBB2

    if component in results_dict.keys():
        # check if result is an expression-based marker
        if results_dict[component][BIOMARKER_TYPE] not in ['mod - promoter region', 'expression']:
            print(f"Component: {component} is not an expression-based marker, skipping...")
            return None
        
        # check expression level
        exp_value = results_dict[component][RESULT]
        print(f"Component: {component}, Expression level: {exp_value}")
        
        if exp_value >= very_low_expression_threshold:
            # methylation is 80%
            # expression is ultra-low
            exp_level = ULTRA_LOW
            return exp_level, exp_value
        elif exp_value >= low_expression_threshold:
            # methylation is 65%-80%
            # expression is low
            exp_level = LOW
            return exp_level, exp_value
        elif component == "ERBB2":
            # only these are relevant for HER2
            # methylation is <65%
            # expression is moderate to high
            # does not meet criteria for HER2-low or HER2-ultralow
            exp_level = "moderately to highly expressed"
            return exp_level, exp_value
        elif exp_value >= high_expression_threshold:
            # methylation is 40%-65%
            # expression is moderate
            exp_level = MODERATE
            return exp_level, exp_value
        elif exp_value >= very_high_expression_threshold:
            # methylation is 10%-40%
            # expression is high
            exp_level = HIGH
            return exp_level, exp_value
        elif exp_value < very_high_expression_threshold:
            # methylation is <10%
            # expression is ultra-high
            exp_level = ULTRA_HIGH
            return exp_level, exp_value
        else:
            raise ValueError(f"Component: {component} expression level is unknown, something went wrong.")
        
    return None, None


def get_HRD_results():
    """Implement external module/program/tool for HRD classifier"""
    print("HRD is not yet implemented.")
    
    return


def get_composite_biomarkers(metadata_df, BM_TYPE_FULL, results_df):
    # Get composite biomarkers from panel metadata file
    
    # identify column  containing the composition of the markers
    COMPOSITE_COMPONENT_HEADER = get_header(metadata_df, 'Combined Marker')
    
    # select rows of biomarker type: composite
    composite_biomarkers_df = filter_to_variant_type(metadata_df, 'composite', BM_TYPE_FULL)
    
    # convert results_df to dict by biomarker name
    results_dict = results_df.set_index(BIOMARKER_ID).to_dict(orient='index')
    
    # get her2 result
    HER2_level, HER2_value = get_gene_exp_level(results_dict, component="ERBB2")
    
    
    for i, marker in composite_biomarkers_df.iterrows():
        marker_id = marker[BIOMARKER_ID]
        marker_name = marker[BIOMARKER_NAME]
        marker_components = marker[COMPOSITE_COMPONENT_HEADER].split('/')
        print(f"Composite biomarker: {marker_id} ({marker_name}) is composed of: {marker_components}")
        
        ######### Set up code for 301, HR+ #########
        if marker_name == "HR+":
            composite_result = get_HR_result(marker_components, results_dict)
            print(f"Composite biomarker: {marker_id} ({marker_name}) result is: {composite_result}")
            # add composites to results_df
            comp_df = pd.DataFrame({BIOMARKER_ID: marker_id, BIOMARKER_NAME: marker_name, SCORING_TYPE: "categorical", BIOMARKER_TYPE: "composite", RESULT_OPTIONS: "Positive (High), Positive, Neutral, Negative, Negative (Low)", RESULT: composite_result}, index=([0]))
            results_df = pd.concat([results_df, comp_df])

        ######### Set up code for 302, 303 HR+ #########
        elif marker_name == "HER2-low":
            if HER2_level == LOW:
                comp_df = pd.DataFrame({BIOMARKER_ID: marker_id, BIOMARKER_NAME: marker_name, SCORING_TYPE: "categorical", BIOMARKER_TYPE: "composite", RESULT_OPTIONS: f"{LOW}, {ULTRA_LOW}", RESULT: HER2_value}, index=([0]))
                results_df = pd.concat([results_df, comp_df])
        
        elif marker_name == "HER2-ultralow":
            if HER2_level == ULTRA_LOW:
                comp_df = pd.DataFrame({BIOMARKER_ID: marker_id, BIOMARKER_NAME: marker_name, SCORING_TYPE: "categorical", BIOMARKER_TYPE: "composite", RESULT_OPTIONS: f"{LOW}, {ULTRA_LOW}", RESULT: HER2_value}, index=([0]))
                results_df = pd.concat([results_df, comp_df])

        ######### Set up code for 304, TNBC #########
        elif marker_name == "TNBC":
            # If any of the markers are high, then breast cancer is not triple negative.
            # Therefore, start with TNBC = True, and if any of the markers are high, then set TNBC = False
            TNBC = True
            for component in marker_components:
                exp_level, exp_value = get_gene_exp_level(results_dict, component)
                if exp_level in [HIGH, ULTRA_HIGH]:
                    TNBC = False
                    # TODO: note that moderate expression is currently not enough to rule out TNBC according to this code, will need review by expert (if moderate expression should rule out TNBC, then add MODERATE to this list)
                    # REVIEW
            if TNBC:
                comp_df = pd.DataFrame({BIOMARKER_ID: marker_id, BIOMARKER_NAME: marker_name, SCORING_TYPE: "categorical", BIOMARKER_TYPE: "composite", RESULT_OPTIONS: "True, False", RESULT: True}, index=([0]))
                results_df = pd.concat([results_df, comp_df])
            else:
                comp_df = pd.DataFrame({BIOMARKER_ID: marker_id, BIOMARKER_NAME: marker_name, SCORING_TYPE: "categorical", BIOMARKER_TYPE: "composite", RESULT_OPTIONS: "True, False", RESULT: False}, index=([0]))
                results_df = pd.concat([results_df, comp_df])

        ######### Code for 305, HRD, is not implemented yet #########
        elif marker_name == "HRD":
            get_HRD_results()
        
        else:
            print(f"Composite biomarker: {marker_id} ({marker_name}) is not yet implemented in the code.")
            continue
        
    results_df.sort_values([BIOMARKER_ID], inplace=True)
    results_df.reset_index(drop=True, inplace=True)
        
    return results_df


def raw_panel_csv(snv_df, sv_df, mod_df, immune_df):
    full_panel_results_df = pd.concat([snv_df, sv_df, mod_df, immune_df])
    full_panel_results_df = full_panel_results_df.sort_values(['ID'])
    full_panel_results_df.reset_index(drop=True, inplace=True)

    # not sure why I output this on its own
    # full_panel_results_df.to_csv(out_raw, index=False)
    
    return full_panel_results_df


def identify_groups(panel_results_df, metadata_df):
    # Identify possible molecular groups
    
    MOL_GROUP_HEADER = get_header(metadata_df, 'Molecular Group')
    FUSION_HEADER = get_header(metadata_df, 'Fusion')
    COMPOSITE_COMPONENT_HEADER = get_header(metadata_df, 'Combined Marker')
    
    # convert metadata_df to dictionary for faster lookup
    metadata_dict = metadata_df.set_index('ID').to_dict(orient='index')
    
    # preclin_stage_panel_result_header = [BIOMARKER_ID, BIOMARKER_NAME, SCORING_TYPE, BIOMARKER_TYPE, RESULT_OPTIONS, RESULT]
    
    for idx, row in panel_results_df.iterrows():
        result_id = row[BIOMARKER_ID]
        result_name = row[BIOMARKER_NAME]
        result_scoring_type = row[SCORING_TYPE]
        result_bm_type = row[BIOMARKER_TYPE]
        result_options = row[RESULT_OPTIONS]

        # get molecular groups for biomarker
        metadata_dict_entry = metadata_dict.get(result_id, {})
        mol_group = metadata_dict_entry.get(MOL_GROUP_HEADER, '') # molecular group, if any, separated by ";"
        
        components = metadata_dict_entry.get(COMPOSITE_COMPONENT_HEADER)  # combined marker components, if any, separated by "/"
        fusion_components = metadata_dict_entry.get(FUSION_HEADER)  # fusion components, if any, separated by "/"
        
    return 


def add_demographic_clinicopath_results(path2panel, panel_results_df, metadata_df):
    """load and process demographic and clinicopath results"""
    
    BM_TYPE_FULL = get_BM_TYPE_FULL(path2panel=path2panel)
    SCORING_TYPE_FULL = get_full_SCORING_TYPE(path2panel=path2panel)

    panel_demog = variant_prep(path2panel, variant_type='demographic')
    panel_clinpath = variant_prep(path2panel, variant_type='clinicopathology')

    demclin_df = pd.concat([panel_demog, panel_clinpath])
    if NOTES not in demclin_df.columns:
        demclin_df[NOTES] = ""
    for i, row in demclin_df.iterrows():
        panel_results_df.loc[i] = [row[BIOMARKER_ID], row[BIOMARKER_NAME], row[SCORING_TYPE_FULL], row[BM_TYPE_FULL], row[RESULT_OPTIONS], row[RESULT], row[NOTES]]  # type: ignore
        
    return panel_results_df


def main(args):

    # Get BM_TYPE header
    BM_TYPE_FULL = get_BM_TYPE_FULL(path2panel=args.panel)
    SCORING_TYPE_FULL = get_full_SCORING_TYPE(path2panel=args.panel)
    
    # Get data
    snv_df = pd.read_csv(args.snv, dtype={"ID": str})
    sv_df = pd.read_csv(args.sv, dtype={"ID": str})
    mod_df = pd.read_csv(args.mod, dtype={"ID": str})
    immune_df = pd.read_csv(args.immune, dtype={"ID": str})
    
    # Get panel
    with open(args.panel, "r") as fpc:
        metadata_df = pd.read_csv(fpc, dtype={"ID": str})

    # Get raw outputs for processing through composite biomarker and molecular group identification
    raw_panel_results_df = raw_panel_csv(snv_df, sv_df, mod_df, immune_df)
    
    # Get all composite biomarkers from panel metadata file and work out composite results
    panel_results_df = get_composite_biomarkers(
        metadata_df, BM_TYPE_FULL, raw_panel_results_df)

    # Identify possible molecular groups    
    groups = identify_groups(raw_panel_results_df, metadata_df)  
    # list for now - in dev. Not sure if needed now that I have the composite df

    # Add in the demographic and clinicopathologic results
    panel_results_df = add_demographic_clinicopath_results(args.panel, panel_results_df, metadata_df)

    # Set up the final result dataframe
    # full_panel_results_df = pd.DataFrame(columns=preclin_stage_panel_result_header)

    # full_panel_results_df = pd.concat([snv_df, sv_df, mod_df, immune_df, panel_results_df])
    full_panel_results_df = panel_results_df.sort_values([BIOMARKER_ID])
    full_panel_results_df.reset_index(drop = True, inplace = True)

    
    # save preclin panel
    full_panel_results_df.to_csv(args.out_final, index=False)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel', required=True)
    parser.add_argument('--snv', required=True)
    parser.add_argument('--sv', required=True)
    parser.add_argument('--mod', required=True)
    parser.add_argument('--immune', required=True)
    parser.add_argument('--out_final', required=True)
    args = parser.parse_args()
    
    main(args)
