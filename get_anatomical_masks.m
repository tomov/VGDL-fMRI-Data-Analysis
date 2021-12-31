function [mask_filenames, region, roi_masks] = get_anatomical_masks(atlas_name)

roi_labels_HarvardOxford = { ...
    {'Frontal Pole'}, ...
    {'Inferior Frontal Gyrus, pars triangularis'}, ...
    {'Inferior Frontal Gyrus, pars opercularis'}, ...
    {'Middle Frontal Gyrus'}, ...
    {'Superior Frontal Gyrus'}, ...
    {'Precentral Gyrus'}, ...
    {'Postcentral Gyrus'}, ...
    {'Supramarginal Gyrus, anterior division'}, ...
    {'Precuneous Cortex'}, ...
    {'Superior Parietal Lobule'}, ...
    {'Parietal Operculum Cortex'}, ...
    {'Lateral Occipital Cortex, superior division'}, ...
    {'Lateral Occipital Cortex, inferior division'}, ...
    {'Occipital Pole'}, ...
    {'Occipital Fusiform Gyrus'}, ...
    {'Intracalcarine Cortex'}, ...
    {'Lingual Gyrus'}, ...
    {'Temporal Occipital Fusiform Cortex'}, ...
    {'Planum Temporale'}, ...
    {'Superior Temporal Gyrus, anterior division'}, ...
    {'Superior Temporal Gyrus, posterior division'}, ...
};

roi_names_HarvardOxford = { ...
    'RLPFC', ...
    'IFG_Tri', ...
    'IFG_Oper', ...
    'Frontal_Mid', ...
    'Frontal_Sup', ...
    'Precentral', ...
    'Postcentral', ...
    'SMG_ant', ...
    'Precuneus', ...
    'Parietal_Sup', ...
    'Parietal_Oper', ...
    'LOC_Sup', ...
    'LOC_Inf', ...
    'Occipital_Pole', ...
    'Occipital_Fusiform', ...
    'Intracalcarine', ...
    'Lingual', ...
    'Temporal_Fusiform', ...
    'Planum', ...
    'STG_ant', ...
    'STG_post', ...
};

roi_labels_AAL2 = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}, ...
    {'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Frontal_Mid_2_L', 'Frontal_Mid_2_R'}, ...
    {'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
    {'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R'}, ...
    {'Precentral_L', 'Precentral_R'}, ...
    {'Postcentral_L', 'Postcentral_R'}, ...
    {'Precuneus_L', 'Precuneus_R'}, ...
    {'Cingulate_Mid_L', 'Cingulate_Mid_R'}, ...
    {'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Rolandic_Oper_L', 'Rolandic_Oper_R'}, ...
    {'Temporal_Sup_L', 'Temporal_Sup_R'}, ...
    {'Temporal_Mid_L', 'Temporal_Mid_R'}, ...
    {'Fusiform_L', 'Fusiform_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R'}, ...
    {'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R'}, ...
    {'Lingual_L', 'Lingual_R'}, ...
    {'Calcarine_L', 'Calcarine_R'}, ...
    {'Cuneus_L', 'Cuneus_R'} ...
};

roi_names_AAL2 = { ...
    'IFG_Tri', ...
    'IFG_Oper', ...
    'Frontal_Mid', ...
    'Frontal_Sup', ...
    'Frontal_Sup_Medial', ...
    'Precentral', ...
    'Postcentral', ...
    'Precuneus', ...
    'Cingulate_Mid', ...
    'SMA', ...
    'Rolandic_Oper', ...
    'Temporal_Sup', ...
    'Temporal_Mid', ...
    'Fusiform', ...
    'Occipital_Sup', ...
    'Occipital_Mid', ...
    'Occipital_Inf', ...
    'Lingual', ...
    'Calcarine', ...
    'Cuneus', ...
};


roi_labels_AAL3 = [ ...
    { ...
    {'OFCmed_L', 'OFCmed_R'}, ...
    {'OFCant_L', 'OFCant_R'}, ...
    {'OFCpost_L', 'OFCpost_R'}, ...
    {'OFClat_L', 'OFClat_R'}, ...
    {'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R'}, ...
    {'Rectus_L', 'Rectus_R'}, ...
    {'Hippocampus_L', 'Hippocampus_R'}, ...
    {'N_Acc_L', 'N_Acc_R'}, ...
    {'VTA_L', 'VTA_R'}, ...
    {'SN_pc_L', 'SN_pc_R'}, ...
    },  ...
    roi_labels_AAL2, ...
];

roi_names_AAL3 = [ ...
    { ...
    'OFCmed', ...
    'OFCant', ...
    'OFCpost', ...
    'OFClat', ...
    'Frontal_Med_Orb', ...
    'Rectus', ...
    'Hippocampus', ...
    'N_Acc', ...
    'VTA', ...
    'SN_pc', ...
    }, ...
    roi_names_AAL2, ...
];

roi_labels_AAL2_grouped = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Mid_2_L', 'Frontal_Mid_2_R', 'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
    {'Precentral_L', 'Precentral_R', 'Postcentral_L', 'Postcentral_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Temporal_Sup_L', 'Temporal_Sup_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R'}, ...
    {'Lingual_L', 'Lingual_R', 'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R'} ...
    {'Cingulate_Mid_L', 'Cingulate_Mid_R', 'Precuneus_L', 'Precuneus_R'}, ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Frontal_Mid_2_L', 'Frontal_Mid_2_R', 'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
};

roi_names_AAL2_grouped = { ...
    'Frontal', ...
    'Sensory_Motor', ...
    'Dorsal_Parietal', ...
    'Ventral_Temporal', ...
    'Early visual', ...
    'Posteromedial', ...
    'Inferior frontal', ...
    'Superior frontal', ...
};

roi_labels_AAL2_grouped2 = { ...
    {'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R'}, ...
    {'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', 'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Lingual_L', 'Lingual_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Lingual_L', 'Lingual_R', 'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Precuneus_L', 'Precuneus_R'}, ...
    {'Precentral_L', 'Precentral_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Precentral_L', 'Precentral_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Mid_2_L', 'Frontal_Mid_2_R', 'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Mid_2_L', 'Frontal_Mid_2_R', 'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
};

roi_names_AAL2_grouped2 = { ...
    'Early visual', ...
    'Early visual + MOC', ...
    'Ventral_Temporal', ...
    'Ventral_Temporal + MOC', ...
    'Dorsal_Parietal', ...
    'Motor', ...
    'Motor + Frontal', ...
    'Frontal',
};


% for GP
roi_labels_AAL2_GP_EMPA_grouped = { ...
    {'Precentral_L', 'Precentral_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Mid_2_L', 'Frontal_Mid_2_R', 'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
    {'Angular_L', 'Angular_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Precuneus_L', 'Precuneus_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Lingual_L', 'Lingual_R', 'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R'}, ...
};

roi_names_AAL2_GP_EMPA_grouped = { ...
    'Motor + Frontal', ...
    'Dorsal_Parietal', ...
    'Ventral_Temporal', ...
    'Early visual', ...
};



% Goes together with AAL2_GP_EMPA_grouped
roi_labels_AAL2_GP_EMPA = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}, ...
    {'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Frontal_Mid_2_L', 'Frontal_Mid_2_R'}, ...
    {'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
    {'Precentral_L', 'Precentral_R'}, ...
    {'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Angular_L', 'Angular_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R'}, ... % ------------
    {'SupraMarginal_L', 'SupraMarginal_R'}, ...
    {'Precuneus_L', 'Precuneus_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R'}, ... % ----------
    {'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Fusiform_L', 'Fusiform_R'}, ...
    {'Lingual_L', 'Lingual_R'}, ...
    {'Calcarine_L', 'Calcarine_R'}, ... % --------------
    {'Cuneus_L', 'Cuneus_R'}, ...
};

roi_names_AAL2_GP_EMPA = { ...
    'IFGtriang', ...
    'IFGoperc', ...
    'MFG', ...
    'SFG', ...
    'PreCG', ...
    'SMA', ...
    'AG', ...
    'SOG', ...
    'SMG', ...
    'PCUN', ...
    'IOG', ...
    'MOG', ...
    'FFG', ...
    'LING', ...
    'CAL', ...
    'CUN', ...
};


% for GLM
roi_labels_AAL2_GLM_102_grouped = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'OFCant_L', 'OFCant_R','Precentral_L', 'Precentral_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Angular_L', 'Angular_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Precuneus_L', 'Precuneus_R', 'Cingulate_Mid_L', 'Cingulate_Mid_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Lingual_L', 'Lingual_R', 'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R'}, ...
};

roi_names_AAL2_GLM_102_grouped = { ...
    'Motor + Frontal', ...
    'Dorsal_Parietal', ...
    'Ventral_Temporal', ...
    'Early visual', ...
};

% for GLM
% Goes with roi_names_AAL2_GLM_102_grouped
roi_labels_AAL2_GLM_102 = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}, ...
    {'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R'}, ...
    {'OFCant_L', 'OFCant_R'}, ...
    {'Precentral_L', 'Precentral_R'}, ...
    {'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Angular_L', 'Angular_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R'}, ... % ------------
    {'Precuneus_L', 'Precuneus_R'}, ...
    {'Cingulate_Mid_L', 'Cingulate_Mid_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R'}, ... % ----------
    {'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Fusiform_L', 'Fusiform_R'}, ...
    {'Lingual_L', 'Lingual_R'}, ...
    {'Calcarine_L', 'Calcarine_R'}, ... % --------------
    {'Cuneus_L', 'Cuneus_R'}, ...
};

roi_names_AAL2_GLM_102 = { ...
    'IFGtriang', ...
    'IFGoperc', ...
    'SFGmedial', ...
    'OFCant', ...
    'PreCG', ...
    'SMA', ...
    'AG', ...
    'SOG', ...
    'PCUN', ...
    'MCC', ...
    'IOG', ...
    'MOG', ...
    'FFG', ...
    'LING', ...
    'CAL', ...
    'CUN', ...
};


% for GP intersect GLM
roi_labels_AAL2_GP_EMPA_GLM_102_grouped = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Angular_L', 'Angular_R', 'Occipital_Sup_L', 'Occipital_Sup_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
};

roi_names_AAL2_GP_EMPA_GLM_102_grouped = { ...
    'Motor + Frontal', ...
    'Dorsal_Parietal', ...
    'Ventral_Temporal', ...
};



% Goes together with AAL2_GP_EMPA_GLM_102_grouped
roi_labels_AAL2_GP_EMPA_GLM_102 = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}, ...
    {'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Angular_L', 'Angular_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R'}, ... % ------------
    {'Occipital_Inf_L', 'Occipital_Inf_R'}, ... % ----------
    {'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Fusiform_L', 'Fusiform_R'}, ...
};

roi_names_AAL2_GP_EMPA_GLM_102 = { ...
    'IFGtriang', ...
    'IFGoperc', ...
    'AG', ...
    'SOG', ...
    'IOG', ...
    'MOG', ...
    'FFG', ...
};

roi_labels_Brodmann = { ...
    {'BA 17'}, ...
    {'BA 18'}, ...
    {'BA 19'}, ...
    {'BA 17', 'BA 18', 'BA 19'}, ...
};

roi_names_Brodmann = { ...
    'V1', ...
    'V2', ...
    'V3+4+5', ...
    'Early visual'
};


% for functional connectivity analysis for theory change flag
roi_labels_AAL2_TETRAD_GLM_109 = {
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}, ...
    {'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Frontal_Mid_2_L', 'Frontal_Mid_2_R'}, ...
    {'Frontal_Sup_2_L', 'Frontal_Sup_2_R'}, ...
    {'Supp_Motor_Area_L', 'Supp_Motor_Area_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R'}, ... % ------------
    {'Occipital_Inf_L', 'Occipital_Inf_R'}, ... % ----------
    {'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Fusiform_L', 'Fusiform_R'}, ...
    {'Lingual_L', 'Lingual_R'}, ...
    {'Calcarine_L', 'Calcarine_R'}, ... % --------------
    {'Cuneus_L', 'Cuneus_R'}, ...
};

roi_names_AAL2_TETRAD_GLM_109 = {
    'IFGtriang', ...
    'IFGoperc', ...
    'MFG', ...
    'SFG', ...
    'SMA', ...
    'SOG', ...
    'IOG', ...
    'MOG', ...
    'FFG', ...
    'LING', ...
    'CAL', ...
    'CUN', ...
};

% for GLM 157
%

roi_labels_AAL2_GLM_157 = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}, ...
    {'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R'}, ...
    {'Frontal_Inf_Orb_2_L', 'Frontal_Inf_Orb_2_R'}, ...
    {'Precentral_L', 'Precentral_R'}, ...
    {'Parietal_Inf_L', 'Parietal_Inf_R'}, ...
    {'Angular_L', 'Angular_R'}, ...
    {'Occipital_Sup_L', 'Occipital_Sup_R'}, ... % ------------
    {'Precuneus_L', 'Precuneus_R'}, ...
    {'Cingulate_Mid_L', 'Cingulate_Mid_R'}, ...
    {'Occipital_Inf_L', 'Occipital_Inf_R'}, ... % ----------
    {'Occipital_Mid_L', 'Occipital_Mid_R'}, ...
    {'Fusiform_L', 'Fusiform_R'}, ...
    {'Lingual_L', 'Lingual_R'}, ...
    {'Calcarine_L', 'Calcarine_R'}, ... % --------------
    {'Cuneus_L', 'Cuneus_R'}, ...
};

roi_names_AAL2_GLM_157 = { ...
    'IFGtriang', ...
    'IFGoperc', ...
    'IFGorb', ...
    'PreCG', ...
    'IPG', ...
    'AG', ...
    'SOG', ...
    'PCUN', ...
    'MCC', ...
    'IOG', ...
    'MOG', ...
    'FFG', ...
    'LING', ...
    'CAL', ...
    'CUN', ...
};

roi_labels_AAL2_GLM_157_all = { ...
    {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', ...
    'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', ...
    'Frontal_Inf_Orb_2_L', 'Frontal_Inf_Orb_2_R', ...
    'Precentral_L', 'Precentral_R', ...
    'Parietal_Inf_L', 'Parietal_Inf_R', ...
    'Angular_L', 'Angular_R', ...
    'Occipital_Sup_L', 'Occipital_Sup_R', ...
    'Precuneus_L', 'Precuneus_R', ...
    'Cingulate_Mid_L', 'Cingulate_Mid_R', ...
    'Occipital_Inf_L', 'Occipital_Inf_R', ...
    'Occipital_Mid_L', 'Occipital_Mid_R', ...
    'Fusiform_L', 'Fusiform_R', ...
    'Lingual_L', 'Lingual_R', ...
    'Calcarine_L', 'Calcarine_R', ...
    'Cuneus_L', 'Cuneus_R'}, ...
};

roi_names_AAL2_GLM_157_all = { ...
    'all', ...
};


switch atlas_name
    case 'AAL2'
        roi_labels = roi_labels_AAL2;
        roi_names = roi_names_AAL2;
    case 'AAL2_grouped'
        roi_labels = roi_labels_AAL2_grouped;
        roi_names = roi_names_AAL2_grouped;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_grouped2'
        roi_labels = roi_labels_AAL2_grouped2;
        roi_names = roi_names_AAL2_grouped2;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GP_EMPA_grouped'
        roi_labels = roi_labels_AAL2_GP_EMPA_grouped;
        roi_names = roi_names_AAL2_GP_EMPA_grouped;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GP_EMPA'
        roi_labels = roi_labels_AAL2_GP_EMPA;
        roi_names = roi_names_AAL2_GP_EMPA;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GLM_102_grouped'
        roi_labels = roi_labels_AAL2_GLM_102_grouped;
        roi_names = roi_names_AAL2_GLM_102_grouped;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GLM_102'
        roi_labels = roi_labels_AAL2_GLM_102;
        roi_names = roi_names_AAL2_GLM_102;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GLM_157'
        roi_labels = roi_labels_AAL2_GLM_157;
        roi_names = roi_names_AAL2_GLM_157;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GLM_157_all'
        roi_labels = roi_labels_AAL2_GLM_157_all;
        roi_names = roi_names_AAL2_GLM_157_all;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GP_EMPA_GLM_102_grouped'
        roi_labels = roi_labels_AAL2_GP_EMPA_GLM_102_grouped;
        roi_names = roi_names_AAL2_GP_EMPA_GLM_102_grouped;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_GP_EMPA_GLM_102'
        roi_labels = roi_labels_AAL2_GP_EMPA_GLM_102;
        roi_names = roi_names_AAL2_GP_EMPA_GLM_102;
        atlas_name = 'AAL2'; % proper atlas name
    case 'AAL2_TETRAD_GLM_109'
        roi_labels = roi_labels_AAL2_TETRAD_GLM_109; 
        roi_names = roi_names_AAL2_TETRAD_GLM_109;
        atlas_name = 'AAL2';% proper atlas name
    case 'AAL3v1'
        roi_labels = roi_labels_AAL3;
        roi_names = roi_names_AAL3;
    case 'HarvardOxford'
        roi_labels = roi_labels_HarvardOxford;
        roi_names = roi_names_HarvardOxford;
    case 'Brodmann'
        roi_labels = roi_labels_Brodmann;
        roi_names = roi_names_Brodmann;
    otherwise
        assert(false, 'no such Atlas');
end

nROIs = length(roi_labels)
assert(length(roi_names) == nROIs);

for m = 1:nROIs
    mask_filenames{m} = fullfile('masks', [roi_names{m}, '.nii']);
    roi_masks{m} = ccnl_create_mask(roi_labels{m}, mask_filenames{m},  atlas_name, true, 'masks/mask.nii');
    %roi_masks_flattened{m} = roi_masks{m}(whole_brain_mask);
end

region = roi_names';
