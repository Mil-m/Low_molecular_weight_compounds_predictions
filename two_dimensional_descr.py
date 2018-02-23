#Group 1----------------------------------------------------------------------------------------------------------------

descriptors_1_group = {'NumRadicalElectrons', 'MinPartialCharge', 'MinEStateIndex', 'MinAbsEStateIndex', 'BalabanJ', 'HallKierAlpha',
               'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA5', 'PEOE_VSA7',
               'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA2', 'SMR_VSA4', 'SMR_VSA8', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11',
               'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA9', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5',
               'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4',
               'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'FractionCSP3', 'NumAromaticHeterocycles','NumHDonors', 'NumHeteroatoms',
               'NumRotatableBonds', 'NumSaturatedRings'}

#Group 2 (Maining)------------------------------------------------------------------------------------------------------

descriptors_2_group = {'BalabanJ', 'BertzCT', 'Chi0n', 'Chi0v', 'Chi1', 'Chi3n', 'HallKierAlpha', 'Kappa2', 'LabuteASA', 'PEOE_VSA1',
                       'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4',
                       'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA3',
                       'SMR_VSA7', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'EState_VSA10', 'VSA_EState8',
                       'VSA_EState9', 'HeavyAtomCount', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'fr_Al_OH', 'fr_Ar_OH',
                       'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_amide', 'fr_aryl_methyl',
                       'fr_azide', 'fr_halogen', 'fr_hdrzone', 'fr_imidazole', 'fr_lactone', 'fr_nitro_arom', 'fr_para_hydroxylation',
                       'fr_phenol', 'fr_phos_acid', 'fr_piperdine', 'fr_sulfone', 'fr_thiazole', 'fr_thiophene', 'fr_unbrch_alkane'}

#Group 3----------------------------------------------------------------------------------------------------------------

descriptors_3_group = {'NumRadicalElectrons', 'MinPartialCharge', 'MinEStateIndex', 'MinAbsEStateIndex', 'BalabanJ', 'HallKierAlpha',
               'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA5', 'PEOE_VSA7',
               'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA2', 'SMR_VSA4', 'SMR_VSA8', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11',
               'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA9', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5',
               'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4',
               'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'FractionCSP3', 'NumAromaticHeterocycles','NumHDonors', 'NumHeteroatoms',
               'NumRotatableBonds', 'NumSaturatedRings', 'BertzCT', 'Chi0n', 'Chi0v', 'Chi1', 'Chi3n', 'Kappa2', 'LabuteASA', 'PEOE_VSA1',
               'PEOE_VSA14', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA3', 'SMR_VSA7', 'SlogP_VSA2', 'SlogP_VSA5',
               'SlogP_VSA6', 'EState_VSA10', 'VSA_EState9', 'HeavyAtomCount', 'fr_Al_OH', 'fr_Ar_OH', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_aldehyde',
               'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_amide', 'fr_aryl_methyl', 'fr_azide', 'fr_halogen', 'fr_hdrzone', 'fr_imidazole',
               'fr_lactone', 'fr_nitro_arom', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phos_acid', 'fr_piperdine', 'fr_sulfone', 'fr_thiazole',
               'fr_thiophene', 'fr_unbrch_alkane'}

#Corr-------------------------------------------------------------------------------------------------------------------

'''corr_descr = {'MolecularFormula', 'Chi0n', 'Chi1n', 'Chi2n', 'Chi3n', 'Kappa3', 'PEOE_VSA7', 'SlogP_VSA12', 'SlogP_VSA4', 'SlogP_VSA7', 'VSA_EState10', 'MolLogP', 'fr_aryl_methyl', 'fr_halogen', 'fr_para_hydroxylation',\
        'fr_phenol_noOrthoHbond', 'SMR_VSA5', 'SlogP_VSA5', 'EState_VSA3', 'EState_VSA4', 'PEOE_VSA8', 'SMR_VSA6', 'SlogP_VSA1', 'SlogP_VSA10', 'fr_ArN', 'fr_NH2', 'fr_aniline', 'EState_VSA6', 'EState_VSA9',\
        'VSA_EState9', 'SMR_VSA7', 'SMR_VSA10', 'SlogP_VSA6', 'SlogP_VSA3', 'EState_VSA5', 'EState_VSA8', 'NumRotatableBonds', 'Chi0v', 'Kappa2', 'LabuteASA', 'PEOE_VSA14', 'PEOE_VSA6', 'EState_VSA1', 'SMR_VSA1',\
        'TPSA', 'EState_VSA2', 'NOCount', 'EState_VSA7', 'fr_Ar_OH', 'fr_phenol', 'SlogP_VSA2', 'NumHeteroatoms', 'PEOE_VSA1', 'HallKierAlpha', 'PEOE_VSA10', 'NumHAcceptors', 'PEOE_VSA3', 'PEOE_VSA9', 'fr_Ar_COO',\
        'fr_COO', 'fr_COO2', 'fr_ester', 'fr_ether', 'NHOHCount', 'NumHDonors', 'fr_C_O_noCOO'}'''
