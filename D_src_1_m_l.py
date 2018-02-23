from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import GraphDescriptors
#from rdkit.Chem import PandasTools
import collections
import xlrd
import csv
import time
import os
import pybel
import openbabel
import xml.etree.ElementTree as ET
from itertools import islice

import two_dimensional_descr
import three_dimensional_descr

#-----------------------------------------------------------------------------------------------------------------------settings
training_group = 1
rewrite_f_training = 0
rewrite_f_prediction = 1
get_sdf_files = 0
get_csv_groups = 0
startTime = time.time()
#input_csv = '/home/mi/PycharmProjects/D_scr_1_m_l/csv_input/AlfaAesarCurated.csv'
#input_csv = '/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_c_O.csv'
#input_csv = '/home/mi/PycharmProjects/D_scr_1_m_l/csv_input/MP_PSA.csv'
input_csv = '/home/mi/PycharmProjects/D_scr_1_m_l/csv_input/predicted_descr.csv'

#-----------------------------------------------------------------------------------------------------------------------descriptors

descr_group = 'Predicted descriptors'

descriptors = []
if (descr_group == 1):
    descriptors = two_dimensional_descr.descriptors_1_group
elif (descr_group == 2):
    descriptors = two_dimensional_descr.descriptors_2_group
elif (descr_group == 3):
    descriptors = two_dimensional_descr.descriptors_3_group
elif (descr_group == '3d'):
    print "Three-dimensional descriptors... Evaluating"
elif (descr_group == 'TPSA'):
    print "TPSA as a descriptor... Evaluating"
elif (descr_group == 'Predicted descriptors'):
    print "Predicted descriptors... Evaluating"
else:
    print "Err in group of descriptors"
    exit()

#-----------------------------------------------------------------------------------------------------------------------input_csv

first_molnames, first_smolecules, first_melting_point, first_logP, \
first_TPSA, first_pred_logP, first_pred_logS, first_pred_PSA, first_pred_refractivity, first_pred_polarizability  = [], [], [], [], [], [], [], [], [], []
second_molnames, second_smolecules, second_melting_point, second_logP, \
second_TPSA, second_pred_logP, second_pred_logS, second_pred_PSA, second_pred_refractivity, second_pred_polarizability = [], [], [], [], [], [], [], [], [], []

csvf = open(input_csv, 'rb')
rows = csv.reader(csvf, delimiter=',', quotechar='"')
rownum = 0
for row in rows:
    rownum += 1
len_csv_file = rownum - 1
csvf.close()

delimiter = len_csv_file - len_csv_file / 3
csvf = open(input_csv, 'rb')
rows = csv.reader(csvf, delimiter=',', quotechar='"')
rownum = 0
for row in rows:
    if rownum != 0:
        if rownum <= delimiter:
            first_molnames.append(row[0].replace(' ', '_'))
            mol = Chem.MolFromSmiles(str(row[1]).replace('"', ''), sanitize=False)
            first_smolecules.append(Chem.MolToSmiles(mol))
            first_melting_point.append(row[2])
            first_logP.append(0)
            if (descr_group == 'TPSA'):
                first_TPSA.append(row[3])
            if (descr_group == 'Predicted descriptors'):
                first_pred_logP.append(row[3])
                first_pred_logS.append(row[4])
                first_pred_PSA.append(row[5])
                first_pred_refractivity.append(row[6])
                first_pred_polarizability.append(row[7])
        else:
            second_molnames.append(row[0].replace(' ', '_'))
            mol = Chem.MolFromSmiles(str(row[1]).replace('"', ''), sanitize=False)
            second_smolecules.append(Chem.MolToSmiles(mol))
            second_melting_point.append(row[2])
            second_logP.append(0)
            if (descr_group == 'TPSA'):
                second_TPSA.append(row[3])
            if (descr_group == 'Predicted descriptors'):
                second_pred_logP.append(row[3])
                second_pred_logS.append(row[4])
                second_pred_PSA.append(row[5])
                second_pred_refractivity.append(row[6])
                second_pred_polarizability.append(row[7])
    rownum += 1

csvf.close()

#-----------------------------------------------------------------------------------------------------------------------write etalon

f_etalon_MP_logP_gr_1 = open('/home/mi/PycharmProjects/D_scr_1_m_l/etalon_MP_logP_gr_1.csv', 'w')
f_etalon_MP_logP_gr_2 = open('/home/mi/PycharmProjects/D_scr_1_m_l/etalon_MP_logP_gr_2.csv', 'w')

f_etalon_MP_logP_gr_1.write('{0};{1}'.format("Num;Name;Smiles;Melting_point;LogP", "\n"))
for i in range(len(first_smolecules)):
    str_first_melting_point = str(first_melting_point[i])
    str_first_logP = str(first_logP[i])
    f_etalon_MP_logP_gr_1.write('{0};{1};{2};{3};{4}{5}'.format(str(i+1), str(first_molnames[i]), str(first_smolecules[i]), str_first_melting_point.replace(".", ","), str_first_logP.replace(".", ","), "\n"))

f_etalon_MP_logP_gr_2.write('{0};{1}'.format("Num;Name;Smiles;Melting_point;LogP", "\n"))
for i in range(len(second_smolecules)):
    str_second_melting_point = str(second_melting_point[i])
    str_second_logP = str(second_logP[i])
    f_etalon_MP_logP_gr_2.write('{0};{1};{2};{3};{4}{5}'.format(str(i+1), str(second_molnames[i]), str(second_smolecules[i]), str_second_melting_point.replace(".", ","), str_second_logP.replace(".", ","), "\n"))

f_etalon_MP_logP_gr_1.close()
f_etalon_MP_logP_gr_2.close()

#-----------------------------------------------------------------------------------------------------------------------training_table

if rewrite_f_training :
    f_training = open('/home/mi/PycharmProjects/D_scr_1_m_l/training_table.csv', 'w')
else:
    f_training = open('/home/mi/PycharmProjects/D_scr_1_m_l/training_table.csv', 'r')

if training_group == 1:
    training_molnames = first_molnames
    control_molnames = second_molnames
    training_smolecules = first_smolecules
    control_smolecules = second_smolecules
    training_melting_point = first_melting_point
    control_melting_point = second_melting_point
    training_logP = first_logP
    control_logP = second_logP
elif training_group == 2:
    training_molnames = second_molnames
    control_molnames = first_molnames
    training_smolecules = second_smolecules
    control_smolecules = first_smolecules
    training_melting_point = second_melting_point
    control_melting_point = first_melting_point
    training_logP = second_logP
    control_logP = first_logP
else:
    print "Err in number group"
    exit()

if rewrite_f_training:

    if (descr_group == 1 or descr_group == 2 or descr_group == 3):
        real_count_of_descriptors = 0

        for descName, descFn in Descriptors._descList:
            if (descName in descriptors):
                real_count_of_descriptors += 1

        descDict = collections.defaultdict(list)
        descDictStr = ""
        descDictStr_for_table = ""

        j = 0
        for descName, descFn in Descriptors._descList:
            descValue = descFn(Chem.MolFromSmiles(training_smolecules[0]))
            descDict[descName] = descValue
            if (descName in descriptors):
                j += 1
                if (j < real_count_of_descriptors):
                    descDictStr += str(descName) + ";"
                    descDictStr_for_table += str(descName) + "+"
                else:
                    descDictStr += str(descName)
                    descDictStr_for_table += str(descName)

        descDictStr = "Num;Smiles;Melting_point;LogP;" + descDictStr
        print descDictStr_for_table

        f_training.write('{0}{1}'.format(descDictStr, "\n"))
        for i in range(len(training_smolecules)):
            f_training.write('{0};{1};{2};{3};'.format(str(i+1), "'" + str(training_smolecules[i]) + "'", str(training_melting_point[i]), str(training_logP[i])))
            j = 0
            for descName, descFn in Descriptors._descList:
                descValue = descFn(Chem.MolFromSmiles(training_smolecules[i]))
                if (descName in descriptors):
                    j += 1
                    if (j < real_count_of_descriptors):
                        f_training.write('{0};'.format(str(descValue)))
                    else:
                        f_training.write('{0}{1}'.format(str(descValue), "\n"))

    elif (descr_group == '3d'):
        descDictStr_for_table = "max_d+min_d+avg_d+med_d+rms_d+" + \
                      "max_v+min_v+avg_v+med_v+rms_v+" + \
                      "max_p+min_p+" + \
                      "mol_wt+box_volume+wt_vol_ratio+mol_volume+volume_ratio+" + \
                      "min_m+max_m+dc_m+" + \
                      "dip_l+dip_pos_l+dip_pos_lq+" + \
                      "dip_neg_l+dip_neg_lq+dip_angle"
        print descDictStr_for_table

        descDictStr = "Num;Smiles;Melting_point;LogP;" + \
                      "max_d;min_d;avg_d;med_d;rms_d;" + \
                      "max_v;min_v;avg_v;med_v;rms_v;" + \
                      "max_p;min_p;" + \
                      "mol_wt;box_volume;wt_vol_ratio;mol_volume;volume_ratio;" + \
                      "min_m;max_m;dc_m;" + \
                      "dip_l;dip_pos_l;dip_pos_lq;" + \
                      "dip_neg_l;dip_neg_lq;dip_angle"

        f_training.write('{0}{1}'.format(descDictStr, "\n"))
        for i in range(len(training_smolecules)):
            f_training.write('{0};{1};{2};{3};'.format(str(i+1), "'" + str(training_smolecules[i]) + "'", str(training_melting_point[i]), str(training_logP[i])))
            arr_3d_descriptors = three_dimensional_descr.count_three_dimensional(training_smolecules[i])
            j = 0
            for descName, descValue in arr_3d_descriptors.iteritems():
                j += 1
                if (j < len(arr_3d_descriptors)):
                    f_training.write('{0};'.format(str(descValue)))
                else:
                    f_training.write('{0}{1}'.format(str(descValue), "\n"))
    elif (descr_group == 'TPSA'):
        descDictStr_for_table = "TPSA"
        print descDictStr_for_table
        descDictStr = "Num;Smiles;Melting_point;LogP;TPSA"
        f_training.write('{0}{1}'.format(descDictStr, "\n"))
        for i in range(len(training_smolecules)):
            f_training.write('{0};{1};{2};{3};{4}{5}'.format(str(i+1), "'" + str(training_smolecules[i]) + "'", str(training_melting_point[i]), str(training_logP[i]), str(first_TPSA[i]), "\n"))
    elif (descr_group == 'Predicted descriptors'):
        descDictStr_for_table = "pred_logP+pred_logS+pred_PSA+pred_refractivity+pred_polarizability"
        print descDictStr_for_table
        descDictStr = "Num;Smiles;Melting_point;LogP;pred_logP;pred_logS;pred_PSA;pred_refractivity;pred_polarizability"
        f_training.write('{0}{1}'.format(descDictStr, "\n"))
        for i in range(len(training_smolecules)):
            f_training.write('{0};{1};{2};{3};{4};{5};{6};{7};{8}{9}'.format(str(i+1), "'" + str(training_smolecules[i]) + "'", str(training_melting_point[i]), str(training_logP[i]),
                                                             str(first_pred_logP[i]), str(first_pred_logS[i]), str(first_pred_PSA[i]),
                                                                     str(first_pred_refractivity[i]), str(first_pred_polarizability[i]), "\n"))
    else:
        print "Err in group of descriptors"
        exit()

f_training.close()

#-----------------------------------------------------------------------------------------------------------------------prediction
if not os.path.exists('/home/mi/PycharmProjects/D_scr_1_m_l/coefficients.txt'):
   print "Err! File 'coefficients.txt' isn't exist"
   exit()

f_control = open('/home/mi/PycharmProjects/D_scr_1_m_l/control.csv', 'w')

if rewrite_f_prediction:
    f_coeff = open('/home/mi/PycharmProjects/D_scr_1_m_l/coefficients.txt', 'r')
    descr_coeff_dict = dict(short='dict', long='dictionary')
    for line in f_coeff:
        coeff_str = line.split(' ')
        descr = coeff_str[0]
        coeff = 0
        i = 0
        for el in coeff_str:
            if (i != 0 and el != ''):
                coeff = el
                break
            i += 1
        descr_coeff_dict[str(descr)] = coeff
    f_coeff.close()

    if (descr_group == 1 or descr_group == 2 or descr_group == 3):
        print "Out . . ."
        f_control.write('{0}{1}'.format("Num;Name;Smiles;Melting_point", "\n"))
        for i in range(len(control_smolecules)):
            regr = 0
            for descName, descFn in Descriptors._descList:
                if descName in descr_coeff_dict:
                    regr += float(descr_coeff_dict[descName]) * float(descFn(Chem.MolFromSmiles(control_smolecules[i])))
            f_control.write('{0};{1};{2};{3}{4}'.format(str(i+1), str(control_molnames[i]), str(control_smolecules[i]), str(regr).replace(".", ","), "\n"))
    elif (descr_group == '3d'):
        print "Out . . ."
        f_control.write('{0}{1}'.format("Num;Name;Smiles;Melting_point", "\n"))
        for i in range(len(control_smolecules)):
            regr = 0
            arr_3d_descriptors = three_dimensional_descr.count_three_dimensional(control_smolecules[i])
            print arr_3d_descriptors
            for descName, descValue in arr_3d_descriptors.iteritems():
                if descName in descr_coeff_dict:
                    regr += float(descr_coeff_dict[descName]) * float(descValue)
            print regr
            f_control.write('{0};{1};{2};{3}{4}'.format(str(i+1), str(control_molnames[i]), str(control_smolecules[i]), str(regr).replace(".", ","), "\n"))
    elif (descr_group == 'TPSA'):
        for i in range(len(control_smolecules)):
            regr = float(descr_coeff_dict['TPSA']) * float(second_TPSA[i])
            f_control.write('{0};{1};{2};{3}{4}'.format(str(i+1), str(control_molnames[i]), str(control_smolecules[i]), str(regr).replace(".", ","), "\n"))
    elif (descr_group == 'Predicted descriptors'):
        for i in range(len(control_smolecules)):
            regr = float(descr_coeff_dict['pred_logP']) * float(second_pred_logP[i]) + float(descr_coeff_dict['pred_logS']) * float(second_pred_logS[i]) + \
                    float(descr_coeff_dict['pred_PSA']) * float(second_pred_PSA[i]) + float(descr_coeff_dict['pred_refractivity']) * float(second_pred_refractivity[i]) + \
                    float(descr_coeff_dict['pred_polarizability']) * float(second_pred_polarizability[i])
            f_control.write('{0};{1};{2};{3}{4}'.format(str(i+1), str(control_molnames[i]), str(control_smolecules[i]), str(regr).replace(".", ","), "\n"))
    else:
        print "Err in group of descriptors"
        exit()

f_control.close()

#-----------------------------------------------------------------------------------------------------------------------write sdf-files

if get_sdf_files:
    f_sdf_gr_1 = pybel.Outputfile('sdf', '/home/mi/PycharmProjects/D_scr_1_m_l/f_sdf_gr_1.sdf', 'w')
    for el in first_smolecules:
        mol = pybel.readstring('smi', el)
        f_sdf_gr_1.write(mol)
    f_sdf_gr_1.close()
    f_sdf_gr_2 = pybel.Outputfile('sdf', '/home/mi/PycharmProjects/D_scr_1_m_l/f_sdf_gr_2.sdf', 'w')
    for el in second_smolecules:
        mol = pybel.readstring('smi', el)
        f_sdf_gr_2.write(mol)
    f_sdf_gr_2.close()

#-----------------------------------------------------------------------------------------------------------------------write csv_groups

if get_csv_groups:
    if not os.path.exists('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups'):
        print "Creating csv groups directory..."
        os.mkdir('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups')

    result_names_arr = first_molnames + second_molnames
    result_smiles_arr = first_smolecules + second_smolecules
    result_m_p = first_melting_point + second_melting_point

    input_gr_C_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C.csv', 'w')
    input_gr_C_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))
    input_gr_C_c_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_c.csv', 'w')
    input_gr_C_c_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))
    input_gr_C_O_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_O.csv', 'w')
    input_gr_C_O_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))
    input_gr_C_c_O_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_c_O.csv', 'w')
    input_gr_C_c_O_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))
    input_gr_C_c_S_O_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_c_S_O.csv', 'w')
    input_gr_C_c_S_O_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))
    input_gr_C_N_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_N.csv', 'w')
    input_gr_C_N_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))
    input_gr_C_N_c_n_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_N_c_n.csv', 'w')
    input_gr_C_N_c_n_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))
    input_gr_C_N_c_n_O_csv = open('/home/mi/PycharmProjects/D_scr_1_m_l/csv_groups/gr_C_N_c_n_O.csv', 'w')
    input_gr_C_N_c_n_O_csv.write('{0}{1}'.format("name,SMILES,mpC", "\n"))

    i = 0
    for el in result_smiles_arr:
        if (el.count('C') == len(el)):
            input_gr_C_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        if (el.count('C') + el.count('c') == len(el)):
            input_gr_C_c_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        if (el.count('C') + el.count('O') == len(el)):
            input_gr_C_O_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        if (el.count('C') + el.count('c') + el.count('O') == len(el)):
            input_gr_C_c_O_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        if (el.count('C') + el.count('c') + el.count('S') + el.count('O') == len(el)):
            input_gr_C_c_S_O_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        if (el.count('C') + el.count('N') == len(el)):
            input_gr_C_N_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        if (el.count('C') + el.count('N') + el.count('c') + el.count('n') == len(el)):
            input_gr_C_N_c_n_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        if (el.count('C') + el.count('N') + el.count('c') + el.count('n') + el.count('O') == len(el)):
            input_gr_C_N_c_n_O_csv.write('{0}{1}{2},{3},{4}{5}'.format('"', result_names_arr[i], '"', el, result_m_p[i], "\n"))
        i += 1

    input_gr_C_csv.close()
    input_gr_C_c_csv.close()
    input_gr_C_O_csv.close()
    input_gr_C_c_O_csv.close()
    input_gr_C_c_S_O_csv.close()
    input_gr_C_N_csv.close()
    input_gr_C_N_c_n_csv.close()
    input_gr_C_N_c_n_O_csv.close()

#-----------------------------------------------------------------------------------------------------------------------

print "Time:", time.time() - startTime