import argparse
from annotate_clusters_consistency import write_csv, read_csv_to_list, create_dict_from_list
import pandas as pd
import os
import os.path



def read_branch_model(branch_file: str) -> list:
    with open(branch_file, 'r') as branch_log:
        found_background = False
        found_1 = False
        branch_log = branch_log.readlines()
        M0_vs_bfree = '-'
        bneut_vs_bfree = '-'
        for line in branch_log:
            if ('M0~' in line) and ('b_free' in line):
                M0_vs_bfree = line.split('|')[2][1:9]
            if ('b_neut' in line) and ('b_free' in line):
                bneut_vs_bfree = line.split('|')[2][1:9]
            if ('background  =>' in line) and (not found_background):
                omega_b = line.split('=>')[1][2:8]
                found_background = True
            if ('#1  =>' in line) and (not found_1):
                omega_1 = line.split('=>')[1][2:8]
                found_1 = True
        if not found_background:
            omega_b = '-'
        if not found_1:
            omega_1 = '-'
        return([omega_b,omega_1,M0_vs_bfree,bneut_vs_bfree])

def read_site_model(site_file: str)-> list:
    with open(site_file, 'r') as site_log:
        site_log = site_log.readlines()
        found_modelM2 = False
        first_m2s = True
        first_m2c = True
        first_m2r = True

        found_modelM8 = False
        first_m8s = True
        first_m8c = True
        first_m8r = True

        found_modelM3 = False
        first_m3s = True
        first_m3c = True
        first_m3r = True

        found_modelM1 = False
        first_m1s = True
        first_m1c = True
        first_m1r = True

        found_modelM7 = False
        first_m7s = True
        first_m7c = True
        first_m7r = True

        m2_sites_c = ''
        m2_sites_s = ''
        m2_sites_r = ''

        m3_sites_c = ''
        m3_sites_s = ''
        m3_sites_r = ''

        m8_sites_c = ''
        m8_sites_s = ''
        m8_sites_r = ''

        m1_sites_c = ''
        m1_sites_s = ''
        m1_sites_r = ''

        m7_sites_c = ''
        m7_sites_s = ''
        m7_sites_r = ''

        m0_vs_m1 = '-'
        m0_vs_m2 = '-'
        m3_vs_m0 = '-'
        m0_vs_m7 = '-'
        m0_vs_m8 = '-'

        m1_vs_m2 = '-'
        m1_vs_m3 = '-'
        m1_vs_m8 = '-'

        m2_vs_m3 = '-'
        m7_vs_m2 = '-'
        m7_vs_m3 = '-'
        m7_vs_m8 = '-'
        m8_vs_m3 = '-'

        for line in site_log:

            if ('M1~' in line) and ('M2~' in line):
                m1_vs_m2 = line.split('|')[2][1:9]
            if ('M7~' in line) and ('M8~' in line):
                m7_vs_m8 = line.split('|')[2][1:9]
            if ('M0~' in line) and ('M3~' in line):
                m3_vs_m0 = line.split('|')[2][1:9]
            if ('M0~' in line) and ('M1~' in line):
                m0_vs_m1 = line.split('|')[2][1:9]
            if ('M0~' in line) and ('M2~' in line):
                m0_vs_m2 = line.split('|')[2][1:9]
            if ('M0~' in line) and ('M7~' in line):
                m0_vs_m7 = line.split('|')[2][1:9]
            if ('M0~' in line) and ('M8~' in line):
                m0_vs_m8 = line.split('|')[2][1:9]
            if ('M8~' in line) and ('M3~' in line):
                m8_vs_m3 = line.split('|')[2][1:9]
            if ('M1~' in line) and ('M3~' in line):
                m1_vs_m3 = line.split('|')[2][1:9]
            if ('M1~' in line) and ('M8~' in line):
                m1_vs_m8 = line.split('|')[2][1:9]
            if ('M2~' in line) and ('M3~' in line):
                m2_vs_m3 = line.split('|')[2][1:9]
            if ('M7~' in line) and ('M2~' in line):
                m7_vs_m2 = line.split('|')[2][1:9]
            if ('M7~' in line) and ('M3~' in line):
                m7_vs_m3 = line.split('|')[2][1:9]


            if 'Model M2' in line:
                found_modelM2 = True
                m2_sites_s = ''
                m2_sites_c = ''
                m2_sites_r = ''
            if ('Model M' in line) and ('Model M2' not in line):
                found_modelM2 = False
            if found_modelM2 and ('Positively-selected' in line):
                m2s_index = line.split('|')[0]
                m2s_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                
                if first_m2s:
                    m2_sites_s = m2s_index + m2s_prob
                    first_m2s = False
                else:
                    m2_sites_s = m2_sites_s + ';' + m2s_index + m2s_prob

            if found_modelM2 and ('Conserved' in line):
                m2c_index = line.split('|')[0]
                m2c_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m2c:
                    m2_sites_c = m2c_index+m2c_prob
                    first_m2c = False
                else:
                    m2_sites_c = m2_sites_c + ';' + m2c_index+m2c_prob


            if found_modelM2 and ('Relaxed' in line):
                m2r_index = line.split('|')[0]
                m2r_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m2r:
                    m2_sites_r = m2r_index+m2r_prob
                    first_m2r = False
                else:
                    m2_sites_r = m2_sites_r + ';' + m2r_index+m2r_prob



            if 'Model M8' in line:
                found_modelM8 = True
                m8_sites_s = ''
                m8_sites_c = ''
                m8_sites_r = ''
            if ('Model M' in line) and ('Model M8' not in line):
                found_modelM8 = False
            if found_modelM8 and ('Positively-selected' in line):
                m8s_index = line.split('|')[0]
                m8s_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m8s:
                    m8_sites_s = m8s_index+m8s_prob
                    first_m8s = False
                else:
                    m8_sites_s = m8_sites_s + ';' + m8s_index+m8s_prob

            if found_modelM8 and ('Conserved' in line):
                m8c_index = line.split('|')[0]
                m8c_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m8c:
                    m8_sites_c = m8c_index+m8c_prob
                    first_m8c = False
                else:
                    m8_sites_c = m8_sites_c + ';' + m8c_index+m8c_prob

            if found_modelM8 and ('Relaxed' in line):
                m8r_index = line.split('|')[0]
                m8r_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m8r:
                    m8_sites_r = m8r_index+m8r_prob
                    first_m8r = False
                else:
                    m8_sites_r = m8_sites_r + ';' + m8r_index+m8r_prob



            if 'Model M3' in line:
                found_modelM3 = True
                m3_sites_s = ''
                m3_sites_c = ''
                m3_sites_r = ''
            if ('Model M' in line) and ('Model M3' not in line):
                found_modelM3 = False

            if found_modelM3 and ('Positively-selected' in line):
                m3s_index = line.split('|')[0]
                m3s_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m3s:
                    m3_sites_s = m3s_index+m3s_prob
                    first_m3s = False
                else:
                    m3_sites_s = m3_sites_s + ';' + m3s_index+m3s_prob

            if found_modelM3 and ('Conserved' in line):
                m3c_index = line.split('|')[0]
                m3c_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m3c:
                    m3_sites_c = m3c_index+m3c_prob
                    first_m3c = False
                else:
                    m3_sites_c = m3_sites_c + ';' + m3c_index+m3c_prob

            if found_modelM3 and ('Relaxed' in line):
                m3r_index = line.split('|')[0]
                m3r_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m3c:
                    m3_sites_r = m3r_index+m3r_prob
                    first_m3r = False
                else:
                    m3_sites_r = m3_sites_r + ';' + m3r_index+m3r_prob


            if 'Model M1' in line:
                found_modelM1 = True
                m1_sites_s = ''
                m1_sites_c = ''
                m1_sites_r = ''
            if ('Model M' in line) and ('Model M1' not in line):
                found_modelM1 = False
            if found_modelM1 and ('Positively-selected' in line):
                m1s_index = line.split('|')[0]
                m1s_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                
                if first_m1s:
                    m1_sites_s = m1s_index + m1s_prob
                    first_m1s = False
                else:
                    m1_sites_s = m1_sites_s + ';' + m1s_index + m1s_prob

            if found_modelM1 and ('Conserved' in line):
                m1c_index = line.split('|')[0]
                m1c_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m1c:
                    m1_sites_c = m1c_index+m1c_prob
                    first_m1c = False
                else:
                    m1_sites_c = m1_sites_c + ';' + m1c_index+m1c_prob


            if found_modelM1 and ('Relaxed' in line):
                m1r_index = line.split('|')[0]
                m1r_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m1r:
                    m1_sites_r = m1r_index+m1r_prob
                    first_m1r = False
                else:
                    m1_sites_r = m1_sites_r + ';' + m1r_index+m1r_prob

            if 'Model M7' in line:
                found_modelM7 = True
                m7_sites_s = ''
                m7_sites_c = ''
                m7_sites_r = ''
            if ('Model M' in line) and ('Model M7' not in line):
                found_modelM7 = False
            if found_modelM7 and ('Positively-selected' in line):
                m7s_index = line.split('|')[0]
                m7s_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                
                if first_m7s:
                    m7_sites_s = m7s_index + m7s_prob
                    first_m7s = False
                else:
                    m7_sites_s = m7_sites_s + ';' + m7s_index + m7s_prob

            if found_modelM7 and ('Conserved' in line):
                m7c_index = line.split('|')[0]
                m7c_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m7c:
                    m7_sites_c = m7c_index+m7c_prob
                    first_m7c = False
                else:
                    m7_sites_c = m7_sites_c + ';' + m7c_index+m1c_prob


            if found_modelM7 and ('Relaxed' in line):
                m7r_index = line.split('|')[0]
                m7r_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_m7r:
                    m7_sites_r = m7r_index+m7r_prob
                    first_m7r = False
                else:
                    m7_sites_r = m7_sites_r + ';' + m7r_index+m7r_prob



        if m2_sites_c == '':
            m2_sites_c = '-'
        if m2_sites_s == '':
            m2_sites_s = '-'
        if m2_sites_r == '':
            m2_sites_r = '-'

        if m8_sites_c == '':
            m8_sites_c = '-'
        if m8_sites_s == '':
            m8_sites_s = '-'
        if m8_sites_r == '':
            m8_sites_r = '-'

        if m3_sites_c == '':
            m3_sites_c = '-'
        if m3_sites_s == '':
            m3_sites_s = '-'
        if m3_sites_r == '':
            m3_sites_r = '-'

        if m1_sites_c == '':
            m1_sites_c = '-'
        if m1_sites_s == '':
            m1_sites_s = '-'
        if m1_sites_r == '':
            m1_sites_r = '-'

        if m7_sites_c == '':
            m7_sites_c = '-'
        if m7_sites_s == '':
            m7_sites_s = '-'
        if m7_sites_r == '':
            m7_sites_r = '-'


        m2_sites_c = m2_sites_c.replace(' ', "")
        m2_sites_c = m2_sites_c.replace('\t', "")
        m2_sites_s = m2_sites_s.replace(' ', "")
        m2_sites_s = m2_sites_s.replace('\t', "")
        m2_sites_r = m2_sites_r.replace(' ', "")
        m2_sites_r = m2_sites_r.replace('\t', "")

        m3_sites_c = m3_sites_c.replace(' ', "")
        m3_sites_c = m3_sites_c.replace('\t', "")
        m3_sites_s = m3_sites_s.replace(' ', "")
        m3_sites_s = m3_sites_s.replace('\t', "")
        m3_sites_r = m3_sites_r.replace(' ', "")
        m3_sites_r = m3_sites_r.replace('\t', "")

        m8_sites_c = m8_sites_c.replace(' ', "")
        m8_sites_c = m8_sites_c.replace('\t', "")
        m8_sites_s = m8_sites_s.replace(' ', "")
        m8_sites_s = m8_sites_s.replace('\t', "")
        m8_sites_r = m8_sites_r.replace(' ', "")
        m8_sites_r = m8_sites_r.replace('\t', "")

        m1_sites_c = m1_sites_c.replace(' ', "")
        m1_sites_c = m1_sites_c.replace('\t', "")
        m1_sites_s = m1_sites_s.replace(' ', "")
        m1_sites_s = m1_sites_s.replace('\t', "")
        m1_sites_r = m1_sites_r.replace(' ', "")
        m1_sites_r = m1_sites_r.replace('\t', "")

        m7_sites_c = m7_sites_c.replace(' ', "")
        m7_sites_c = m7_sites_c.replace('\t', "")
        m7_sites_s = m7_sites_s.replace(' ', "")
        m7_sites_s = m7_sites_s.replace('\t', "")
        m7_sites_r = m7_sites_r.replace(' ', "")
        m7_sites_r = m7_sites_r.replace('\t', "")


    return([m0_vs_m1,m0_vs_m2, m3_vs_m0,m0_vs_m7, m0_vs_m8,m1_vs_m2,m1_vs_m3,m1_vs_m8,m2_vs_m3,m7_vs_m2, m7_vs_m3,m7_vs_m8, m8_vs_m3,
           m1_sites_s, m1_sites_c, m1_sites_r, m2_sites_s, m2_sites_c, m2_sites_r, 
           m3_sites_s, m3_sites_c, m3_sites_r, m7_sites_s, m7_sites_c, m7_sites_r, m8_sites_s, m8_sites_c, m8_sites_r])


def read_branch_site_model(branch_site_file: str)-> list:
    with open(branch_site_file, 'r') as branchsite_log:
        branchsite_log = branchsite_log.readlines()
        found_model_bsA = False
        first_mAs = True
        first_mAc = True
        first_mAr = True

        found_model_bsB = False
        first_mBs = True
        first_mBc = True
        first_mBr = True

        found_model_bsA1 = False
        first_mA1s = True
        first_mA1c = True
        first_mA1r = True

        mA_sites_c = ''
        mA_sites_s = ''
        mA_sites_r = ''

        mB_sites_c = ''
        mB_sites_s = ''
        mB_sites_r = ''

        mA1_sites_c = ''
        mA1_sites_s = ''
        mA1_sites_r = ''


        m0_vs_mA = '-'
        m0_vs_mA1 = '-'
        m0_vs_mB = '-'

        mA1_vs_mA = '-'
        mA1_vs_mB = '-'
        mA_vs_mB = '-'

        for line in branchsite_log:

            if ('bsA.' in line) and ('bsA1.' in line) and ('|' in line):
                mA1_vs_mA = line.split('|')[2][1:9]
            mA1_vs_mA = mA1_vs_mA.replace(' ', "")
            mA1_vs_mA = mA1_vs_mA.replace('\t', "")

            if ('bsB.' in line) and ('bsA1.' in line) and ('|' in line):
                mA1_vs_mB = line.split('|')[2][1:9]
            mA1_vs_mB = mA1_vs_mB.replace(' ', "")
            mA1_vs_mB = mA1_vs_mB.replace('\t', "")

            if ('bsA.' in line) and ('bsB.' in line) and ('|' in line):
                mA_vs_mB = line.split('|')[2][1:9]
            mA_vs_mB = mA_vs_mB.replace(' ', "")
            mA_vs_mB = mA_vs_mB.replace('\t', "")

            if ('M0~' in line) and ('bsA.' in line) and ('|' in line):
                m0_vs_mA = line.split('|')[2][1:9]
            m0_vs_mA = m0_vs_mA.replace(' ', "")
            m0_vs_mA = m0_vs_mA.replace('\t', "")

            if ('M0~' in line) and ('bsA1.' in line) and ('|' in line):
                m0_vs_mA1 = line.split('|')[2][1:9]
            m0_vs_mA1 = m0_vs_mA1.replace(' ', "")
            m0_vs_mA1 = m0_vs_mA1.replace('\t', "")

            if ('M0~' in line) and ('bsB.' in line) and ('|' in line):
                m0_vs_mB = line.split('|')[2][1:9]
            m0_vs_mB = m0_vs_mB.replace(' ', "")
            m0_vs_mB = m0_vs_mB.replace('\t', "")


            if 'Model bsA.' in line:
                found_model_bsA = True
                mA_sites_s = ''
                mA_sites_c = ''
                mA_sites_r = ''


            if ('Model bs' in line) and ('Model bsA.' not in line):
                found_model_bsA = False

            if found_model_bsA and ('Positively-selected' in line):
                 mAs_index = line.split('|')[0]

                 mAs_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                 if first_mAs:
                     mA_sites_s = mAs_index+mAs_prob
                     first_mAs = False
                 else:
                     mA_sites_s = mA_sites_s + ';' + mAs_index+mAs_prob

            if found_model_bsA and ('Conserved' in line):
                mAc_index = line.split('|')[0]
                mAc_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_mAc:
                    mA_sites_c = mAc_index+mAc_prob
                    first_mAc = False
                else:
                    mA_sites_c = mA_sites_c + ';' + mAc_index+mAc_prob

            if found_model_bsA and ('Relaxed' in line):
                print(line)
                mAr_index = line.split('|')[0]
                mAr_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_mAr:
                    mA_sites_r = mAr_index+mAr_prob
                    first_mAr = False
                else:
                    mA_sites_r = mA_sites_r + ';' + mAr_index+mAr_prob


            if 'Model bsB.' in line:
                found_model_bsB = True
                mB_sites_s = ''
                mB_sites_c = ''
                mB_sites_r = ''

            if ('Model bs' in line) and ('Model bsB.' not in line):
                found_model_bsB = False

            if found_model_bsB and ('Positively-selected' in line):
                mBs_index = line.split('|')[0]
                mBs_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_mBs:
                    mB_sites_s = mBs_index+mBs_prob
                    first_mBs = False
                else:
                    mB_sites_s = mB_sites_s + ';' + mBs_index+mBs_prob

            if found_model_bsB and ('Conserved' in line):
                mBc_index = line.split('|')[0]
                mBc_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_mBc:
                    mB_sites_c = mBc_index+mBc_prob
                    first_mBc = False
                else:
                    mB_sites_c = mB_sites_c + ';' +mBc_index + mBc_prob

            if found_model_bsB and ('Relaxed' in line):
                print(line)
                mBr_index = line.split('|')[0]
                mBr_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_mBc:
                    mB_sites_r = mBr_index+mBc_prob
                    first_mBr = False
                else:
                    mB_sites_r = mB_sites_r + ';' +mBr_index+ mBr_prob


            if 'Model bsA1.' in line:
                found_model_bsA1 = True
                mA1_sites_s = ''
                mA1_sites_c = ''
                mA1_sites_r = ''

            if ('Model bs' in line) and ('Model bsA1.' not in line):
                found_model_bsA1 = False

            if found_model_bsA1 and ('Positively-selected' in line):
                 mA1s_index = line.split('|')[0]

                 mA1s_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                 if first_mA1s:
                     mA1_sites_s = mA1s_index+mA1s_prob
                     first_mAs = False
                 else:
                     mA1_sites_s = mA1_sites_s + ';' + mA1s_index+mA1s_prob

            if found_model_bsA1 and ('Conserved' in line):
                mA1c_index = line.split('|')[0]
                mA1c_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_mA1c:
                    mA1_sites_c = mA1c_index+mA1c_prob
                    first_mA1c = False
                else:
                    mA1_sites_c = mA1_sites_c + ';' + mA1c_index+mA1c_prob

            if found_model_bsA1 and ('Relaxed' in line):
                mA1r_index = line.split('|')[0]
                mA1r_prob = '(' + line.split('|')[1].split()[3].replace('probability > ','')
                if first_mA1r:
                    mA1_sites_r = mA1r_index+mA1r_prob
                    first_mA1r = False
                else:
                    mA1_sites_r = mA1_sites_r + ';' + mA1r_index+mA1r_prob


        if mA_sites_c == '':
            mA_sites_c = '-'
        if mA_sites_s == '':
            mA_sites_s = '-'
        if mA_sites_r == '':
            mA_sites_r = '-'

        if mB_sites_c == '':
            mB_sites_c = '-'
        if mB_sites_s == '':
            mB_sites_s = '-'
        if mB_sites_r == '':
            mB_sites_r = '-'

        if mA1_sites_c == '':
            mA1_sites_c = '-'
        if mA1_sites_s == '':
            mA1_sites_s = '-'
        if mA1_sites_r == '':
            mA1_sites_r = '-'

        mA_sites_c = mA_sites_c.replace(' ', "")
        mA_sites_c = mA_sites_c.replace('\t', "")
        mA_sites_s = mA_sites_s.replace(' ', "")
        mA_sites_s = mA_sites_s.replace('\t', "")
        mA_sites_r = mA_sites_r.replace(' ', "")
        mA_sites_r = mA_sites_r.replace('\t', "")

        mA1_sites_c = mA1_sites_c.replace(' ', "")
        mA1_sites_c = mA1_sites_c.replace('\t', "")
        mA1_sites_s = mA1_sites_s.replace(' ', "")
        mA1_sites_s = mA1_sites_s.replace('\t', "")
        mA1_sites_r = mA1_sites_r.replace(' ', "")
        mA1_sites_r = mA1_sites_r.replace('\t', "")

        mB_sites_c = mB_sites_c.replace(' ', "")
        mB_sites_c = mB_sites_c.replace('\t', "")
        mB_sites_s = mB_sites_s.replace(' ', "")
        mB_sites_s = mB_sites_s.replace('\t', "")
        mB_sites_r = mB_sites_r.replace(' ', "")
        mB_sites_r = mB_sites_r.replace('\t', "")


    return([m0_vs_mA, m0_vs_mA1,m0_vs_mB, mA1_vs_mA, mA_vs_mB, mA1_vs_mB, mA_sites_s, mA_sites_c, mA_sites_r, 
           mA1_sites_s, mA1_sites_c, mA1_sites_r, mB_sites_s,mB_sites_c, mB_sites_r])



def return_domain_names(domain: str) -> list:
    """
        Returns domain names for a recombination event either for domains exposed to recombination or the intact ones
    """
    rtype = domain
    if rtype.startswith('domain1'):
        return ['domain2', 'domain3']
    elif rtype.startswith('domain2'):
        return ['domain1','domain3']
    elif rtype.startswith('domain3'):
        return ['domain1','domain2']
    else:
        return None

def get_LnL_res(event_dir, tox_type):
    model_type_dir = os.path.join(os.path.realpath(event_dir),tox_type)
    if 'branch_site_model' in tox_type:
        model_dict = {'M0':-1, 'bsA':-1, 'bsA1':-1, 'bsB':-1}
    elif 'branch_model' in tox_type:
        model_dict = {'M0':-1, 'b_free':-1, 'b_neut':-1}
    elif 'site_model' in tox_type:
        model_dict = {'M0':-1,'M1':-1,'M2':-1,'M3':-1,'M7':-1,'M8':-1}

    if not os.path.isdir(model_type_dir):
        return(['-']*len(model_dict))

    for subdir in os.listdir(model_type_dir):
        model_ind = subdir.split('~')[0].split('.')[0]

        out_file = os.path.join(os.path.realpath(model_type_dir),os.path.join(subdir,'out'))
        with open(out_file, 'r') as model_res:
            for line in model_res:
                if 'lnL' in line:
                    ln_res = line.split()[4]
                    model_dict[model_ind] = ln_res

    return(list(model_dict.values()))


def read_log_files(row, evol_dir, event_type):
    event_dir=os.path.join(os.path.realpath(args.evol_dir),event_type+str(row[1]['ID']))


    if os.path.isdir(event_dir):

        branch_rec=os.path.join(os.path.realpath(event_dir),'branch_model_rec.log')
        branch_par=os.path.join(os.path.realpath(event_dir),'branch_model_par.log')

        branch_site_rec=os.path.join(os.path.realpath(event_dir),'branch_site_model_rec.log')
        branch_site_par=os.path.join(os.path.realpath(event_dir),'branch_site_model_par.log')

        site=os.path.join(os.path.realpath(event_dir),'site_model.log')

        if os.path.isfile(branch_rec):
            rec_branch=read_branch_model(branch_rec) + get_LnL_res(event_dir, 'branch_model_rec')
        else:
            rec_branch=['-']*7

        if os.path.isfile(branch_par):
            par_branch=read_branch_model(branch_par) + get_LnL_res(event_dir, 'branch_model_par')
        else:
            par_branch=['-']*7

        if os.path.isfile(branch_site_rec):
            rec_branch_site=read_branch_site_model(branch_site_rec) + get_LnL_res(event_dir, 'branch_site_model_rec')
        else:
            rec_branch_site=['-']*19

        if os.path.isfile(branch_site_par):
            par_branch_site=read_branch_site_model(branch_site_par) + get_LnL_res(event_dir, 'branch_site_model_par')

        else:
            par_branch_site=['-']*19
       
        if os.path.isfile(site):
            site_res=read_site_model(site) + get_LnL_res(event_dir, 'site_model')
        else:
            site_res=['-']*34

        return(rec_branch+site_res+rec_branch_site+par_branch+par_branch_site)

    else:
        return(None)
    

            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarises evolutionary selection results for all recombination events')
    parser.add_argument('-r', '--rec', dest='rec_tab', help='the path to the preparatory table for selection evaluation',
                        type=str)
    parser.add_argument('-e', '--evol', dest='evol_dir', help='the path to the selection evaluation results',
                        type=str)


    args = parser.parse_args()

    rec_table=read_csv_to_list(args.rec_tab, headless=False)
    rec_df = pd.DataFrame(rec_table[1:], columns=rec_table[0])
    head_row = ['Id', 'Rec_flag','Rec', 'Min_par', 'Maj_par1','Maj_par2', 'Type','Current_domain' ,'len'] 

    dummy_row = ['omega_b','omega_1','M0_vs_bfree','bneut_vs_bfree', 'M0_lnL', 'b_bfree_L', 'b_neut_L', 
                'm0_vs_m1','m0_vs_m2', 'm3_vs_m0','m0_vs_m7', 'm0_vs_m8','m1_vs_m2','m1_vs_m3','m1_vs_m8','m2_vs_m3','m7_vs_m2', 'm7_vs_m3','m7_vs_m8', 'm8_vs_m3',
                'm1_sites_s', 'm1_sites_c', 'm1_sites_r', 'm2_sites_s', 'm2_sites_c', 'm2_sites_r', 
                'm3_sites_s', 'm3_sites_c', 'm3_sites_r', 'm7_sites_s', 'm7_sites_c', 'm7_sites_r', 'm8_sites_s', 'm8_sites_c', 'm8_sites_r', 'M0_lnL', 'M1_lnL', 'M2_lnL', 'M3_lnL', 'M7_lnL', 'M8_lnL',
                'm0_vs_mA', 'm0_vs_mA1', 'm0_vs_mB', 'mA1_vs_mA', 'mA_vs_mB', 'mA1_vs_mB', 'mA_sites_s', 'mA_sites_c', 'mA_sites_r', 
                'mA1_sites_s', 'mA1_sites_c', 'mA1_sites_r', 'mB_sites_s','mB_sites_c', 'mB_sites_r', 'M0_lnL', 'bsA_lnL', 'bsA1_lnL', 'bsB_lnL']

    row_rec = ['rec_'+el for el in dummy_row]
    row_par = ['par_'+el for el in dummy_row[0:7]+dummy_row[41:]]
    head_row = head_row + row_rec + row_par


    event_rows=[]
    event_rows.append(head_row)
    for row in rec_df.iterrows():
        non_rec_domains=return_domain_names(row[1]['Type'])
        min_par_row=[row[1]['ID'],'Rec', row[1]['Rec_min'],row[1]['Min_par'],row[1]['Maj1_par'],row[1]['Maj2_par'],row[1]['Type'],row[1]['Type'],row[1]['Min_evol']]
        min_log_results=read_log_files(row, args.evol_dir,'minor_event_')

        if min_log_results:
            event_rows.append(min_par_row+min_log_results)

        maj1_par_row=[row[1]['ID'],'Maj1', row[1]['Rec_min'],row[1]['Min_par'],row[1]['Maj1_par'],row[1]['Maj2_par'],row[1]['Type'],non_rec_domains[0],row[1]['Maj1_evol']]
        maj1_log_results=read_log_files(row, args.evol_dir,'major1_event_')

        if maj1_log_results:
            event_rows.append(maj1_par_row+maj1_log_results)

        maj2_par_row=[row[1]['ID'],'Maj2', row[1]['Rec_min'],row[1]['Min_par'],row[1]['Maj1_par'],row[1]['Maj2_par'],row[1]['Type'],non_rec_domains[1],row[1]['Maj2_evol']]
        maj2_log_results=read_log_files(row, args.evol_dir,'major2_event_')

        if maj2_log_results:
            event_rows.append(maj2_par_row+maj2_log_results)
        
    write_csv(event_rows,'selection_results_summary_with_relaxed.csv')
         



