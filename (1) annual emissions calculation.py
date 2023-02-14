# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 10:46:59 2022

@author: Daria
"""

import os
import pandas as pd
import numpy as np

# In[1]: define the directories

input_directory_path = os.path.join('input')
output_directory_path = os.path.join('output')

os.makedirs(input_directory_path, exist_ok=True)
os.makedirs(output_directory_path, exist_ok=True)

# In[2]: define functions, e.g. for uploading raw data (Energy Balances and Emission factors) and modifying dataframes

def load_nrg_bal(path, fn):
    df = pd.read_csv(os.path.join(path, fn), sep = ',', header = 0) #the 1st line contains header, the 1st row does not contain indexes
    
    #rename long names of the column "nrg_bal" content the original file   
    df.replace({
            'Transformation input - electricity and heat generation - main activity producer electricity only - energy use': 'TI_EHG_MAPE_E',
            'Transformation input - electricity and heat generation - main activity producer combined heat and power - energy use': 'TI_EHG_MAPCHP_E', 
            'Transformation input - electricity and heat generation - autoproducer electricity only - energy use': 'TI_EHG_APE_E', 
            'Transformation input - electricity and heat generation - autoproducer combined heat and power - energy use': 'TI_EHG_APCHP_E',
            'Energy sector - electricity and heat generation - energy use':'NRG_EHG_E',
            'Distribution losses': 'DL',
            'Gross electricity production': 'GEP',
            'Gross heat production': 'GHP',
            'Gross heat production - main activity producer combined heat and power': 'GHP_MAPCHP',
            'Gross heat production - autoproducer combined heat and power': 'GHP_APCHP',
            'Germany (until 1990 former territory of the FRG)': 'Germany'}, inplace = True) #inplace = false - dos not overwrite the existing dataframe (returns copy of the initial string)
    return df

def load_SIEC_codes(path, fn):
    df = pd.read_csv(os.path.join(path, fn), sep = r"[", header = None, names = ['SIEC', 'Code'])
    df.loc[:,"SIEC"] = df.loc[:,"SIEC"].str[:-1] #remove " " at the end
    df.loc[:,"Code"] = df.loc[:,"Code"].str[:-1] #remove "]"
    return df

def match_nrg_bal_SIEC(path, fn):
    df = pd.read_csv(os.path.join(path, fn))
    return df

def load_factors(path, fn):
    df = pd.read_csv(os.path.join(path, fn), sep = ';')# decimal = ',')
    df = df.replace({np.nan: 0}).replace('-', 0.)
    col_to_numeric = df[['direct [kg CO2/TJ]', 'direct [kg CH4/TJ]', 'direct [kg N2O/TJ]']] #'upstream [kg CO2/kWh]']]
    for col in col_to_numeric:
        df[col] = [float(str(i).replace(",",".")) for i in df[col]]
    col_to_kWh = df[['direct [kg CO2/TJ]', 'direct [kg CH4/TJ]', 'direct [kg N2O/TJ]']]
    for col in col_to_kWh:
        df[col] = [i * 3.6 / 1e6 for i in df[col]]#convert TJ values to kWh. 1 kWh = 3.6e-6 TJ
    df = df.rename({'direct [kg CO2/TJ]':'cf_CO2_d',
                    'direct [kg CH4/TJ]':'cf_CH4_d',
                    'direct [kg N2O/TJ]':'cf_N2O_d',}, axis = 1)#the upstream values are in CO2e, not CO2 (error in the csv.) 'upstream [kg CO2/kWh]':'cf_CO2e_up'
    df['cf_CO2e_d'] = df["cf_CO2_d"] + df['cf_CH4_d']*28 + df["cf_N2O_d"]*256 #IPCC AR5 default values for converting GHG to CO2e [kg CO2e/kWh]
    return df
    
def extend_SIEC(df, df_codes):
    SIEC_codes_dict = df_codes.set_index('SIEC').to_dict()['Code']
    df_new = df.copy()
    df_new['Code'] = df_new['SIEC'].map(SIEC_codes_dict)
    SIEC_col_index = df_new.columns.get_loc("SIEC") #returns index of the newly added column "SIEC"
    col_Code = df_new.pop("Code")
    df_new.insert(SIEC_col_index,"Code",col_Code)
    del SIEC_col_index
    del col_Code
    return df_new
    
def aggregate_fuels(df_total,col_heading, fuel_str, fuel_type):
    df_agg = df_total.loc[df_total[col_heading].str.fullmatch(fuel_str)]
    df_agg = df_agg.groupby(["TIME", "GEO"], as_index = False).sum()
    df_agg['Fuel_type'] = fuel_type
    #place the new column with the fuel type after 'GEO' 
    get_col_index = df_agg.columns.get_loc('GEO')
    add_col_new = df_agg.pop('Fuel_type')
    df_agg.insert(get_col_index + 1,'Fuel_type',add_col_new)
    return df_agg

def extract_annual_emissions(annual_emissions):
    # create an empty dictionary
    dict_result = {}
    
    # derive unique area and year values
    areas = annual_emissions['geo'].drop_duplicates().to_list()
    years = annual_emissions['TIME_PERIOD'].drop_duplicates().to_list()
    
    # create a list of keys for each country and year
    keys = ['annual_emissions_' + area + '_' + str(year) for area in areas for year in years]
    
    for key in keys:
        
        # derive area and year
        area = str(key)[len('annual_emissions_') : -len('_20XX') :]
        year = int(str(key)[-len('20XX') : :])
        
        # derive a dataframe corresponding to area and year
        df = annual_emissions[(annual_emissions['geo'] == area) & (annual_emissions['TIME_PERIOD'] == year)]
        
        # insert the dataframe in the dictionary
        dict_result[key] = df
    return(dict_result)

# In[3]: download raw data as DataFrames
    
# load energy balance data in GWh from file
# required parameters are listed in Table S1 of Supplementary Information
nrg_bal = load_nrg_bal(input_directory_path, 'nrg_bal_c_v_14_04_2022.csv') 

#load SIEC codes from file
SIEC_codes = load_SIEC_codes(input_directory_path, 'SIEC_code_label.csv') 

#add a column with SIEC codes
nrg_bal = extend_SIEC(nrg_bal,SIEC_codes) 
factors = load_factors(input_directory_path, 'multipl_factors.csv')

# In[4]:filter out dataframes with parameter-specific data from Energy Balances (nrg_bal). 

# TI - transformation input
# AP - autoproducers
# MAP - main activity producers
# CHP - combined heat and power
# ele - electricity

# Procedure for combining data on MAP and AP:
#   1) define string containing nrg_bal parameters to filter (refers to TI and GHP)
#   2) derive a dataframe matching the string and replace ":" and ":z" by zeros in the column 'OBS_VALUE'
#   3) convert string to a float and remove empty spaces in the column 'OBS_VALUE' (necessary for the old version of nrg_bal)
#   4) calculate the sum of parameters AP and MAP within each dataframe and remove unrelevant columns
#   5) delete the string (refers to TI and GHP) / reset index in the newly created dataframe

# TI to ele plants from AP and MAP
ti_ele_plants_str = 'TI_EHG_MAPE_E|TI_EHG_APE_E'
ti_ele_plants = nrg_bal.loc[nrg_bal['NRG_BAL'].str.contains(ti_ele_plants_str, regex = True)].replace(':',0.).replace(': z',0.) 
ti_ele_plants['Value'] = [float(str(i).replace(' ','')) for i in ti_ele_plants['Value']] 
ti_ele_plants = ti_ele_plants.groupby(['TIME','GEO','Code', 'SIEC'],as_index = False).sum(['Value'])
del ti_ele_plants_str

# eurostat 2023 update verison
# ti_ele_plants_str = 'TI_EHG_MAPE_E|TI_EHG_APE_E'
# ti_ele_plants = nrg_bal.loc[nrg_bal['nrg_bal'].str.contains(ti_ele_plants_str, regex = True)].replace(':',0.).replace(': z',0.) 
# ti_ele_plants['OBS_VALUE'] = [float(str(i).replace(" ","")) for i in ti_ele_plants['OBS_VALUE']]
# ti_ele_plants = ti_ele_plants.groupby(['TIME_PERIOD', 'geo', 'siec'], as_index = False).sum(['OBS_VALUE'])
# del ti_ele_plants_str

# TI to CHP plants for heat and ele generation from AP and MAP
ti_CHP_str = 'TI_EHG_APCHP_E|TI_EHG_MAPCHP_E'
ti_CHP_plants = nrg_bal.loc[nrg_bal['NRG_BAL'].str.contains(ti_CHP_str, regex = True)].replace(':',0.).replace(': z',0.)
ti_CHP_plants['Value'] = [float(str(i).replace(' ','')) for i in ti_CHP_plants['Value']]
ti_CHP_plants = ti_CHP_plants.groupby(["TIME","GEO","Code", "SIEC"], as_index = False).sum(['Value'])
del ti_CHP_str 

# eurostat 2023 update verison
# ti_CHP_str = 'TI_EHG_APCHP_E|TI_EHG_MAPCHP_E'
# ti_CHP_plants = nrg_bal.loc[nrg_bal['nrg_bal'].str.contains(ti_CHP_str, regex = True)].replace(':',0.).replace(': z',0.)
# ti_CHP_plants['OBS_VALUE'] = [float(str(i).replace(" ","")) for i in ti_CHP_plants['OBS_VALUE']]
# ti_CHP_plants = ti_CHP_plants.groupby(['TIME_PERIOD', 'geo', 'siec'], as_index = False).sum(['OBS_VALUE'])
# del ti_CHP_str 

# group gross heat output of CHP plants from AP and MAP
heat_output_CHP_str='GHP_MAPCHP|GHP_APCHP'
heat_output_CHP = nrg_bal.loc[nrg_bal['NRG_BAL'].str.contains(heat_output_CHP_str, regex = True)].replace(':',0.).replace(': z',0.)
heat_output_CHP['Value'] = [float(str(i).replace(' ','')) for i in heat_output_CHP["Value"]]
heat_output_CHP = heat_output_CHP.groupby(['TIME', 'GEO', 'Code', 'SIEC'],as_index=False).sum(['Value'])
del heat_output_CHP_str

# eurostat 2023 update version
# heat_output_CHP_str='GHP_MAPCHP|GHP_APCHP'
# heat_output_CHP = nrg_bal.loc[nrg_bal['nrg_bal'].str.contains(heat_output_CHP_str, regex = True)].replace(':',0.).replace(': z',0.)
# heat_output_CHP['OBS_VALUE'] = [float(str(i).replace(" ","")) for i in heat_output_CHP['OBS_VALUE']]
# heat_output_CHP = heat_output_CHP.groupby(['TIME_PERIOD', 'geo', 'siec'],as_index=False).sum(['OBS_VALUE'])
# del heat_output_CHP_str

# gross ele production, from ele and CHP plants, from AP anf MAP
GEP = nrg_bal.loc[nrg_bal['NRG_BAL'].str.fullmatch('GEP')].replace(':',0.).replace(': z',0.)
GEP['Value'] = [float(str(i).replace(' ','')) for i in GEP['Value']]
GEP_total = GEP.loc[GEP['SIEC'].str.fullmatch('Total')]
GEP_total = GEP_total.drop(['NRG_BAL', 'Code', 'SIEC', 'UNIT'],axis = 1)
GEP_total.reset_index(drop = True, inplace = True)

# eurostat 2023 update version
# GEP = nrg_bal.loc[nrg_bal['nrg_bal'].str.fullmatch('GEP')].replace(':',0.).replace(': z',0.)
# GEP['OBS_VALUE'] = [float(str(i).replace(' ','')) for i in GEP['OBS_VALUE']]
# GEP_total = GEP.loc[GEP['siec'].str.fullmatch('TOTAL')]
# GEP_total = GEP_total[['TIME_PERIOD', 'geo', 'OBS_VALUE']]
# GEP_total.reset_index(drop = True, inplace = True)

# gross heat production, from AP and MAP
GHP = nrg_bal.loc[nrg_bal['NRG_BAL'].str.fullmatch('GHP')].replace(':',0.).replace(': z',0.)
GHP['Value'] = [float(str(i).replace(' ','')) for i in GHP['Value']]
GHP.reset_index(drop = True, inplace = True)
GHP_total = GHP.loc[GHP["SIEC"].str.fullmatch('Total')]
GHP_total = GHP_total.drop(['NRG_BAL', 'Code', 'SIEC', 'UNIT'],axis = 1)
GHP_total.reset_index(drop = True, inplace = True)

# eurostat 2023 update version
# GHP = nrg_bal.loc[nrg_bal['nrg_bal'].str.fullmatch('GHP')].replace(':',0.).replace(': z',0.)
# GHP['OBS_VALUE'] = [float(str(i).replace(" ","")) for i in GHP['OBS_VALUE']]
# GHP_total = GHP.loc[GHP['siec'].str.fullmatch('TOTAL')]
# GHP_total = GHP_total[['TIME_PERIOD', 'geo', 'OBS_VALUE']]
# GHP_total.reset_index(drop = True, inplace = True)

# part of elecreicity directly utilzed in ele and CHP plants for ele generation
own_use_total = nrg_bal.loc[nrg_bal['NRG_BAL'].str.contains('NRG_EHG_E', regex = True)].replace(':',0.).replace(': z',0.)
own_use_total['Value'] = [float(str(i).replace(' ','')) for i in own_use_total['Value']]
own_use_total = own_use_total.loc[own_use_total['SIEC'].str.contains('Electricity')]
own_use_total = own_use_total.drop(['NRG_BAL', 'Code', 'SIEC', 'UNIT'],axis = 1)
own_use_total.reset_index(drop = True, inplace = True)

# eurostat 2023 update version
# own_use_total = nrg_bal.loc[nrg_bal['nrg_bal'].str.contains('NRG_EHG_E', regex = True)].replace(':',0.).replace(': z',0.)
# own_use_total['OBS_VALUE'] = [float(str(i).replace(" ","")) for i in own_use_total['OBS_VALUE']]
# own_use_total = own_use_total.loc[own_use_total['siec'].str.contains('E700')] # there is no division of own use by fuel types -> E700 ('electricity')
# own_use_total = own_use_total[['TIME_PERIOD', 'geo', 'OBS_VALUE']]
# own_use_total.reset_index(drop = True, inplace = True)

# In[5]: create a dataframe with data on transformation input to CHP plants for electricity generation only.
    
#calulate the transformation input (TI) to CHP plants for electricity gneration only
#merge two dataframes - ti_CHP_plants nd heat_output_CHP
ti_CHP_plants_ele = pd.merge(ti_CHP_plants, heat_output_CHP, on = ['TIME', 'GEO', 'Code', 'SIEC'], how = 'inner', suffixes = ('_in', '_out')) # for eurostat 2023 update version: on = ['TIME_PERIOD', 'geo', 'siec']

# calculate TI to CHP plants for ele production using fixed heat approach (IEA) with a constant efficiency of 90%
# the method does not work with Finland and Norway (see explanation in IPCC) --> generated values for these countries are incorrect
ti_CHP_plants_ele["Value"] = ti_CHP_plants_ele["Value_in"] - (ti_CHP_plants_ele['Value_out']/0.9) 

# remove two columns with the values from the two orginal dataframes
ti_CHP_plants_ele = ti_CHP_plants_ele.drop(columns = ['Value_in', 'Value_out'])

# In[6]: derive the average annual share of electricity consumed by all plants (ele and CHP) for electricity generation
# own use share = % of internally consumed ele from the total ele generated (refers to ele only, not ele and heat!)

# merge the 'ti_CHP_plants' and 'heat_output_CHP' dataframes
dfs = [GEP_total, GHP_total, own_use_total]
from functools import reduce
own_use_ele = reduce(lambda left, right: pd.merge(left, right, on = ["TIME", "GEO"]), dfs)
del dfs

# derive the share of electricity own use allocated to electricity production 
own_use_ele["Value"] = own_use_ele["Value"] / (own_use_ele["Value_y"] + own_use_ele["Value_x"])
own_use_ele = own_use_ele.drop(columns = ["Value_x", "Value_y"])

# In[7]: create a separate dataframe for collecting data relevant for the EF calculation

input_total = ti_ele_plants.rename({'Value': 'ti_ele'}, axis = 1)

# add a column with fuel-specific parameters to the 'input_total' dataframe
# Procedure:
#   1) create an iterator for checking 3 columns per row
#   2) extract data filtered by the iterator in a separate dataframe ('df_...') in a list
#   3) convert list of values to 1d array 
#   4) delete the intermediate list

# add a column with fuel-specific GEP
iterator = zip(input_total['TIME'], input_total['GEO'], input_total['Code'])
df_GEP = [GEP[(GEP['TIME'] == time) & (GEP['GEO'] == geo) & (GEP['Code'] == code)].Value for time, geo, code in iterator]
input_total['GEP'] = np.array(df_GEP)[:,0] #convert to 1d array
del df_GEP
    
# add a column with own use for ele generation
iterator = zip(input_total['GEP'], input_total['TIME'], input_total['GEO'])#create an iterator from 3 columns per row
df_own_use = [GEP * own_use_ele[(own_use_ele['TIME'] == time) & (own_use_ele['GEO'] == geo)].Value for GEP,time,geo in iterator]
input_total['own_use_ele'] = np.array(df_own_use)[:,0] 
del df_own_use

# add a column with fuel-specific GHP
iterator = zip(input_total['TIME'], input_total['GEO'], input_total['Code'])
df_GHP = [GHP[(GHP['TIME'] == time) & (GHP['GEO'] == geo) & (GHP['Code'] == code)].Value for time, geo, code in iterator]
input_total['GHP'] = np.array(df_GHP)[:,0] 
del df_GHP

# add a column with fuel-specific input to CHP
iterator = zip(input_total['TIME'], input_total['GEO'], input_total['Code'])
df_CHP_ele = [ti_CHP_plants_ele[(ti_CHP_plants_ele['TIME'] == time) & (ti_CHP_plants_ele['GEO'] == geo) & (ti_CHP_plants_ele['Code'] == code)].Value for time, geo, code in iterator]
input_total['ti_CHP'] = np.array(df_CHP_ele)[:,0] 
del df_CHP_ele

# calculate total electricity output as sum of ti_ele_plants, own_use_ele, and ti_CHP_plants_ele
col_input_total = ['ti_ele', 'own_use_ele', 'ti_CHP']
input_total['ti_ele_total'] = input_total[col_input_total].sum(axis = 1)
del col_input_total

# In[8]: calculate annual combustion emissions per SIEC fuel type

#calculate annual direct combustion emissions
annual_emissions = pd.merge(input_total, factors, on = ['Code', 'SIEC'])

# eurostat 2023 update version
# rename column headings in the 'factors' dataframe to comply with the updated data structure
# factors = factors.rename(columns = {'Code': 'siec', 'SIEC': 'label'}) 

# annual_emissions = pd.merge(input_total, factors, on = 'siec')

# # change the position of the 'labels' column
# index = annual_emissions.columns.get_loc('siec')
# col = annual_emissions.pop('label')
# annual_emissions.insert(index, 'siec_label', col)

annual_emissions['Mt_CO2_d'] = annual_emissions['ti_ele_total'] * annual_emissions['cf_CO2_d'] / 1e3
annual_emissions['Mt_CH4_d'] = annual_emissions['ti_ele_total'] * annual_emissions['cf_CH4_d'] / 1e3
annual_emissions['Mt_N2O_d'] = annual_emissions['ti_ele_total'] * annual_emissions['cf_N2O_d'] / 1e3
annual_emissions['Mt_CO2e_d'] = annual_emissions['ti_ele_total'] * annual_emissions['cf_CO2e_d'] / 1e3

annual_emissions = annual_emissions.drop(columns = ['GHP','cf_CO2_d', 'cf_CH4_d', 'cf_N2O_d', 'cf_CO2e_d'])

# extract dataframes for each year and country (eurostat 2023 update version)
# dict_annual_emissions_test = extract_annual_emissions(annual_emissions)

# In[9]: aggregate annual emissions into larger SIEC groups and ENTSO-E generation types 
# eurostat 2023 update version: use 'siec' instead of 'Code'

hard_coal = 'C0110|C0121|C0129'
an_em_hard_coal = aggregate_fuels(annual_emissions, 'Code', hard_coal, 'hard_coal')

brown_coal = 'C0210|C0220'
an_em_brown_coal = aggregate_fuels(annual_emissions, 'Code', brown_coal, 'brown_coal')

coal_products = 'C0311|C0312|C0313|C01314|C0320|C0330|C0340'
an_em_coal_products = aggregate_fuels(annual_emissions, 'Code',coal_products, 'coal_products')

coal_derived_gas = 'C0350|C0360|C0371|C0372|C0379|C0390'
an_em_coal_derived_gas = aggregate_fuels(annual_emissions, 'Code', coal_derived_gas, 'coal_derived_gas')

peat = 'P1100|P1200'
an_em_peat = aggregate_fuels(annual_emissions, 'Code', peat, 'peat')

oil_shale = 'S2000'
an_em_oil_shale = aggregate_fuels(annual_emissions, 'Code', oil_shale, 'oil_shale')

oil = 'O4410|O4200|O4300|O4400|O4500|O4610|O4620|O4630|O4640|O4651XR5210B|O4652|O4653|O4661XR5230B|O4669|O4671XR5220B|O4672|O4680|O4691|O4692|O4693|O4694|O4695|O4699'
an_em_oil = aggregate_fuels(annual_emissions, 'Code', oil, 'oil')

natural_gas = 'G3000'
an_em_nat_gas = aggregate_fuels(annual_emissions, 'Code', natural_gas, 'natural_gas')

nonR_waste = 'W6100|W6220'
an_em_nonR_waste = aggregate_fuels(annual_emissions, 'Code', nonR_waste, 'nonR_waste')

biofuels = 'R5110-5150_W6000RI|W6210|R5300|R5210P|R5210B|R5220P|R5220B|R5230P|R5230B|R5290'
an_em_biofuels = aggregate_fuels(annual_emissions, 'Code', biofuels, 'biofuels')

an_em_total = aggregate_fuels(annual_emissions, 'Code', 'TOTAL', 'total') # total refers to electricity only, but has no emissions

# derive three additional categories following the ENTSO-E structure (split coal_products between hard and brown coal)
hard_coal_ENTSOE = 'C0110|C0121|C0129|C0311|C0312|C0320|C0340'
an_em_hard_coal_ENTSOE = aggregate_fuels(annual_emissions, 'Code', hard_coal_ENTSOE, 'hard_coal')

brown_coal_ENTSOE = 'C0210|C0220|C0330'
an_em_brown_coal_ENTSOE = aggregate_fuels(annual_emissions, 'Code', brown_coal_ENTSOE, 'brown_coal')

fossil_gas_ENTSOE = 'C0350|C0360|C0371|C0372|C0379|C0390|G3000'
an_em_fossil_gas_ENTSOE = aggregate_fuels(annual_emissions, 'Code', fossil_gas_ENTSOE, 'fossil_gas')

# In[10]: create two separate tables with data aggregated by fuel types and ENTSO-E generation types

fuel_type_agg = [an_em_total,
                 an_em_hard_coal, 
                 an_em_brown_coal,
                 an_em_coal_derived_gas,
                 an_em_coal_products,
                 an_em_peat,
                 an_em_oil_shale,
                 an_em_oil,
                 an_em_nat_gas,
                 an_em_nonR_waste,
                 an_em_biofuels]

gen_type_agg = [an_em_total,
                an_em_hard_coal_ENTSOE, 
                an_em_brown_coal_ENTSOE,
                an_em_fossil_gas_ENTSOE,
                an_em_peat,                                     
                an_em_oil_shale,
                an_em_oil,
                an_em_nonR_waste,
                an_em_biofuels]

annual_data_SIEC = pd.concat(fuel_type_agg, axis = 0, ignore_index = True) # emissions columns not filled yet, contains electricity data only

annual_data_ENTSOE = pd.concat(gen_type_agg, axis = 0, ignore_index = True) # emissions columns not filled yet, contains electricity data only
annual_data_ENTSOE = annual_data_ENTSOE.replace({'biofuels': 'biomass'}).rename(columns = {'Fuel_type': 'Generation_type'}) 


del [an_em_total,
     an_em_hard_coal, 
     an_em_brown_coal,
     an_em_coal_derived_gas,
     an_em_peat,
     an_em_oil_shale,
     an_em_oil,
     an_em_nat_gas,
     an_em_nonR_waste,
     an_em_biofuels,
     an_em_hard_coal_ENTSOE, 
     an_em_brown_coal_ENTSOE,
     an_em_fossil_gas_ENTSOE]

del fuel_type_agg
del gen_type_agg

# In[11]: calculate annual emissions for the SIEC fuel types and ENTSO-E generation types

# extract dataframe with the total annual electricity generation
# columns with emissions remain empty
annual_total = annual_data_SIEC.loc[annual_data_SIEC['Fuel_type'].str.fullmatch('total')]

# remove the "total" rows
annual_data_SIEC = annual_data_SIEC.drop(annual_data_SIEC[annual_data_SIEC['Fuel_type'] == 'total'].index)
annual_data_ENTSOE = annual_data_ENTSOE.drop(annual_data_ENTSOE[annual_data_ENTSOE["Generation_type"] == 'total'].index)

# aggregate coal-related fuel types into one group
coal_total = 'hard_coal|brown_coal|coal_products|coal_derived_gas|peat|oil_shale'
annual_data_SIEC_coal = aggregate_fuels(annual_data_SIEC, 'Fuel_type', coal_total, 'coal_total')
annual_data_SIEC = pd.concat([annual_data_SIEC, annual_data_SIEC_coal], axis = 0, ignore_index = False) # add coal_total for validation purposes

# define all fuel types relevant for the emissions calculation
fossil_fuels_str = 'coal_total|oil|natural_gas|nonR_waste' 

# derive a dataframe with fossil fuels emissions (excl. biofuels) relevant for the EF calculation
fossil_fuels_SIEC = annual_data_SIEC.loc[annual_data_SIEC['Fuel_type'].str.fullmatch(fossil_fuels_str)]

# sum up all values for each country and year in the dataframe with relevant fossil fuels emissions (without biofuels)
fossil_emissions_sum = fossil_fuels_SIEC.groupby(['TIME', 'GEO'], as_index = False).sum().drop(columns = ['ti_ele','GEP','own_use_ele','ti_CHP','ti_ele_total' ]) # eurostat 2023 update: use 'TIME_PERIOD' instead of 'TIME'

# fill the column correposnding to the total fossil-based emissions in the dataframe with the total electricity data 
emissions_types = ['Mt_CO2e_d', 'Mt_CO2_d', 'Mt_CH4_d', 'Mt_N2O_d']
for et in emissions_types:
    annual_total[et] = [fossil_fuels_SIEC.loc[(fossil_fuels_SIEC['TIME'] == time) & (fossil_fuels_SIEC['GEO'] == geo), et].sum() for time, geo in zip(annual_total['TIME'], annual_total['GEO'])]

annual_data_SIEC = pd.concat([annual_total, annual_data_SIEC], axis = 0, ignore_index = True) 
annual_data_ENTSOE = pd.concat([annual_total.rename(columns = {'Fuel_type': 'Generation_type'}), annual_data_ENTSOE], axis = 0, ignore_index = True) 

# In[12]: export dataframes to csv. format

# annual combustion emissions in the SIEC format (aggregtaed by fuel types)
annual_data_SIEC.to_csv(os.path.join(output_directory_path,'annual_data_SIEC.csv'))

# annual combustion emissions aggergated according to the ENTSO-E generation types
annual_data_ENTSOE.to_csv(os.path.join(output_directory_path, 'annual_data_ENTSOE.csv'))



