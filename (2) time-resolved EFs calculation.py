# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:28:47 2022

@author: Daria Blizniukova
"""

import os
import pandas as pd
from datetime import datetime
import numpy as np

# In[1]: define the directory

input_directory_path = os.path.join('input')
output_directory_path = os.path.join('output')

# In[2]: define functions, e.g. for uploading raw data (annual emissions and ENTSO-E raw data) and modifying dataframes

def load_entsoe(path, fn):
    df = pd.read_csv(os.path.join(path, fn), sep = ',', header = 0, parse_dates = True, ).replace("n/e", 0) #the 1st line contains header, the 1st row does not contain indexes
        
    # convert the datetime string to the datettime format
    df['Datetime'] = [datetime.strptime(x.split(" ")[0] + " " + x.split(" ")[1], '%d.%m.%Y %H:%M') for x in df.MTU]
    
    # insert the datetime string
    df.insert(2, 'Datetime', df.pop('Datetime'))
    del(df['MTU'])
    
    # rename the headings
    df.columns = df.columns.str.replace("  - Actual Aggregated", "").str.replace(" " , "_").str.replace("Hydro_Pumped_Storage__-_Actual_Consumption" , "Pumping_consumption")
    
    return df    

# In[3]: load time-resolved electricity generation data from ENTSO-E and convert the unit from MW to GWh

entsoe_15_min = load_entsoe(input_directory_path, 'entsoe_gene_DE_2019.csv') 

for col in entsoe_15_min.loc[:,'Biomass_[MW]':]:
    entsoe_15_min[col] = [float(i) * 0.25 / 1000 for i in entsoe_15_min.loc[:,"Biomass_[MW]":][col]]
entsoe_15_min.columns = entsoe_15_min.columns.str.replace('MW', 'GWh')

# In[4]: aggregate the ENTSO-E values to the annual level

entsoe_year = entsoe_15_min.groupby(['Area'], as_index = True).sum().T

# rename heading of a column with values
year = entsoe_15_min['Datetime'][0].year
col_name = entsoe_year.columns[0]
entsoe_year = entsoe_year.rename({col_name: year}, axis = 1)

#add a row with the annual total
entsoe_year.loc['Total_[GWh]'] = entsoe_year.sum() - entsoe_year.loc["Pumping_consumption_[GWh]"]

# In[5]: upload annual net generation values for scaling (case study: from Energiebilanzen AG (AGEB) updated in December 2021

AGEB_net_generation = pd.read_csv(os.path.join(input_directory_path, 'AGEB_net_gene_DE.csv'), sep = ';', float_precision = 'round_trip' )

# nonR waste is calculated as the sum of nonR household waste and industrial nonR waste
# the value for the generation type "other" is calculated by subtracting Hydro pumped storage (HPS) and nonR waste from "other" provided by AGEB
# pumping consumption referes to the total amount of electricity consumed by Hydro Pumped storage (including stored energy being realeased as generation from HPS)

# in AGEB the category "Fossil gas" includes both natural and coal-derived gas. Therefore, both categories are eggregated in the ENTSO-E data
entsoe_year.loc['Fossil_Gas_[GWh]'] += entsoe_year.loc['Fossil_Coal-derived_gas_[GWh]'] 
entsoe_year.drop(['Fossil_Coal-derived_gas_[GWh]'], inplace = True)

entsoe_15_min.loc[:, 'Fossil_Gas_[GWh]'] += entsoe_15_min.loc[:,'Fossil_Coal-derived_gas_[GWh]'] 
entsoe_15_min = entsoe_15_min.drop('Fossil_Coal-derived_gas_[GWh]', axis = 1)

# Hydro Run-of_river and poundage in AGEB is assumed to include the missing category "Water reservoir"
# such category is not mentioned in AGEB, but the annual value is significantly higher than the aggreagted ENTSO-E value
# --> Water reservoir and Run-of-river and pountage are aggregated to "Hydro" in the ENTSO-E raw data
# for further explanation of the hydroelectric plants classification see 
# ENTSOE-E (2019): "Hydropower modelling – New database complementing PECD", page 21
# https://docstore.entsoe.eu/Documents/SDC%20documents/MAF/2019/Hydropower_Modelling_New_database_and_methodology.pdf

entsoe_year.loc['Hydro_Run-of-river_and_poundage_[GWh]'] += entsoe_year.loc['Hydro_Water_Reservoir_[GWh]'] 
entsoe_year.drop(['Hydro_Water_Reservoir_[GWh]'], inplace = True)
entsoe_year = entsoe_year.rename({'Hydro_Run-of-river_and_poundage_[GWh]' : 'Hydro_excl._HPS_[GWh]'},axis = 0)

entsoe_15_min.loc[:,'Hydro_Run-of-river_and_poundage_[GWh]'] += entsoe_15_min.loc[:,'Hydro_Water_Reservoir_[GWh]'] 
entsoe_15_min = entsoe_15_min.drop(['Hydro_Water_Reservoir_[GWh]'], axis = 1)
entsoe_15_min = entsoe_15_min.rename({'Hydro_Run-of-river_and_poundage_[GWh]' : 'Hydro_excl._HPS_[GWh]'},axis = 1)

#reset index
entsoe_year = entsoe_year.reset_index()
entsoe_year = entsoe_year.rename({'index' : 'Generation_type'},axis = 1)

# In[6]: estimate scaling factors for each generation type following the approach by Agora Energiewende

# derive the year
year = str(entsoe_year.columns[1])

# extract column with the values correponsindg to the year together with generation types from the AGEB dataframe
AGEB_year = AGEB_net_generation[['Generation_type', year]].copy().fillna(0).set_index('Generation_type')
AGEB_year[year] = pd.to_numeric(AGEB_year[year])

# convert the units from TWh to GWh
AGEB_year = AGEB_year * 1000
AGEB_year.index = AGEB_year.index.str.replace('TWh', 'GWh')

# rename the hydropower category
AGEB_year = AGEB_year.rename({'Hydro_Run-of-river_and_poundage_[GWh]' : 'Hydro_excl._HPS_[GWh]'}, axis = 0)
AGEB_year = AGEB_year.drop(['Hydro_Water_Reservoir_[GWh]'], axis = 0)

# reset index
AGEB_year = AGEB_year.reset_index()

# merge the AGEB and the ENTSO-E aggeragted annual data per generation type
scaling_factors = pd.merge(entsoe_year, AGEB_year, on = ['Generation_type'], how = 'inner')
scaling_factors = scaling_factors.rename(columns = {scaling_factors.columns[1]:'ENTSOE', scaling_factors.columns[2] : 'AGEB' })

# calculate the scaling factors (sf)
scaling_factors['sf'] = scaling_factors['AGEB'] / scaling_factors['ENTSOE']

# In[7]: scale up time-resolved net generation values for each generation type

log = True
gen_type = entsoe_15_min.columns[3:] 
scaled_entsoe_15_min = entsoe_15_min.copy()
for g in gen_type:
    if log: print(f"{g}")
    factor = scaling_factors[scaling_factors.Generation_type == g].sf.iloc[0]
    
    if log: print(f"Scaling factor {factor}")
    scaled_entsoe_15_min.loc[:,g] = scaled_entsoe_15_min.loc[:,g] * factor 

# In[8]: check whether the scaling was done correctly

# aggergare the scaled values to one annual value per generation type. Remove the generation type "other renewable".
scaled_entsoe_year = scaled_entsoe_15_min.groupby(['Area'], as_index = True).sum().T 
scaled_entsoe_year = scaled_entsoe_year.drop("Other_renewable_[GWh]").rename({'Germany (DE)': 2019}, axis = 1)

# add a row with the annual total value
scaled_entsoe_year.loc['Total_[GWh]'] = scaled_entsoe_year.sum() - scaled_entsoe_year.loc["Pumping_consumption_[GWh]"]

# reset index
scaled_entsoe_year = scaled_entsoe_year.reset_index()
scaled_entsoe_year = scaled_entsoe_year.rename({'index' : 'Generation_type'},axis = 1)

log = True
for g in scaled_entsoe_year.Generation_type:
    expected_value = AGEB_year[AGEB_year.Generation_type == g].iloc[:, 1].iloc[0]# iloc[:, 1] = column with the 1st index ("year"); iloc[0] is the first and only element of the row
    calc_value = scaled_entsoe_year[scaled_entsoe_year.Generation_type == g].iloc[:, 1].iloc[0]
    difference = expected_value - calc_value
    if difference > 0.0001: #there is no specific reason behind this treshold
        if log: print(f"The expected value for {g} is: {expected_value}, the calculated value {calc_value}, difference is {difference} GWh")
    del (difference, expected_value, calc_value)

# There will be a differnce between the "Total" value calculated here and published by AGEB due to the rounding in the AGEB file. 
# The scaled sum (calculated here) should be more precise.

# In[9]: upload direct annual emissions claculated in the "(1)annual_emissions_calculation" code
    
# upload direct annual emissions calculated in the code "(1)annual_emissions_calculation"
annual_emissions = pd.read_csv(os.path.join(output_directory_path,'annual_data_ENTSOE.csv'))

# extract data on direct emissions for the country (Germany) and year from the "annual_emissions" dataframe
annual_emissions_year = annual_emissions.loc[(annual_emissions['GEO'] == 'Germany') & (annual_emissions['TIME'] == int(year))] # Warning: 'TIME' was renamed to 'TIME_PERIOD' in the omst recent eurostat update (January 2023)
annual_emissions_year = annual_emissions_year.drop(columns = ['Unnamed: 0', 'TIME', 'GEO', 'ti_ele', 'GEP', 'own_use_ele', 'ti_CHP', 'ti_ele_total'])
annual_emissions_year.reset_index(drop = True, inplace = True)

annual_emissions_year = annual_emissions_year.replace({'hard_coal' : 'Fossil_Hard_coal',
                                                       'brown_coal' : 'Fossil_Brown_coal/Lignite',
                                                       'fossil_gas' : 'Fossil_Gas',
                                                       'peat' : 'Fossil_Peat',
                                                       'oil_shale' : 'Fossil_Oil_shale',
                                                       'nonR_waste' : 'Waste',
                                                       'biomass' : 'Biomass',
                                                       'oil' : 'Fossil_Oil',
                                                       'total': 'Total'})

# In[10]: create a new dataframe for estimating annual net generation EF (EF_net,g) for calculating direct quarter-hourly emissions

annual_net_EF_year = scaled_entsoe_year.copy()

# rename the generation types for matching two dataframes
annual_net_EF_year['Generation_type'] = annual_net_EF_year['Generation_type'].str.replace(r'[][]', '', regex = True).str.replace('_GWh', '')
annual_net_EF_year.rename(columns = {annual_net_EF_year.columns[1]: 'net_gene_[GWh]'}, inplace = True)

# merge two dataframes based on the generation types relevant for direct emissions calculation
annual_net_EF_year = annual_net_EF_year.merge(annual_emissions_year, how = 'outer').set_index('Generation_type')

# remove the rows not relevant for calculating direct combustion emissions
rows_to_del = ['Total', 'Geothermal', 'Hydro_Pumped_Storage', 'Pumping_consumption', 'Hydro_excl._HPS', 'Marine', 'Nuclear', 'Solar', 'Wind_Offshore', 'Wind_Onshore']
annual_net_EF_year = annual_net_EF_year.drop(index = rows_to_del).fillna(0)
del(rows_to_del)

# In[11]: calculate annual net generation fmission factors (EF_net,g) based on the scaled annual net generation and combustion emissions
# Catrgory 'Other' consists of non specified fossil fuels. In this step, the correpoding EF is not yet calculated since the annual combustion emissions are assigned to the specific categories.

annual_net_EF_year['kg_CO2e/kWh_d'] = annual_net_EF_year['Mt_CO2e_d'] / annual_net_EF_year['net_gene_[GWh]'] * 10**3
annual_net_EF_year['kg_CO2/kWh_d'] = annual_net_EF_year['Mt_CO2_d'] / annual_net_EF_year['net_gene_[GWh]'] * 10**3
annual_net_EF_year['kg_CH4/kWh_d'] = annual_net_EF_year['Mt_CH4_d'] / annual_net_EF_year['net_gene_[GWh]'] * 10**3
annual_net_EF_year['kg_N2O/kWh_d'] = annual_net_EF_year['Mt_N2O_d'] / annual_net_EF_year['net_gene_[GWh]'] * 10**3

annual_net_EF_year = annual_net_EF_year.fillna(0)

# In[12]: calculate annual net generation EF for 'Other' as a weighted average of the EF of fossil fuels based on their net generation share

emission_factors = ['kg_CO2e/kWh_d', 'kg_CO2/kWh_d', 'kg_CH4/kWh_d', 'kg_N2O/kWh_d']
for e in emission_factors:
    
    # net generation by all fossil fuels including 'Other', but excl. 'Biomass'
    net_gene_biomass = annual_net_EF_year.loc['Biomass', 'net_gene_[GWh]']
    net_gene_fossil = sum(annual_net_EF_year['net_gene_[GWh]']) - net_gene_biomass

    # net generation by 'Other' only
    net_gene_other = annual_net_EF_year.loc['Other', 'net_gene_[GWh]']
    net_gene_fossilXother = net_gene_fossil - net_gene_other
    
    # the weighted average = sum(net generation * EF_net) / sum(net generation)
    productXbiomass = sum(annual_net_EF_year['net_gene_[GWh]'] * annual_net_EF_year[e]) - net_gene_biomass * annual_net_EF_year.loc['Biomass', e]
    annual_net_EF_year.loc['Other', e] = productXbiomass / net_gene_fossilXother    

del(e, net_gene_biomass, net_gene_fossil, net_gene_other, net_gene_fossilXother, productXbiomass)
    
# In[13]: upload upsteam emissions factor from Umweltbundesamt (UBA) and combine with the direct EF

annual_EF_CO2e_up = pd.read_csv(os.path.join(input_directory_path, 'upstream_EF.csv'), sep = ';')

# In[14]: calculate the upstream factor (EF_up) for 'Biomass' based on the GEP shares of from 'Primary Solid biofuels' and 'Biogases'.

# According to SIEC, Biomass in Germany = 'Primary Solid biofuels' + 'Biogases' + 'Other liquid biofuels'
# The share of 'Other liquid biofuels' is ignored due to its low share on the gross electricity production (GEP) from biomass (< 1 %) and the inabilty to determine the exact fuel type and finding an approproate upstream EF.
# The reason for using GEP instead of net electricity generation is the lack of data (neither AGEB not ENTSO-E differentiate between the subtypes of biomass).

# upload Energy Balances
nrg_bal = pd.read_csv(os.path.join(input_directory_path, 'nrg_bal_c_v_14_04_2022.csv'), sep = ',', header = 0)

# derive the annual GEP from Biogases
biogases_GEP = float(nrg_bal.loc[(nrg_bal['NRG_BAL'] =='Gross electricity production') 
                                & (nrg_bal['TIME'] == int(year)) 
                                & (nrg_bal['SIEC'] == 'Biogases') 
                                & nrg_bal['GEO'].str.contains('Germany')]['Value'].iloc[0].replace(" ", "")) #filter out the row and take the first value from the "Value" column

solid_biofuels_GEP = float(nrg_bal.loc[(nrg_bal['NRG_BAL'] =='Gross electricity production') 
                                & (nrg_bal['TIME'] == int(year)) 
                                & (nrg_bal['SIEC'] == 'Primary solid biofuels') 
                                & nrg_bal['GEO'].str.contains('Germany')]['Value'].iloc[0].replace(" ", ""))

# for checking the share of liquid biofuels 
# liquid_biofuels_GEP = float(nrg_bal.loc[(nrg_bal['NRG_BAL'] =='Gross electricity production') 
#                                 & (nrg_bal['TIME'] == int(year)) 
#                                 & (nrg_bal['SIEC'] == 'Other liquid biofuels') 
#                                 & nrg_bal['GEO'].str.contains('Germany')]['Value'].iloc[0].replace(" ", ""))

# calculate total biomass GEP
biomass_GEP = biogases_GEP + solid_biofuels_GEP # + liquid_biofuels_GEP

# calculate the corresponding GEP shares on the GEP of Biomass (liquid biomass is excluded)
biogases_share = biogases_GEP / biomass_GEP
solid_biofuels_share = solid_biofuels_GEP / biomass_GEP
# liquid_biofuels_share = liquid_biofuels_GEP/ biomass_GEP

# derive the upstream emission factors EF_up for biogas and solid biofuels
biogases_EF_up = annual_EF_CO2e_up.loc[annual_EF_CO2e_up['Generation_type'] == 'Biogases']['kg_CO2e/kWh_up'].iloc[0]
solid_biofuels_EF_up = annual_EF_CO2e_up.loc[annual_EF_CO2e_up['Generation_type'] == 'Primary solid biofuels']['kg_CO2e/kWh_up'].iloc[0]

# calculate the upstream EF for biomass as a weighted average of biogases and solid biofuels
biomass_EF_up = biogases_share * biogases_EF_up + solid_biofuels_share * solid_biofuels_EF_up

# add the aggregated EF to the dataframe as the last row
annual_EF_CO2e_up.loc[len(annual_EF_CO2e_up.index)] = {'Generation_type': 'Biomass',
                                                       'kg_CO2e/kWh_up': biomass_EF_up,
                                                       'Source' : 'calculated in the code'}

del(biomass_GEP, biogases_share, solid_biofuels_share, biogases_EF_up, solid_biofuels_EF_up)

# In[15]: calculate the upstream factor for 'Other' as weighted average of the EF_up fossil fuels based on the annual net generation

# create a dataframe with the annual net generation values and upstream EF for the fossil fuels 
fossils = pd.DataFrame(columns = ['Generation_type', 'net_generation', 'EF_up'])

# define the generation types relevant for the calculation
fossil_types = ['Fossil_Brown_coal/Lignite', 'Fossil_Gas', 'Fossil_Hard_coal', 'Fossil_Oil', 'Other', 'Waste']

for f in fossil_types:
    net_generation = scaled_entsoe_year.loc[scaled_entsoe_year['Generation_type'] == str(f)+'_[GWh]'].iloc[:, 1].iloc[0]
    try:
        EF_up = annual_EF_CO2e_up.loc[annual_EF_CO2e_up['Generation_type'] == f]['kg_CO2e/kWh_up'].iloc[0]
   
    # there is no EF_up for the category 'Other' --> put a 0 instead
    except IndexError:
        EF_up = 0

    # put the values as a new line at the end of the 'fossils' dataframe
    fossils.loc[len(fossils.index)] = [f, net_generation, EF_up]

del(fossil_types)

# calculate  weighted average of the upstream EF to be assigned as the EF_up for 'Other'
net_gene_fossil = sum(fossils['net_generation'])
net_gene_other = fossils.loc[fossils['Generation_type'] == 'Other']['net_generation'].iloc[0]
net_gene_fossilXother = net_gene_fossil - net_gene_other
other_EF_up = sum(fossils['net_generation'] * fossils['EF_up'])/net_gene_fossilXother

# add the value to the dataframe as the last row
annual_EF_CO2e_up.loc[len(annual_EF_CO2e_up.index)] = {'Generation_type': 'Other',
                                                       'kg_CO2e/kWh_up': other_EF_up,
                                                       'Source' : 'calculated in the code'}

# In[16]: merge two dataframes with the direct and upstream annual EF

annual_net_EF_year = annual_net_EF_year.reset_index()
annual_EF = pd.merge(annual_net_EF_year, annual_EF_CO2e_up, how = 'outer').fillna(0)

annual_EF = annual_EF.drop(columns = ['net_gene_[GWh]','Description', 'Source', 'Link'])

# In[17]: prepare data for the calculation of short-term emissions
    
#create a copy dataframe for short-term emissions calculation without units in the heading
template = scaled_entsoe_15_min.copy()
template.columns = template.columns.str.replace(r'[][]', '', regex = True).str.replace('_GWh' , '')

# In[18]: calculate short-term emissions for each time interval

# filter out generation types and create a temporary dataframe with intersecting columns only
common_gen_types = np.intersect1d(annual_EF.Generation_type, template.columns) 
generation = template[common_gen_types]

emission_types = ['Mt_CO2e_d', 'Mt_CO2_d', 'Mt_CH4_d', 'Mt_N2O_d', 'Mt_CO2e_up']
emission_factors = ['kg_CO2e/kWh_d', 'kg_CO2/kWh_d', 'kg_CH4/kWh_d', 'kg_N2O/kWh_d', 'kg_CO2e/kWh_up']

#create a dataframe for the emissions excl. biomass
emissions_15_minXbio = scaled_entsoe_15_min[['Area', 'Datetime']] 

#create a dataframe for Biomass emissions only for transparency
emissions_15_min_bio = scaled_entsoe_15_min[['Area', 'Datetime']]

log = True
for e, ef in zip(emission_types, emission_factors):
    emissions = generation.copy()
    for g in generation:      
        if log: print(f"{g}")
        
        # derive the annual EFs for the generation type g
        EF_value = annual_EF[annual_EF.Generation_type == g][ef].iloc[0] 
        if log: print("The EF is %.4f" % EF_value)     
        
        # rewrite the existing values
        emissions.loc[:,g] = generation.loc[:, g] * EF_value 
    emissions_15_minXbio[e] = emissions.sum(axis = 1) - emissions['Biomass']
    emissions_15_min_bio[e] = emissions['Biomass']

# In[19]: derive net electricity generation and pumping consumption for each time interval from the scaled ENTSO-E data

# add a column with the total net generation in each time interval
scaled_entsoe_15_min['Net_generation_[GWh]'] = scaled_entsoe_15_min.drop('Pumping_consumption_[GWh]', axis = 1).sum(axis = 1)

# create a new dataframe with time-resolved values for net generation and pumping consumption
ele_supply_15_min = scaled_entsoe_15_min[['Area', 'Datetime', 'Net_generation_[GWh]', 'Pumping_consumption_[GWh]']]

# In[20]: calculate distribution losses (DL) and final consumption for each time interval

# derive the annual distribution losses (DL) value from Energy Balances for the country and year
dl_year = nrg_bal.loc[(nrg_bal['NRG_BAL'] =='Distribution losses') 
                 & (nrg_bal['TIME'] == int(year)) 
                 & (nrg_bal['SIEC'] == 'Electricity') 
                 & nrg_bal['GEO'].str.contains('Germany')]['Value'].iloc[0] #filter out the row and take the first value from the "Value" column

dl_year = float(dl_year.replace(' ', '')) 

# calculate the share of DL relative to the annual net generation
dl_year_share = dl_year / scaled_entsoe_year.loc[scaled_entsoe_year['Generation_type'] == 'Total_[GWh]'][int(year)].iloc[0]

# add a column with the values
ele_supply_15_min['T&D_losses_[GWh]'] = dl_year_share * (ele_supply_15_min['Net_generation_[GWh]'] - ele_supply_15_min['Pumping_consumption_[GWh]'])

# calculate the final amount of electricity available (= consumption)
ele_supply_15_min['Consumption_[GWh]'] = ele_supply_15_min['Net_generation_[GWh]'] - ele_supply_15_min['Pumping_consumption_[GWh]'] - ele_supply_15_min['T&D_losses_[GWh]'] 

# In[21]: calculate time-resolved EFs for direct GHG emissions at the stage of net electricity generation

# merge time-resolved electricty data and emissions
EF_15_min = pd.merge(emissions_15_minXbio, ele_supply_15_min, on = ['Area', 'Datetime'])

# calculate EF_net for direct GHG emissions
emission_types = ['Mt_CO2e_d', 'Mt_CO2_d', 'Mt_CH4_d', 'Mt_N2O_d']
for e in emission_types:
    EF_15_min['EF_net_[kg_' + e[3:-2] + '/kWh]'] = EF_15_min[e] * 1000 / EF_15_min['Net_generation_[GWh]']

# In[22]: calculate delta EFs for the share of emissions allocated to the pumping consumption
 #formula: emissions / (net_generation - pumping_consumption) - EF_net
    
for e in emission_types:
    d = EF_15_min['Net_generation_[GWh]'] - EF_15_min['Pumping_consumption_[GWh]']
    EF_net = EF_15_min['EF_net_[kg_' + e[3:-2] + '/kWh]']
    EF_15_min['EF_p_[kg_' + e[3:-2] + '/kWh]'] = EF_15_min[e] * 1000 / d - EF_net
    del(d, EF_net)

# In[23]: calculate delta EFs for the share of emissions allocated to T&D
 #formula: emissions / (net_generation - pumping_consumption - T&D) - EF_net - EF_p

for e in emission_types:
    d = EF_15_min['Net_generation_[GWh]'] - EF_15_min['Pumping_consumption_[GWh]'] - EF_15_min['T&D_losses_[GWh]']
    EF_net = EF_15_min['EF_net_[kg_' + e[3:-2] + '/kWh]']
    EF_p = EF_15_min['EF_p_[kg_' + e[3:-2] + '/kWh]']
    EF_15_min['EF_l_[kg_' + e[3:-2] + '/kWh]'] = EF_15_min[e] * 1000 / d - EF_net - EF_p
    del(d, EF_net, EF_p)
 
# In[24]: calculate consumption EFs (EF_c) for direct and upstream GHG emissions
# formula: emissions / (net_generation - pumping_consumption - T&D)

for e in emission_types:
    d = EF_15_min['Net_generation_[GWh]'] - EF_15_min['Pumping_consumption_[GWh]'] - EF_15_min['T&D_losses_[GWh]']
    EF_15_min['EF_c_[kg_' + e[3:-2] + '/kWh]'] = EF_15_min[e] * 1000 / d
    del(d)

EF_15_min['EF_up_[kg_CO2e/kWh]'] = EF_15_min["Mt_CO2e_up"] * 1000 / EF_15_min['Consumption_[GWh]']

# In[25]: check whether EF_c = EF_net + EF_p + EF_l for CO2e

log = True
for row in EF_15_min:
    EF_net = EF_15_min['EF_net_[kg_CO2e/kWh]']
    EF_p = EF_15_min['EF_p_[kg_CO2e/kWh]']
    EF_l = EF_15_min['EF_l_[kg_CO2e/kWh]']
    EF_c = EF_15_min['EF_c_[kg_CO2e/kWh]']
    EF_15_min['EF_c_check_[kg_CO2e/kWh]'] = EF_net + EF_p + EF_l 
    EF_c_check = EF_15_min['EF_c_check_[kg_CO2e/kWh]']
    difference = EF_c_check.iloc[0] - EF_c.iloc[0]
    if difference > 0.0001:
        if log:
            print(f"The sum of EF is: {EF_c_check}, the aggregated EF_c is {EF_c}, the difference is {difference} [kg CO2e/kWh]")
    del(EF_net, EF_p, EF_l, EF_c, EF_c_check, row)

# In[26]: extract the EF in CO2e only

EF_15_min_CO2e = EF_15_min[['Area', 'Datetime', 'EF_c_[kg_CO2e/kWh]','EF_net_[kg_CO2e/kWh]', 'EF_p_[kg_CO2e/kWh]', 'EF_l_[kg_CO2e/kWh]', 'EF_up_[kg_CO2e/kWh]']]

country = EF_15_min_CO2e['Area'].iloc[0]
EF_15_min_CO2e.to_csv(os.path.join(output_directory_path, str(country + '_' + year + '_' + 'EF_15_min_CO2e.csv')))
 
















