# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:09:23 2022

@author: Daria Blizniukova
"""
# In[1]: import the modules
    
import os
import pandas as pd
from datetime import datetime
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
import statistics
from matplotlib.lines import Line2D

# In[2]: define the functions

def change_fontsize(size, xticklabel, xlabel, ylabel):
    ax.tick_params(axis = 'y', labelsize = int(size))
    ax.set_xticklabels(xticklabel, fontsize = int(size))
    ax.set_ylabel(str(ylabel), fontsize = int(size))
    ax.set_xlabel(str(xlabel), fontsize = int(size))

def calc_annual_grid_mix_detailed(generation_DE):
    # aggregate less releveant generation type groups
    other_fossil_agg = ['Fossil_Oil', 'Fossil_Oil_shale',  'Other', 'Waste',  'Fossil_Peat']
    other_renewable_agg = ['Other_renewable', 'Geothermal', 'Hydro_Pumped_Storage', 'Hydro_excl._HPS', 'Marine']
    coal_agg = ['Fossil_Brown_coal/Lignite', 'Fossil_Hard_coal']
    wind_agg = ['Wind_Onshore', 'Wind_Offshore']
    
    # create an empty dataframe for storing the results
    df_result = pd.DataFrame(columns = ['Year', 'Area', 'Biomass', 'Coal', 'Fossil_Gas', 'Other', 'Nuclear', 'Wind', 'Solar', 'Other_renewable'])

    for df in list(generation_DE):
        area = df['Area'].iloc[0]
        year = df['Datetime'].dt.year.unique()

        # refine the dataframe
        df = df.set_index(['Area', 'Datetime'])
        df = df.drop(columns = ['Pumping_consumption'])

        # caculate net generation in GWh per 
        biomass_gwh = df['Biomass'].sum()
        fossil_gas_gwh = df['Fossil_Gas'].sum()
        nuclear_gwh = df['Nuclear'].sum()
        solar_gwh = df['Solar'].sum()
        coal_gwh = np.array([df[g] for g in coal_agg]).sum()
        wind_gwh = np.array([df[g] for g in wind_agg]).sum()
        other_fossil_gwh = np.array([df[g] for g in other_fossil_agg]).sum()
        other_renewable_gwh = np.array([df[g] for g in other_renewable_agg]).sum()

        # calculate total net generation
        df['Total'] = df.sum(axis = 1)
        total = df['Total'].sum()
        
        # calculate the shares
        biomass = biomass_gwh / total
        coal = coal_gwh / total
        fossil_gas = fossil_gas_gwh / total
        nuclear = nuclear_gwh / total
        wind = wind_gwh / total
        solar = solar_gwh / total
        other_fossil = other_fossil_gwh / total
        other_renewable = other_renewable_gwh / total
        
        # put the values into the results dataframe
        index = len(df_result.index)
        df_result.loc[index] = {'Area': area, 'Year': year, 
                                'Biomass': biomass, 'Coal': coal, 'Fossil_Gas': fossil_gas, 
                                'Other': other_fossil, 'Nuclear':nuclear, 'Other_renewable': other_renewable,
                                'Wind': wind, 'Solar': solar}
    
    # adjust the formats
    df_result['Year'] = df_result['Year'].astype(int)
    
    del(df, biomass, nuclear, coal, fossil_gas, other_fossil, other_renewable, wind, solar, total, index)
    return(df_result)

def calculate_monthly_Scope_shares(area, years, monthly_mean):
    results_df = pd.DataFrame(columns = ['Year','Month', 'Area', 'Scope_2_[g_CO2e/kWh]', 'Scope_3_[g_CO2e/kWh]', 'Scope_2_share', 'Scope_3_share'])
    for y, year in enumerate(years):
        df_area_year = monthly_mean[(monthly_mean['Year'] == year) & (monthly_mean['Area'] == area)].reset_index(drop = True)
        months = df_area_year['Month'].unique()
    
        for m in months:
            df_area_year_month = df_area_year[df_area_year['Month'] == m]
            
            # derive the Scope 2 and 3 EF
            Scope_2 = df_area_year_month['EF_net'].iloc[0]
            
            EF_p = df_area_year_month['EF_p'].iloc[0]
            EF_l = df_area_year_month['EF_l'].iloc[0]
            EF_up = df_area_year_month['EF_up'].iloc[0]
            Scope_3 = df_area_year_month['Scope_3'].iloc[0]
            EF_total = df_area_year_month['Total'].iloc[0]
            
            # calculate the shares
            Scope_2_share = Scope_2 / EF_total
            Scope_3_share = Scope_3 / EF_total       
        
            # put the values in the dataframe
            index = len(results_df.index) 
            results_df.loc[index] = {'Year': year, 'Month': m, 'Area': area,
                                     'Scope_2_[g_CO2e/kWh]': Scope_2, 'Scope_3_[g_CO2e/kWh]': Scope_3,
                                     'Scope_2_share': Scope_2_share * 100, 'Scope_3_share': Scope_3_share * 100}
    del(df_area_year, df_area_year_month, Scope_2, Scope_3, EF_p, EF_l, EF_up, EF_total)
    return(results_df)

def calculate_annual_Scope_shares(areas, years, annual_mean):
    results_df = pd.DataFrame(columns = ['Year', 'Area', 'Scope_2_[g_CO2e/kWh]', 'Scope_3_[g_CO2e/kWh]', 'Scope_2_share', 'Scope_3_share'])
    for year in years:
        for area in areas:    
            df_area_year = annual_mean[(annual_mean['Year'] == year) & (annual_mean['Area'] == area)].reset_index(drop = True)
             
            # derive the Scope 2 and 3 EF
            Scope_2 = df_area_year['EF_net'].iloc[0]
            
            EF_p = df_area_year['EF_p'].iloc[0]
            EF_l = df_area_year['EF_l'].iloc[0]
            EF_up = df_area_year['EF_up'].iloc[0]
            Scope_3 = df_area_year['Scope_3'].iloc[0]
            EF_total = df_area_year['Total'].iloc[0]
            
            # calculate the shares
            Scope_2_share = Scope_2 / EF_total
            Scope_3_share = Scope_3 / EF_total       
        
            # put the values in the dataframe
            index = len(results_df.index) 
            results_df.loc[index] = {'Year': year, 'Area': area,
                                      'Scope_2_[g_CO2e/kWh]': Scope_2,
                                      'Scope_3_[g_CO2e/kWh]': Scope_3,
                                      'Scope_2_share': Scope_2_share * 100, 
                                      'Scope_3_share': Scope_3_share * 100}
    del(df_area_year, Scope_2, Scope_3, Scope_2_share, Scope_3_share, EF_p, EF_l, EF_up, EF_total)
    return(results_df)

def calc_Scope_3_share(data_dict, areas, years):
    df_result = pd.DataFrame() 
    for key in list(data_dict):
        for year in years:
            for area in areas:
                if area in key and str(year) in key:
                    # print(key)
                    df = data_dict[key].copy()      
                    df['EF_total_[g_CO2e/kWh]'] = df['EF_c_[g_CO2e/kWh]'] + df['EF_up_[g_CO2e/kWh]']
                    df['Scope_3_[g_CO2e/kWh]'] = df['EF_total_[g_CO2e/kWh]'] - df['EF_net_[g_CO2e/kWh]']
                    df['Scope_3_[%]'] = df['Scope_3_[g_CO2e/kWh]'] / df['EF_total_[g_CO2e/kWh]'] * 100
                    df_short = df[['Area', 'Datetime', 'Scope_3_[g_CO2e/kWh]', 'Scope_3_[%]']]
                    df_result = pd.concat([df_result, df_short])
    return(df_result)

def calc_st_vs_annual_Scope_3(annual_Scope_3_shares, st_Scope_3_shares, areas, years):
    results_df = pd.DataFrame(columns = ['Area', 'Year', 'annual_Scope_3_[g_CO2e/kWh]', 'st_Scope_3_[g_CO2e/kWh]', 'annual_Scope_3_share_[%]', 'st_Scope_3_share_[%]']) 
    for area in areas:
        annual_df_area = annual_Scope_3_shares[annual_Scope_3_shares['Area'] == area]
        st_df_area = st_Scope_3_shares[st_Scope_3_shares['Area'] == area]
        
        for year in years:
            annual_df_area_year = annual_df_area[annual_df_area['Year'] == year]
            st_df_area_year = st_df_area[st_df_area['Datetime'].dt.year == year]
    
            annual_value = annual_df_area_year['Scope_3_[g_CO2e/kWh]'].iloc[0]
            annual_share = annual_df_area_year['Scope_3_share'].iloc[0]
            
            st_value = st_df_area_year['Scope_3_[g_CO2e/kWh]'].mean()
            st_share = st_df_area_year['Scope_3_[%]'].mean()
            
            # put the values in the dataframe
            index = len(results_df.index) 
            results_df.loc[index] =   {'Area': area, 'Year': year,
                                      'annual_Scope_3_[g_CO2e/kWh]': annual_value,
                                      'st_Scope_3_[g_CO2e/kWh]': st_value,
                                      'annual_Scope_3_share_[%]': annual_share, 
                                      'st_Scope_3_share_[%]': st_share}
    return(results_df)

def calc_short_term_total(data_dict, areas, years):
    df_result = pd.DataFrame() 
    for key in list(data_dict):
        for year in years:
            for area in areas:
                if area in key and str(year) in key:
                    df = data_dict[key].copy()      
                    df_short = df[['Area', 'Datetime', 'Total']]
                    df_result = pd.concat([df_result, df_short])
    return(df_result)

def calc_annual_EF_by_zone(df_list, result_list, year, results_annual):
    for df in df_list:
        area = df.Area[0]
        # retrive the annual value from results_annual
        df = results_annual[(results_annual['Area'] == area) & (results_annual['Datetime'] == year)]
        df = df.reset_index(drop = True)
        EF_total = (df['EF_c_[g_CO2e/kWh]'] + df['EF_up_[g_CO2e/kWh]']).iloc[0]
        # put the value in the list
        result_list.append(EF_total)
        
def calc_monthly_EF_by_zone(data_df, result_list, results_monthly, year, month):
    for df in data_df:
        area = df.Area[0]
        # derive monthly values from results_monthly
        df_m = results_monthly[(results_monthly['Area'] == area) & (results_monthly['Year'] == year) & (results_monthly['Month'] == month)]
        EF_total = (df_m['EF_c_[g_CO2e/kWh]'] + df_m['EF_up_[g_CO2e/kWh]']).iloc[0]
        result_list.append(EF_total)
        
def extract_monthly_data(data_df, month):
    for df in data_df:
        df['Month'] = df.Datetime.dt.month
        df.drop(df[df["Month"] != month].index, inplace = True)
        df.reset_index(drop = True, inplace = True)
        df['Total'] = df['EF_up_[g_CO2e/kWh]'] + df['EF_c_[g_CO2e/kWh]']
        cols_to_keep = ['Area', 'Total']
        df.drop(columns = [col for col in df if col not in cols_to_keep], inplace = True)

def calculate_annual_EF_p_shares(areas, years, results_annual):
    results_df = pd.DataFrame(columns = ['Year', 'Area', 'EF_p_[g_CO2e/kWh]', 'EF_p_share_from_total', 'EF_p_share_from_Sc_3'])
    for year in years:
        for area in areas:    
            df_area_year = results_annual[(results_annual['Datetime'] == year) & (results_annual['Area'] == area)].reset_index(drop = True)
             
            # derive the Scope 2 and 3 EF
            Scope_2 = df_area_year['EF_net_[g_CO2e/kWh]'].iloc[0]
            EF_p = df_area_year['EF_p_[g_CO2e/kWh]'].iloc[0]
            EF_l = df_area_year['EF_l_[g_CO2e/kWh]'].iloc[0]
            EF_up = df_area_year['EF_up_[g_CO2e/kWh]'].iloc[0]
            Scope_3 = EF_p + EF_l + EF_up
            
            # calculate the shares
            EF_total = Scope_2 + Scope_3
            EF_p_share = EF_p / EF_total
            EF_p_share_Sc_3 = EF_p / Scope_3
        
            # put the values in the dataframe
            index = len(results_df.index) 
            results_df.loc[index] = {'Year': year, 'Area': area,
                                      'EF_p_[g_CO2e/kWh]': EF_p,
                                      'EF_p_share_from_total': EF_p_share * 100,
                                      'EF_p_share_from_Sc_3': EF_p_share_Sc_3 * 100}
            del(df_area_year, Scope_2, Scope_3, EF_p)
    return(results_df)

def calc_abs_EF_by_type(data_dict, areas, years):
    df_result = pd.DataFrame()#columns = ['Area', 'Year', 'EF_net_[g_CO2e/kWh]', 'EF_p_[g_CO2e/kWh]', 'EF_l_[g_CO2e/kWh]', 'EF_up_[g_CO2e/kWh]']) 
    for key in list(data_dict):
        for year in years:
            for area in areas:
                if area in key and str(year) in key:
                    df = data_dict[key].copy()                          
                    df_short = df[['Area', 'Datetime', 'EF_net_[g_CO2e/kWh]', 'EF_p_[g_CO2e/kWh]', 'EF_l_[g_CO2e/kWh]', 'EF_up_[g_CO2e/kWh]']]
                    df_result = pd.concat([df_result, df_short])     
    return(df_result)

def calc_mean_sd_EF_total(areas, years, dict_EF_15_min_CO2e, annual_mean):
    result_df = pd.DataFrame(columns = ['Area', 'Year', 'min_value', 'max_value', 'range', 'annual_mean', 'st_mean', 'sd'])
    
    keys = ['{}_{}_EF_15_min_CO2e'.format(area, year) for area in areas for year in years]
    for key in keys:
        df = dict_EF_15_min_CO2e[key]  
        area = df.Area.iloc[0]
        year = df.Datetime.dt.year.iloc[0]
            
        if key == 'TransnetBW_2017_EF_15_min_CO2e':
            # start = datetime(2017, 4, 10, 2, 30)
            # end = datetime(2017, 4, 10, 3, 45)
            dt1 = datetime(2017, 4, 10, 2, 30)
            dt2 = datetime(2017, 4, 10, 2, 45)
            dt3 = datetime(2017, 4, 10, 3, 00)
            dt4 = datetime(2017, 4, 10, 3, 30)
            dt5 = datetime(2017, 4, 10, 3, 45)
            dt6 = datetime(2017, 4, 10, 3, 15)
            
            dts = [dt1, dt2, dt3, dt4, dt5, dt6]
            ind_dt = []
            
            for dt in dts:
                ind = df[df['Datetime'] == dt].index.values.astype(int)[0]
                ind_dt.append(ind)
            
            df = df.drop(index = ind_dt)
        
        df = df.dropna()
        st_mean = statistics.mean(df['Total'])
        sd = df['Total'].std()
        min_value = min(df['Total'])
        max_value = max(df['Total'])
        range_val = max_value - min_value
        
        annual_EF_df = annual_mean[(annual_mean['Year'] == year) & (annual_mean['Area'] == area)]
        annual_mean_value = annual_EF_df['Total'].iloc[0]
        
        # put the values into the results dataframe
        index = len(result_df.index)
        result_df.loc[index] = {'Area': area, 'Year': year, 'annual_mean': annual_mean_value, 'min_value': min_value, 'max_value': max_value,
                                 'range': range_val, 'st_mean': st_mean, 'sd': sd}
    return(result_df)

def calc_mean_sd_Scope_3(areas, years, dict_EF_15_min_CO2e, annual_mean):
    result_df = pd.DataFrame(columns = ['Area', 'Year', 'min_value', 'max_value', 'range', 'annual_mean', 'st_mean', 'sd'])
    
    keys = ['{}_{}_EF_15_min_CO2e'.format(area, year) for area in areas for year in years] 
    for key in keys:
        df = dict_EF_15_min_CO2e[key].copy()
        df['Scope_3_[g_CO2e/kWh]'] = df['EF_p_[g_CO2e/kWh]'] + df['EF_l_[g_CO2e/kWh]'] + df['EF_up_[g_CO2e/kWh]']
        area = df.Area.iloc[0]
        year = df.Datetime.dt.year.iloc[0]
            
        if key == 'TransnetBW_2017_EF_15_min_CO2e':
            # start = datetime(2017, 4, 10, 2, 30)
            # end = datetime(2017, 4, 10, 3, 45)
            dt1 = datetime(2017, 4, 10, 2, 30)
            dt2 = datetime(2017, 4, 10, 2, 45)
            dt3 = datetime(2017, 4, 10, 3, 00)
            dt4 = datetime(2017, 4, 10, 3, 30)
            dt5 = datetime(2017, 4, 10, 3, 45)
            dt6 = datetime(2017, 4, 10, 3, 15)
            
            dts = [dt1, dt2, dt3, dt4, dt5, dt6]
            ind_dt = []
            
            for dt in dts:
                ind = df[df['Datetime'] == dt].index.values.astype(int)[0]
                ind_dt.append(ind)
            
            df = df.drop(index = ind_dt)
            
        df = df.dropna()
        
        st_mean = statistics.mean(df['Scope_3_[g_CO2e/kWh]'])
        sd = df['Scope_3_[g_CO2e/kWh]'].std()
        min_value = min(df['Scope_3_[g_CO2e/kWh]'])
        max_value = max(df['Scope_3_[g_CO2e/kWh]'])
        range_val = max_value - min_value
        
        annual_EF_df = annual_mean[(annual_mean['Year'] == year) & (annual_mean['Area'] == area)]
        annual_mean_value = annual_EF_df['Total'].iloc[0]
        
        # put the values into the results dataframe
        index = len(result_df.index)
        result_df.loc[index] = {'Area': area, 'Year': year, 'annual_mean': annual_mean_value, 'min_value': min_value, 'max_value': max_value,
                                 'range': range_val, 'st_mean': st_mean, 'sd': sd}
        
    del(dt1, dt2, dt3, dt4, dt5, dt6, dts, ind_dt)
    return(result_df)

def calc_annual_mean_sd_by_EF_type(areas, years, dict_EF_15_min_CO2e):
    result_df = pd.DataFrame(columns = ['Area', 'Year', 'mean_EF_p', 'sd_EF_p', 'mean_EF_up', 'sd_EF_up', 'mean_EF_l', 'sd_EF_l'])
    keys = ['{}_{}_EF_15_min_CO2e'.format(area, year) for area in areas for year in years]     
    for key in keys:
        df = dict_EF_15_min_CO2e[key].copy()
        area = df.Area.iloc[0]
        year = df.Datetime.dt.year.iloc[0]
        
        # remove incomplete /error rows from the raw data
        if key == 'TransnetBW_2017_EF_15_min_CO2e':
            # start = datetime(2017, 4, 10, 2, 30)
            # end = datetime(2017, 4, 10, 3, 45)
            dt1 = datetime(2017, 4, 10, 2, 30)
            dt2 = datetime(2017, 4, 10, 2, 45)
            dt3 = datetime(2017, 4, 10, 3, 00)
            dt4 = datetime(2017, 4, 10, 3, 30)
            dt5 = datetime(2017, 4, 10, 3, 45)
            dt6 = datetime(2017, 4, 10, 3, 15)
            
            dts = [dt1, dt2, dt3, dt4, dt5, dt6]
            ind_dt = []
            
            for dt in dts:
                ind = df[df['Datetime'] == dt].index.values.astype(int)[0]
                ind_dt.append(ind)
            
            df = df.drop(index = ind_dt)
        
        # drop empty rows (case of 50Hertz 2019)
        df = df.dropna()
        
        # derive the annual mean for EF_p, EF_l and EF_up
        mean_p = statistics.mean(df['EF_p_[g_CO2e/kWh]'])
        mean_up = statistics.mean(df['EF_up_[g_CO2e/kWh]'])
        mean_l = statistics.mean(df['EF_l_[g_CO2e/kWh]'])
        
        # derive the SD for short-term EF_p, EF_l and EF_up
        sd_p = df['EF_p_[g_CO2e/kWh]'].std()
        sd_up = df['EF_up_[g_CO2e/kWh]'].std()
        sd_l = df['EF_l_[g_CO2e/kWh]'].std()
        
        # put the values into the results dataframe
        index = len(result_df.index)
        result_df.loc[index] = {'Area': area, 'Year': year, 'mean_EF_p': mean_p, 'sd_EF_p': sd_p,
                                         'mean_EF_up': mean_up, 'sd_EF_up': sd_up, 'mean_EF_l': mean_l, 'sd_EF_l': sd_l}    
    # del(dt1, dt2, dt3, dt4, dt5, dt6, dts, ind_dt, mean_p, mean_up, mean_l, sd_up, sd_p, sd_l)
    return(result_df)
        
def calc_annual_mean_EF_by_type(data_dict, areas, years, EF_types):
    result_df = pd.DataFrame()
    
    for key in list(data_dict):
        for year in years:
            for area in areas:
                if area in key and str(year) in key:
                    for EF_type in EF_types:
                        df = data_dict[key].copy()    
                        df['Scope_3_[g_CO2e/kWh]'] = df['EF_p_[g_CO2e/kWh]'] + df['EF_l_[g_CO2e/kWh]'] + df['EF_up_[g_CO2e/kWh]']
               
                        EF_df = df[['Area', 'Datetime', f'{EF_type}_[g_CO2e/kWh]']]
                        EF_df.insert(2, 'EF_type', EF_type)
                        
                        # change column name
                        EF_df = EF_df.rename(columns = {f'{EF_type}_[g_CO2e/kWh]': 'value_[g_CO2e/kWh]'})

                        # put the values into the results dataframe
                        result_df = pd.concat([result_df, EF_df])
    return(result_df)

def extract_st_grid_mix_shares(area, year, month):
    # load raw data for the area and year
    key = 'gene_15_min_{}_{}.csv'.format(area, year)
    raw_data = pd.read_csv(os.path.join(output_directory_path, key), index_col = 0, parse_dates= ['Datetime'])
    # filter out the month
    if month == None:
        raw_data_m = raw_data
    else:
        raw_data_m = raw_data[raw_data['Datetime'].dt.month == month]
    # calculate short-term totel net generation
    raw_data_m['Total'] = raw_data_m.drop(columns = ['Datetime', 'Pumping_consumption', 'Area']).sum(axis = 1)
    # create a dataframe for the results
    results_df = raw_data_m[['Datetime', 'Total']]
    
    # define generation type groups 
    other_fossil_agg = ['Fossil_Oil', 'Fossil_Oil_shale',  'Other', 'Waste', 'Fossil_Peat']
    other_renewable_agg = ['Other_renewable', 'Geothermal', 'Hydro_Pumped_Storage', 'Hydro_excl._HPS', 'Marine']
    coal_agg = ['Fossil_Brown_coal/Lignite', 'Fossil_Hard_coal']
    wind_agg = ['Wind_Onshore', 'Wind_Offshore']
    
    # calculate short-term shares of each generation type from the grid mix
    total_share = results_df['Total'] / 100
    results_df['biomass_%'] = raw_data_m['Biomass'] / total_share
    results_df['coal_%'] = sum(raw_data_m[g] for g in coal_agg)/ total_share
    results_df['gas_%'] = raw_data_m['Fossil_Gas']/ total_share
    results_df['nuclear_%'] = raw_data_m['Nuclear'] / total_share
    results_df['other_fossil_%'] = sum(raw_data_m[g] for g in other_fossil_agg) / total_share
    results_df['other_RES_%'] = sum(raw_data_m[g] for g in other_renewable_agg) / total_share
    results_df['solar_%'] = raw_data_m['Solar'] / total_share
    results_df['wind_%'] = sum(raw_data_m[g] for g in wind_agg) / total_share
    
    # check the sum of the single shares
    # results_df['check'] = results_df.drop(columns = 'Total').sum(axis = 1)

    return(results_df)

def aggregate_st_grid_mix_shares(df):
    
    # define generation type groups 
    fossil_agg = ['coal_%', 'gas_%', 'other_fossil_%']
    RES_agg = ['other_RES_%', 'solar_%', 'wind_%']
    
    df['total_fossil_%'] = sum(df[g] for g in fossil_agg)
    df['total_RES_%'] = sum(df[g] for g in RES_agg)
    
    # check the sum of the single shares
    df['check'] = df['total_fossil_%'] + df['total_RES_%'] + df['nuclear_%'] + df['biomass_%']

def calculate_sd(dict_shares, areas, years):
    # create a new dataframe
    results_df = pd.DataFrame(columns = ['Area', 'Year', 'sd'])
    
    # derive the SD
    for area in areas:
        for year in years:        
            key = 'EF_15_min_CO2e_{}{}'.format(area, year)
            sd = dict_shares[key]['Total'].std()
            
            # fille out the results' dataframe
            results_df.loc[len(results_df.index)] = {'Area': area, 'Year': year, 'sd': sd}
    return(results_df)

def aggregate_annual_grid_mix_shares(df):
    
    # define generation type groups 
    fossil_agg = ['Coal', 'Fossil_Gas', 'Other']
    RES_agg = ['Wind', 'Solar', 'Other_renewable']
    
    df['Total_Fossil'] = sum(df[g] for g in fossil_agg)
    df['Total_RES'] = sum(df[g] for g in RES_agg)
    
    # check the sum of the single shares
    df['check'] = df['Total_Fossil'] + df['Total_RES'] + df['Nuclear'] + df['Biomass']
    

# In[3]: define the directories

input_directory_path = os.path.join('input')
output_directory_path = os.path.join('output')
plots_main_dir = os.path.join('plots main')
plots_SI_dir = os.path.join('plots SI')
plots_temp_dir = os.path.join('plots temp')
output_paper_dir = os.path.join('output paper')

# In[4]: load the EF results (raw data for plotting)

scaled_entsoe_year = pd.read_csv(os.path.join(output_directory_path,'scaled_entsoe_year.csv'), index_col = 0)

# annual resolution level
annual_average = pd.read_csv(os.path.join(output_directory_path,'annual_average_results.csv'), index_col = 0)
annual_mean = pd.read_csv(os.path.join(output_directory_path,'annual_mean_results.csv'), index_col = 0)
annual_grid_mix = pd.read_csv(os.path.join(output_directory_path,'annual_grid_mix.csv'), index_col = 0)

# monthly resolution level
monthly_average = pd.read_csv(os.path.join(output_directory_path,'monthly_average_results.csv'), index_col = 0)
monthly_mean = pd.read_csv(os.path.join(output_directory_path,'monthly_short-term_mean.csv'), index_col = 0)

monthly_grid_mix_DE_detailed = pd.read_csv(os.path.join(output_directory_path,'monthly_grid_mix_DE_detailed.csv'), index_col = 0)
monthly_grid_mix_2020 = pd.read_csv(os.path.join(output_directory_path,'monthly_grid_mix_2020_detailed.csv'), index_col = 0)

# In[5]: create a dictionary with the zone- and year-specific data

areas = ['Germany', '50Hertz', 'Amprion', 'TenneT GER', 'TransnetBW']
years = [2017, 2018, 2019, 2020]

keys = [str(area) + "_" + str(year) + "_" + 'EF_15_min_CO2e' for area in areas for year in years]

dict_EF_15_min_CO2e = {}
for key in keys:
    dict_EF_15_min_CO2e[key] = pd.read_csv(os.path.join(output_directory_path,'{}.csv'.format(key)), index_col = 0)
    dict_EF_15_min_CO2e[key]['Datetime'] = pd.to_datetime(dict_EF_15_min_CO2e[key]['Datetime'], errors='coerce')
    
    # remove two rows with the missing data from the TenneT 2020 dataframe
    if key == 'TenneT GER_2020_EF_15_min_CO2e':
        dt1 = datetime(2020, 10, 25, 2, 45)
        dt2 = datetime(2020, 3, 29, 1, 45)
        
        ind1 = dict_EF_15_min_CO2e[key][dict_EF_15_min_CO2e[key]['Datetime'] == dt1].index.values.astype(int)[0]
        ind2 = dict_EF_15_min_CO2e[key][dict_EF_15_min_CO2e[key]['Datetime'] == dt2].index.values.astype(int)[0]
        ind = [ind1, ind2]
        
        dict_EF_15_min_CO2e[key] = dict_EF_15_min_CO2e[key].drop(index = [ind1, ind2])
        
# In[6]: calculate shares and absolute Scope 3 values for all TSO zones and years (annual average and short-term mean) for Tables S5 and S6

years = [2017, 2018, 2019, 2020]
areas = ['Amprion','50Hertz','Germany', 'TenneT GER', 'TransnetBW']

# calculate annual average Scope 3 emissions
annual_Scope_3_shares = calculate_annual_Scope_shares(areas, years, annual_mean)

#calculate short-term Scope 3 emissions
st_Scope_3_shares = calc_Scope_3_share(dict_EF_15_min_CO2e, areas, years)

# combine both dataframes
st_vs_annual_Scope_3 = calc_st_vs_annual_Scope_3(annual_Scope_3_shares, st_Scope_3_shares, areas, years)

annual_EF_p_shares = calculate_annual_EF_p_shares(areas, years, annual_average)

# In[7]: calculate short-term mean and SD values for EF_total, Scope_3 and single EF types by zone and year (Table S5 and S6)

# define the parameters
areas = ['Germany', '50Hertz', 'Amprion', 'TenneT GER', 'TransnetBW']
years = [2017, 2018, 2019, 2020]

# calculate the values for EF_total
total_mean_sd = calc_mean_sd_EF_total(areas, years, dict_EF_15_min_CO2e, annual_mean)

# extract the results as .csv file
total_mean_sd.to_csv(os.path.join(output_paper_dir,'short-term_mean_sd.csv'))

# calculate the values for Scope 3
scope_3_st_mean_sd = calc_mean_sd_Scope_3(areas, years, dict_EF_15_min_CO2e, annual_mean)

# caclulate the values for EF_p, EF_up and EF_l (for plotting)
scope_3_EF_types_mean_sd = calc_annual_mean_sd_by_EF_type(areas, years, dict_EF_15_min_CO2e)
        
# In[8]: define area and year specific colors and hatches
    
hatches_area = {'Germany': '',
                '50Hertz': '/',
                'Amprion': '///',
                'TenneT GER': '..' ,
                'TransnetBW': '\\'}  

hatches_year = {2017: '',
                2018: '///',
                2019: '..',
                2020: '\\'}

# define a color for each year
colors_year = {2017: '#4c72b0',
               2018: '#dd8452',
               2019: '#55a868',
               2020: '#c44e52'}

# for plotting grid mixes
gen_types = ['bio', 'o', 'c', 'fg', 'nuc', 'o_RES', 'w', 'sol']

col_names = {'bio': "Biomass",
              'o': "Other",
              'c': 'Coal',
              'fg': 'Fossil_Gas',
              'nuc': 'Nuclear',
              'o_RES': 'Other_renewable',
              'w': 'Wind',
              'sol': 'Solar'}

colors_g = {'bio': '#8ca252',
            'c':'#5f657f',
            'fg': '#8091a0',
            'nuc': '#E8C07D', #'#F0B27A'
            'o': '#9fbcbf',
            'w': '#89bedc',
            'sol': '#F1C40F',
            'o_RES': '#45B39D'}

labels_g = {'bio': 'biomass',
            'o': 'other fossil',
            'c':'coal',
            'fg': 'fossil gas',
            'nuc': 'nuclear', #'#F0B27A'
            'o_RES': 'other renewable',
            'w': 'wind',
            'sol': 'solar'}

colors_EF = {'EF_p': '#c44e52',
             'EF_l': '#55a868',
             'EF_up': '#8172b3'}

sns.set_theme(style="whitegrid")

# In[9]: Figure 4: Annual average consumption EF for Germany and its TSO zones from 2017 to 2020 divided into net generation (EF_net),
# pumping consumption (∆EF_P), T&D losses (EF_l), and upstream (EF_( up)). Single EF are differentiated by colors and the years by hatches. 
# The annotated values correspond to EF_Total.
    
# define parameters for plotting
areas = scaled_entsoe_year['GEO'].unique()
years = scaled_entsoe_year['TIME'].unique()
width = 0.3 # the width of the bars
xtra_space = 0.05
gap = 0.5
ind = (np.arange(len(areas)) + .15) * 2 - 0.8 # the x locations for the groups

# initialize the matplotlib figure
f, ax = plt.subplots(figsize=(10, 5))

for a, area in enumerate(areas):  
    for y, year in enumerate(years):
        df = annual_average
        annual_average_area_year = df[(df['Area'] == area) & (df['Datetime'] == int(year))]
        
        # derive data for plotting
        EF_net = annual_average_area_year['EF_net_[g_CO2e/kWh]'].iloc[0]
        EF_p = annual_average_area_year['EF_p_[g_CO2e/kWh]'].iloc[0]
        EF_l = annual_average_area_year['EF_l_[g_CO2e/kWh]'].iloc[0]
        EF_up = annual_average_area_year['EF_up_[g_CO2e/kWh]'].iloc[0]
        
        # plot stacked bars  
        x = y * (width + xtra_space) + a * (len(years) * (width + xtra_space) + gap)
        
        EF_net_bar = ax.bar(x, EF_net, width, bottom = 0, label = 'EF_net (Scope 2)', color = '#4c72b0', hatch = hatches_year[year])
        EF_p_bar = ax.bar(x, EF_p, width, bottom = EF_net, label = 'EF_p (Scope 3)', color = '#c44e52', hatch = hatches_year[year])
        EF_l_bar = ax.bar(x, EF_l, width, bottom = EF_p + EF_net, label = 'EF_l (Scope 3)', color = '#55a868', hatch = hatches_year[year])
        EF_up_bar = ax.bar(x, EF_up, width, bottom = EF_p + EF_net + EF_l, label = 'EF_up (Scope 3)', color = '#8172b3', hatch = hatches_year[year])
      
        ax.bar_label(ax.containers[-1], fmt = '%.0f', padding = 1.5, fontsize = 11)# 'medium') #11
plt.show()
    
# shift x tick labels
step = len(years) * (width + xtra_space) + gap
loc = np.array([a * step for a in range(0,len(areas))])
ax.set_xticks(loc + len(years) * (width + xtra_space) / 2.5)

# adjust grid lines
ax.xaxis.grid(False)
# sns.despine(left=True, bottom=True)

# change fontsize
change_fontsize(12, ['Germany', '50Hertz', 'Amprion', 'TenneT', 'TransnetBW'], '', r'g CO$_2$e/kWh') #12
    
# add a legend for EF types
l_EF_net = mpatches.Patch(facecolor = '#4c72b0', label = r'EF$_{net}$ (Scope 2)') 
l_EF_p = mpatches.Patch(facecolor = '#c44e52', label = r'ΔEF$_{p}$ (Scope 3)') #c44e52
l_EF_l = mpatches.Patch(facecolor = '#55a868', label = r'ΔEF$_{l}$ (Scope 3)') #55a868
l_EF_up = mpatches.Patch(facecolor = '#8172b3', label = r'EF$_{up}$ (Scope 3)') #8172b3

# dummy legend for distributing labels between two legend columns
l_dummy = mpatches.Patch(facecolor = 'w', label = '') 

# legend for the years
l_year = []
for year in years:
    l = mpatches.Patch(edgecolor = 'w', facecolor = 'grey', label = year, hatch = hatches_year[year])
    l_year.append(l)

# define parameters for placing th legend
kw = dict(bbox_to_anchor=(0.5, -0.3), columnspacing = 1.3, fontsize = 12)

# add the upper legend with EF types
leg1 = plt.legend(handles = [l_EF_net, l_EF_p, l_EF_l, l_EF_up], ncol=4, frameon=False, loc='lower center', **kw)
plt.gca().add_artist(leg1)

# add the lower legend with the year labels
leg2 = plt.legend(handles = l_year, ncol=4, frameon=False, loc='lower center', **kw)
leg2.remove()
leg1._legend_box._children.append(leg2._legend_handle_box)
leg1._legend_box.stale = True
   
plt.subplots_adjust(bottom=0.22)

plt.savefig('{}/Figure_4_annual_EF.png'.format(plots_main_dir), dpi = 300)

# In[10]: Figure 5: Annual grid mix composition based on net generation in Germany and its TSO zones between 2017 and 2020.

sns.set_style('whitegrid')

# define parameters for plotting
areas = scaled_entsoe_year['GEO'].unique()
width = 0.3 # the width of the bars
xtra_space = 0.05
gap = 0.5

# start the plot
df = annual_grid_mix
fig, ax = plt.subplots(figsize=(10,5))

for a, area in enumerate(areas):
    for y, year in enumerate(years):
        df_area_year = df[(df['Area'] == area) & (df['Year'] == int(year))]
                
        # define the location of a bar^
        x = y * (width + xtra_space) + a * (len(years) * (width + xtra_space) + gap)
        
        for g in gen_types:
            col = col_names[g]
            # value correponding to the generation type g
            value = df_area_year[col].iloc[0]
            # index of the generation type in the list of all generation types
            index = gen_types.index(g)
            # list of generation types prior to g in the list of all generation types
            index_list = [i for i in range(0, index)]
            g_list = [gen_types[i] for i in index_list]
            # sum of all values for generation types prior to g (= the hight of the new bar wto build upon the previous gen. types)
            bottom_val = sum(df_area_year[col_names[i]].iloc[0] for i in g_list)
            ax.bar(x, value, width, bottom = bottom_val, label = labels_g[g], color = colors_g[g], hatch = hatches_year[year])
plt.show()

# shift x tick labels
step = len(years) * (width + xtra_space) + gap
loc = np.array([a * step for a in range(0,len(areas))])
ax.set_xticks(loc + len(years) * (width + xtra_space) / 2.5)

# adjust grid lines
ax.xaxis.grid(False)

# change fontsize
change_fontsize(12, ['Germany', '50Hertz', 'Amprion', 'TenneT', 'TransnetBW'], '', '')# r'g CO$_2$e/kWh')
    
# convert y axis as percents axis
vals = [i for i in np.arange(0, 1.1, 0.1)]
plt.yticks(np.arange(0, 1.1, 0.1))
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
  
# define labels for the legend by generation type
g_legend_order = ['o', 'bio', 'fg', 'c', 'o_RES', 'nuc', 'sol', 'w']
l_gen = []
for g in g_legend_order:
    l = mpatches.Patch(edgecolor = 'w', facecolor = colors_g[g], label = labels_g[g])
    l_gen.append(l)
    
# define labels for the legend by area
l_year = []
for year in years:
    l = mpatches.Patch(edgecolor = 'w', facecolor = 'grey', label = year, hatch = hatches_year[year])
    l_year.append(l)

# arrange the legend location
kw = dict(bbox_to_anchor=(0.5, -0.36), fontsize = 12, columnspacing = 1.3)

leg1 = plt.legend(handles = l_gen, ncol = 4, frameon=False, loc='lower center', **kw)
plt.gca().add_artist(leg1)
# add the lower legend with the zone labels
leg2 = plt.legend(handles = l_year, ncol=4, frameon=False, loc='lower center', **kw)
leg2.remove()
leg1._legend_box._children.append(leg2._legend_handle_box)
leg1._legend_box.stale = True

plt.subplots_adjust(bottom=0.25)
plt.savefig('{}/Figure_5_annual_grid_mix.png'.format(plots_main_dir), bbox_inches = 'tight', dpi = 300)

# In[11]: Figure S1: Monthly average EF_Total including direct and indirect emissions in Germany from 2017 to 2020. 

# isolate data for DE only
monthly_mean_DE = monthly_mean[monthly_mean['Area'] == 'Germany']
# monthly_mean_DE['Total'] = monthly_mean_DE['EF_c_[g_CO2e/kWh]'] + monthly_mean_DE['EF_up_[g_CO2e/kWh]']

# plot monthly average total EF (upstream + consumption emissions)
g = sns.catplot(x= "Month", y= "Total", hue = "Year", kind="point", data = monthly_mean_DE, markers = '.', legend = False, frameon = True, height = 5, aspect = 1)#legend_out = False)
g.set_ylabels(r'g CO$_2$e/kWh', fontsize = 12)
g.set(xlabel = None)
g.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
sns.despine(right=False, top = False)

# derive total annual average (upstream + consumption emissions). Equivlent to EF_up_bar
annual_mean_DE = annual_mean[annual_mean['Area'] == 'Germany']

# add horizontal lines correponsind to the annual mean
for year in years:
    g.axes[0][0].axhline(y = annual_mean_DE[annual_mean_DE['Year'] == year]['Total'].iloc[0], color = colors_year[year], ls = '--')# xmin = xmin, xmax = xmax,
 
# add annotation with the annual average EF
l_year = []
for year in years:
    l = Line2D([0], [0], color = colors_year[year], label = str(year), lw = 2)
    l_year.append(l)

l_annual = Line2D([0], [0], color = 'gray', label = 'annual mean', linestyle = '--', lw = 1.5)

# arrange the legend location
kw = dict(bbox_to_anchor=(0.5, -0.25), fontsize = 11, columnspacing = 0.9)
leg1 = plt.legend(handles = l_year, ncol=4, frameon=False, loc='lower center', **kw)
plt.gca().add_artist(leg1)
leg2 = plt.legend(handles = [l_annual], ncol=1, frameon=False, loc='lower center', **kw)
leg2.remove()
leg1._legend_box._children.append(leg2._legend_handle_box)
leg1._legend_box.stale = True

plt.subplots_adjust(bottom=0.2)
plt.savefig('{}/Figure_S1_monthly_EF_total.png'.format(plots_SI_dir), bbox_inches = 'tight', dpi = 300)

# In[12]: Figure S2: Monthly average share of Scope 3 emissions from EF_Total in Germany from 2017 to 2020. 

years = [2017, 2018, 2019, 2020]
area = 'Germany'

# calculate the Scopes shares
area_Scope_shares_monthly = calculate_monthly_Scope_shares(area, years, monthly_mean)
area_Scope_shares_annual = calculate_annual_Scope_shares([area], years, annual_mean)

# start the plot
g = sns.catplot(x="Month", y="Scope_3_share", hue ='Year', kind="point", data = area_Scope_shares_monthly, markers = '.', legend = False, frameon = True, height = 5, aspect = 1)#legend_out = False)
g.set(xlabel = None, ylabel = None)
g.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])

fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
yticks = mtick.FormatStrFormatter(fmt)
for ax in g.axes.flat:
    ax.yaxis.set_major_formatter(yticks)

sns.despine(right=False, top = False)

# add horizontal lines correponsind to the annual average
for year in area_Scope_shares_annual['Year'].unique():
    g.axes[0][0].axhline(y = area_Scope_shares_annual[area_Scope_shares_annual['Year'] == year]['Scope_3_share'].iloc[0], color = colors_year[year], ls = '--')# xmin = xmin, xmax = xmax,

# add annotation for each year
l_year = []
for year in years:
    l = Line2D([0], [0], color = colors_year[year], label = year, lw = 2)
    l_year.append(l)

# add annotation for the annual average EF
l_annual = Line2D([0], [0], color = 'gray', label = 'annual mean', linestyle = '--', lw = 1.5)

# arrange the legend location
kw = dict(bbox_to_anchor=(0.5, -0.25), fontsize = 11, columnspacing = 0.9)
leg1 = plt.legend(handles = l_year, ncol=4, frameon=False, loc='lower center', **kw)
plt.gca().add_artist(leg1)
leg2 = plt.legend(handles = [l_annual], ncol=1, frameon=False, loc='lower center', **kw)
leg2.remove()
leg1._legend_box._children.append(leg2._legend_handle_box)
leg1._legend_box.stale = True

plt.subplots_adjust(bottom=0.2)
plt.savefig('{}/Figure_S2_monthy_mean_Scope_3_shares_{}.png'.format(plots_SI_dir, area), bbox_inches = 'tight', dpi = 300)

# In[13]: Figure 6: Development of the short-term EF_Total between 2017 and 2020 in Germany and its TSO zones. 
# The data for the TransnetBW zone in 2017 is excluded due to data gaps on a short-term level.

# derive the data
areas = ['Germany', '50Hertz', 'Amprion', 'TenneT GER', 'TransnetBW']
years = [2017, 2018, 2019, 2020]
data_boxplot = calc_short_term_total(dict_EF_15_min_CO2e, areas, years)

# delete TransnetBW 2017 values
data_boxplot = data_boxplot[(data_boxplot['Area'] != 'TransnetBW') | (data_boxplot['Datetime'].dt.year != 2017)]
data_boxplot.reset_index(drop = True)

sns.set_palette('deep')

fig, ax = plt.subplots(figsize=(9, 5))
df = data_boxplot[(data_boxplot['Area'] != 'TransnetBW') | (data_boxplot['Datetime'].dt.year != 2017)]
year = df.Datetime.dt.year
ax = sns.boxplot(y = 'Area', x = 'Total', data = df, hue = year, showfliers = True, saturation = 0.7, flierprops = {'marker': '|'}, showmeans = True, meanprops = {'marker':'o', 'markerfacecolor' : 'w','markeredgecolor':'black', 'markersize': 5})

# adjust the labels
ax.set_ylabel('')
ax.set_xlabel(r'g CO$_2$e/kWh', fontsize = 11)
ax.set_yticklabels(['Germany', '50Hertz', 'Amprion', 'TenneT', 'TransnetBW'])

# add annotation for each year
l_year = []
for year in years:
    l = mpatches.Patch(edgecolor = 'w', facecolor = colors_year[year], label = year)
    l_year.append(l)

l_mean, = [Line2D([0], [0], marker = 'o',color = 'white', markeredgecolor = 'black', label = 'short-term mean', linestyle = '')]
l_outlier, = [Line2D([0], [0], marker = '|',color = 'white', markeredgecolor = 'black', label = 'outlier', linestyle = '')]

# arrange the legend 
plt.legend([],[])
kw = dict(bbox_to_anchor=(0.5, -0.35), fontsize = 11, columnspacing = 1.3)
# add the legend for the year colors
leg1 = plt.legend(handles = l_year, ncol=4, frameon=False, loc='lower center', **kw)
plt.gca().add_artist(leg1)

# add the legend for the mean value
leg2 = plt.legend(handles = [l_mean, l_outlier], ncol = 2, frameon=False, loc='lower center', **kw)
leg2.remove()
leg1._legend_box._children.append(leg2._legend_handle_box)
leg1._legend_box.stale = True

plt.subplots_adjust(bottom=0.22)

plt.savefig('{}/Figure_6_short_term_EF_total.png'.format(plots_main_dir), bbox_inches = 'tight', dpi = 300)

# In[14]: Figure 7: Distribution of short-term Scope 3 emissions per kWh consumed electricity comprising pumping consumption (∆EF_P), T&D losses (∆EL_l), and upstream emissions (EF_up) in German and its TSO zones in 2020.
# does not include outliers
    
# derive the data
areas = ['Germany', '50Hertz', 'Amprion', 'TenneT GER', 'TransnetBW']
years = 2020
data_boxplot = calc_Scope_3_share(dict_EF_15_min_CO2e, areas, [years])
sns.set_style('whitegrid')

fig, ax = plt.subplots(figsize=(5, 4))
year = data_boxplot.Datetime.dt.year
area = data_boxplot['Area']
ax = sns.boxplot(x = area, y = 'Scope_3_[g_CO2e/kWh]', data = data_boxplot, color = '#d1d1d1', showfliers = False, flierprops = {'marker': 'x'}, showmeans = False)#, height = 5, aspect = 1)#, meanprops = {'marker':'o','markerfacecolor':'white','markeredgecolor':'black'})

# adjust the labels
ax.set_xlabel('')
ax.set_ylabel(r'g CO$_2$e/kWh', fontsize = 13)
plt.rc('ytick', labelsize=12)
ax.set_xticklabels(['Germany', '50Hertz', 'Amprion', 'TenneT', 'TransnetBW'], fontsize = 12)

plt.savefig('{}/Figure_7_short-term_Scope_3_distirbution_{}.png'.format(plots_main_dir, years), bbox_inches = 'tight', dpi = 300)

# In[15]: Figure 8: Distribution of short-term Scope 3 emissions per kWh consumed electricity attributed to pumping consumption (a),
# T&D losses (b), and upstream emissions (c) in German and its TSO zones in 2020.

# define the parameters
areas = ['Germany', '50Hertz', 'TenneT', 'Amprion', 'TransnetBW']
year = 2020
EF_types = ['EF_p', 'EF_l', 'EF_up'] #'Scope_3',

data_boxplot = calc_annual_mean_EF_by_type(dict_EF_15_min_CO2e, areas, [year], EF_types)

# define color per box (hue and hue_kws = {'color': ['#c44e52', '#55a868', '#8172b3']} does not work in the current matplotlib version)
# pal = {'color': ['#c44e52', '#55a868', '#8172b3']}

g = sns.FacetGrid(data_boxplot, col = 'EF_type', height = 3, aspect = 0.9) #, hue = EF_type, palette = colors_EF)
g.map(sns.boxplot, 'Area', 'value_[g_CO2e/kWh]', showfliers = False, showmeans = False, color = '#d1d1d1', width = 0.7) #meanprops = {'marker':'o', 'markerfacecolor' : 'w','markeredgecolor':'black', 'markersize': 5})#, palette = colors_EF)

# adjust axis labels
axes = g.axes.flatten() 
for ax in axes:
    ax.set_xlabel('')
    ax.set_ylabel(r'g CO$_2$e/kWh', fontsize = 12)
    ax.set_xticklabels(['Germany', '50Hertz', 'Amprion', 'TenneT', 'TransnetBW']) # rename TenneT_GER to TenneT
    ax.tick_params(axis = 'x', rotation = 25)
    ax.set_ylim(bottom = 0)

# adjust subplot titles 
# axes[0].set_title('Scope 3')
axes[0].set_title(r'ΔEF$_{p}$')
axes[1].set_title(r'ΔEF$_{l}$')
axes[2].set_title(r'EF$_{up}$')

# add letter labels to each of the subplots
axes[0].text(0, 54, 'a', bbox=dict(facecolor='none', edgecolor='black', pad=5.0))
axes[1].text(0, 54, 'b', bbox=dict(facecolor='none', edgecolor='black', pad=5.0))
axes[2].text(0, 54, 'c', bbox=dict(facecolor='none', edgecolor='black', pad=5.0))

plt.savefig('{}/Figure_8_short-term_Scope_3_EF_distribution_{}.png'.format(plots_main_dir, year), bbox_inches = 'tight', dpi = 300)

# In[16]: Figure 9: Mean of short-term electricity consumption EF divided into net generation (EF_net), pumping consumption (∆EF_P), 
# T&D losses (∆EF_l), and upstream emissions (EF_( up)) for Germany and its TSO zones in 2020, aggregated by season. 
# Winter season corresponds to the periods Jan – Mar and Oct – Dec, and summer season to Apr – Sep.

#  define parameters for plotting
year = 2020
labels = [r'EF$_{net}$ (Scope 2)', r'ΔEF$_{p}$ (Scope 3)', r'ΔEF$_{l}$ (Scope 3)', r'EF$_{up}$ (Scope 3)']
colors = ['#4c72b0', '#c44e52', '#55a868', '#8172b3']
seasons = ['winter', 'summer']

y_bottom_lim = {'winter': 80, 'summer': 50}
y_top_lim = {'winter': 550, 'summer': 650}

# load raw data
daily_mean_season = pd.read_csv(os.path.join(output_directory_path,'daily_mean_by_season.csv'), index_col = 0)
df = daily_mean_season[daily_mean_season['Year'] == year]
df['Area'] = df['Area'].replace('TenneT GER', 'TenneT')

areas = df.Area.unique()
times = np.unique(df['Time'].values)

# start the plot
sns.set_theme(style ='whitegrid')
fig, ax = plt.subplots(nrows = len(areas), ncols = len(seasons), sharey = 'row', figsize = (10,13))#

# begin iterating by row
for a, area in enumerate(areas):
    df_a = df[df['Area'] == area]
    ax[a,0].set_ylabel(r'g CO$_2$e/kWh')
    ax[a,0].set_title(f'{area}', loc = 'left') 
    
    # iterate by area within one row
    for s, season in enumerate(seasons):
        df_s = df_a[df_a['Season'] == season]
    
        # derive the EF
        EF_net = df_s['EF_net']
        EF_p = df_s['EF_p']
        EF_l = df_s['EF_l']
        EF_up = df_s['EF_up']
        EF_values = np.vstack([EF_net, EF_p, EF_l, EF_up])
    
        ax[a,s].stackplot(times, EF_values, labels = labels, colors = colors, alpha = 0.75)

        # adjust ticks and labels
        ticks = times[0:-1:4*6] # step: hours only (*4), each 6th element (*6) 
        ticks = np.append(ticks, '24:00:00')
        ticklabels = [f"{x}:00" for x in range(0,23,6)] + ['00:00']
        ticklabels[0] = '00:00' # otherwise '0:00'
        ax[a,s].set_xticks(ticks)
        ax[a,s].set_xticklabels(ticklabels)
        
        ax[a,s].set_xlim(left = 0, right = 96)
fig.tight_layout(h_pad = 1)

# adjust the y axes
ax[0,0].set_ylim(bottom = 220)
ax[1,0].set_ylim(bottom = 300)
ax[2,0].set_ylim(bottom = 80)
ax[3,0].set_ylim(bottom = 300)
ax[4,0].set_ylim(bottom = 50)

# add the season notation
for s, season in enumerate(seasons):
    ax[0,s].text(40, 470, f'{season}', fontsize = 13)

# add a legend for EF types
l_EF_net = mpatches.Patch(facecolor = '#4c72b0', label = r'EF$_{net}$ (Scope 2)') 
l_EF_p = mpatches.Patch(facecolor = '#c44e52', label = r'ΔEF$_{p}$ (Scope 3)') #c44e52
l_EF_l = mpatches.Patch(facecolor = '#55a868', label = r'ΔEF$_{l}$ (Scope 3)') #55a868
l_EF_up = mpatches.Patch(facecolor = '#8172b3', label = r'EF$_{up}$ (Scope 3)') #8172b3  
l_EF = [l_EF_net, l_EF_p, l_EF_l, l_EF_up]

kw = dict(bbox_to_anchor=(0, -0.45), columnspacing = 1.3, fontsize = 12)
plt.legend(handles = l_EF, ncol=4, frameon=False, loc='lower center', **kw)
plt.subplots_adjust(bottom = 0.10, hspace = 0.4)

plt.savefig('{}/Figure_9_short_term_EF_by_season_{}.png'.format(plots_main_dir, year), bbox_inches = 'tight', dpi = 300) 

