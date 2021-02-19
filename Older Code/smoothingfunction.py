# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 12:20:48 2020

@author: lpb20
"""

def smooth_df(df,smooth_points = 50):

    columns_num = (len(df.columns)-1) * smooth_points
    index_num =len(df.index)
    start_en = df.columns
    start_array = np.zeros([index_num , columns_num])
    final_df = pd.DataFrame(start_array)
    final_df.index = df.index
    for j in range(110):
        final_en = np.empty([])
        final_vals = np.empty([])
        
        for i in range(len(df.columns)-1):
            x1 = start_en[i]
            x2 = start_en[i+1]
            y1 = df.iloc[j,i] 
            y2 = df.iloc[j,i+1]
            temp_en = np.linspace(x1,x2,smooth_points)
            b = np.log(y2/y1) / np.log(x2/x1)
            a = y1 / (x1**b)
            
            final_vals_i =  a * temp_en ** b    
            final_vals = np.append(final_vals,final_vals_i)
            final_en = np.append(final_en , temp_en)
            
        final_en = final_en[1:]
        final_vals = final_vals[1:]  
        final_df.iloc[j,:] = final_vals
    final_df.columns = final_en
    final_df.fillna(0,inplace = True)
    
    data_p = final_df
    return data_p
