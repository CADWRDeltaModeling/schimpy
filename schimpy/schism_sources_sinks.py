# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import yaml
from schimpy.schism_setup import create_schism_setup

def read_source_sink_yaml(yaml_fn):
    """ Read source sink yaml file
    
    Parameters
    ----------
    yaml_fn : str
        source_sink.yaml filename.

    Returns
    -------
    df_sources : PANDAS dataframe
        For sources
        
    df_sinks : PANDAS dataframe
        For sinks
    """
    with open(yaml_fn, 'r') as file:
        data = yaml.safe_load(file)
    if 'sources' in data.keys():
        df_sources = pd.DataFrame.from_dict(data['sources'],orient='index',
                                            columns=['x','y'])
    else:
        df_sources = pd.DataFrame() 
    if 'sinks' in data.keys():
        df_sinks = pd.DataFrame.from_dict(data['sinks'],orient='index',
                                          columns=['x','y'])
    else:
        df_sinks = pd.DataFrame()    
    return df_sources, df_sinks

def write_source_sink_yaml(df_sources,df_sinks,yaml_fn):
    """ Write source sink to yaml
    
    Parameters
    ----------
    df_sources : PANDAS dataframe
        with colums:['x','y']
    df_sinks : TYPE
        DESCRIPTION.
    yaml_fn : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    dict_sources = {}
    for r in df_sources.iterrows():
        dict_sources[r[0]] = [float(r[1].x),float(r[1].y)]    
    dict_file = {'sources':dict_sources}
    
    dict_sinks = {}
    for r in df_sinks.iterrows():
        dict_sinks[r[0]] = [float(r[1].x),float(r[1].y)]
    dict_file.update({'sinks':dict_sinks})
            
    with open(yaml_fn, 'w') as file:
        documents = yaml.dump(dict_file, file)    
    
def concat_source_sink_yaml(yaml_fn1,yaml_fn2,new_yaml):
    df_sources1,df_sinks1 = read_source_sink_yaml(yaml_fn1)
    df_sources2,df_sinks2 = read_source_sink_yaml(yaml_fn2)
    df_sources = pd.concat([df_sources1,df_sources2])
    df_sinks = pd.concat([df_sinks1,df_sinks2])
    write_source_sink_yaml(df_sources,df_sinks,new_yaml)

def read_source_sink_csv(csv_fn):
    th = pd.read_csv(csv_fn,sep=' ')
    th= th.set_index(pd.to_datetime(th.datetime))
    th = th.asfreq(pd.infer_freq(th.index))
    return th

def write_source_sink_csv(th, csv_fn):
    th.to_csv(csv_fn,sep=' ',index=False)

def unsorted_unique(a):
    """ Find unique elements with no sorting
    
    `np.unique` automatically sorts the elements. This function performs unique
    without sorting. 
    
    Parameters
    ----------
    a : List or array
        Input array

    Returns
    -------
    a_unique : array
        Unsorted unique array
    """

    indexes = np.unique(a, return_index = True)[1]
    a_unique = np.array(a)[sorted(indexes)]
    return a_unique

def concat_msource_csv(csv_fn1,csv_fn2,merged_source_sink_in,
                       csv_merged,freq='infer',how='left'):
    th1 = read_source_sink_csv(csv_fn1)
    th2 = read_source_sink_csv(csv_fn2)
    var1 = unsorted_unique([s.split('_')[0] for s in th1.columns[1:]])
    var2 = unsorted_unique([s.split('_')[0] for s in th2.columns[1:]])
    # check if there are variables other than T and S
    var1 = [v for v in var1 if v not in ['T','S']]
    var2 = [v for v in var2 if v not in ['T','S']]
    # if any(var1) or any(var2): check if var2 or var1 includes the other array. 
    # if not, I will need to specify the order of the modules. 
    if set(var1).issubset(var2):
        var = np.append(['T','S'],var2)
    elif set(var2).issubset(var1):
        var = np.append(['T','S'],var1)
    else:
        raise Exception("Multiple modules exist and the order of the modules is unclear")

    # merged_source_sink_in: the merged source_sink.in or source_sink.yaml file 
    # where the data sources are from csv_fn1, csv_fn2. 
    if merged_source_sink_in.endswith('yaml'):
        df_sources,df_sinks = read_source_sink_yaml(merged_source_sink_in)
    elif merged_source_sink_in.endswith('in'):
        df_sources,df_sinks = read_source_sink_in(merged_source_sink_in)
    else:
        raise NotImplementedError(
            'merged_source_sink_in can either be .yaml or .in file')
    # generate cols based on the modules and the order of sources defined in merged_source_sink_in.   
    sites = df_sources.index
  
    cols = [["%s_%s"%(v,s) for s in sites] for v in var]
    cols = np.concatenate(cols)
    # frequency check
    if freq=='infer':
        if th1.index.freq!=th2.index.freq:
            print("th1 and th2 has different frequency")
    else:
        th1 = th1.asfreq(freq)
        th2 = th2.asfreq(freq)    
    # perform merging 
    th_merged = th1.join(th2,how=how,rsuffix='r').drop(columns=['datetimer'])

    # if an item of cols is not in the columns of th. 
    colm = np.asarray(list(set(cols) - set(th_merged.columns.values)))
    colm = [str(s) for s in colm]    
    th_merged[colm] = pd.DataFrame(np.ones([len(th_merged),len(colm)])*-9999.0)   
    
    th_merged = th_merged.fillna(-9999.0)
    
    cols = np.append(['datetime'],cols)
    th_merged = th_merged[cols] #rearrange the array to have the same order as defined in merged_source_sink_in
    th_merged['datetime'] = np.datetime_as_string(th_merged.index.values,'h')
    write_source_sink_csv(th_merged,csv_merged)    

def concat_vsource_sink_csv(csv_fn1,csv_fn2,merged_source_sink_in,
                            csv_type,csv_merged,freq='infer',how='left'):
    
    """ Concatenate source and sink files
    
    Parameters
    ----------
    csv_fn1 : STR
        Filename of the 1st dated vsource.csv or vsink.csv
    csv_fn2 : STR
        Filename of the 2nd dated vsource.csv or vsink.csv
    merged_source_sink_in : STR
        Filename of source_sink.in or source_sink.yaml
    csv_type : STR
        Indicate if the files are "sources" or "sinks"
    csv_merged : STR
        Output merged csv filename
    freq : STR, optional
        Frequency for the output csv file. The default is 'infer'.
    how : STR, optional
        Options see pandas.join(). The default is 'left'.

    Raises
    ------
    NotImplementedError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # merged_source_sink_in: the merged source_sink.in or source_sink.yaml file 
    # where the data sources are from csv_fn1, csv_fn2. 
    if merged_source_sink_in.endswith('yaml'):
        df_sources,df_sinks = read_source_sink_yaml(merged_source_sink_in)
    elif merged_source_sink_in.endswith('in'):
        df_sources,df_sinks = read_source_sink_in(merged_source_sink_in)
    else:
        raise NotImplementedError(
            'merged_source_sink_in can either be .yaml or .in file')
    if csv_type == 'sources':
        sites = df_sources.index
    elif csv_type == 'sink':
        sites = df_sinks.index
    else:
        raise NotImplementedError('csv_type can either be sources or sinks')
    th1 = read_source_sink_csv(csv_fn1)
    th2 = read_source_sink_csv(csv_fn2)
    if freq=='infer':
        if th1.index.freq!=th2.index.freq:
            print("th1 and th2 has different frequency")
    else:
        th1 = th1.asfreq(freq)
        th2 = th2.asfreq(freq)
    th_merged = th1.join(th2,how=how,rsuffix='r').drop(columns=['datetimer'])
    th_merged = th_merged.fillna(-9999.0)
    cols = np.append(['datetime'],sites)
    th_merged = th_merged[cols] #rearrange the array to have the same order as defined in merged_source_sink_in
    th_merged['datetime'] = np.datetime_as_string(th_merged.index.values,'h')
    write_source_sink_csv(th_merged,csv_merged)   
    
def read_source_sink_th(th_fn, source_or_sink_df):
    """ Read source_sink.th file
    
    Parameters
    ----------
    th_fn : STR
        \*.th file
    source_or_sink_df : PANDAS DATAFRAME
        source_df or sink_df produced by read_source_sink_in.

    Returns
    -------
    th : PANDAS DATAFRAME
        source or sink data frame with time series data of T, S and other variables
    """
    th = pd.read_csv(th_fn,header=None, delimiter=r"\s+") 
    if len(th.columns) == len(source_or_sink_df) +1:
        col_names = ['seconds'] + list(source_or_sink_df.name.values)         
    elif len(th.columns) == len(source_or_sink_df)*2 +1:  # T and S
        T_col = ['T_%s'%s for s in source_or_sink_df.name.values]
        S_col = ['S_%s'%s for s in source_or_sink_df.name.values]
        col_names = ['seconds'] + T_col + S_col
    else: # if more sources are involved
        nsources = int((len(th.columns)-1)/len(source_or_sink_df))-2
        T_col = ['T_%s'%s for s in source_or_sink_df.name.values]
        S_col = ['S_%s'%s for s in source_or_sink_df.name.values]     
        B_col = []
        for n in range(nsources):
            B_col +=['B%d_%s'%(n+1,s) for s in source_or_sink_df.name.values]
        col_names = ['seconds'] + T_col + S_col + B_col
    th.columns = col_names
    return th  

def read_source_sink_in(source_sink_in):  
    """ Parse source sink.in file 
    
    Parameters
    ----------
    source_sink_in : STR
        source_sink.in
    Returns
    -------
    source_df : PANDAS DATAFRAME
        a dataframe of source names
    sink_df : TYPE
        a dataframe of sink names
    """      
    
    with open(source_sink_in,'r') as file:
        lines = file.readlines()        
    nsource = 0
    nsink=0
    for l in lines:
        if 'total # of elems with sources' in l:
            nsource = int(l.split()[0])
            source_ele = []
            source_name = []
            source_from = []
        elif 'total # of elems with sinks' in l:
            nsink = int(l.split()[0])
            sink_ele = []
            sink_name = []
            sink_from = []
        elif nsource >0 and len(source_ele)<nsource:
            source_ele.append(int(l.split()[0]))
            source_name.append(l.split()[2])
            if 'delta' in source_name[-1]:
                source_from.append('delta')
            elif 'suisun' in source_name[-1]:
                source_from.append('suisun')
            elif 'dicu' in source_name[-1]:
                source_from.append('dicu')
            elif 'potw' in source_name[-1]:
                source_from.append('potw')
            else:
                source_from.append('unknown')
        elif nsink>0 and len(sink_ele)<nsink:
            sink_ele.append(int(l.split()[0]))
            sink_name.append(l.split()[2])
            if 'delta' in sink_name[-1]:
                sink_from.append('delta')
            elif 'suisun' in sink_name[-1]:
                sink_from.append('suisun')   
    
    assert(nsource==len(source_ele))
    assert(nsink==len(sink_ele))
    
    source_df = pd.DataFrame({'element': source_ele,
                              'name': source_name,
                              'source': source_from})
    sink_df = pd.DataFrame({'element': sink_ele,
                             'name': sink_name,
                             'source': sink_from})
    # setting index will reorder the df, so this should not be implemented
    #source_df = source_df.set_index('name')
    #sink_df = sink_df.set_index('name')
    return source_df, sink_df             

    
def write_source_sink_in(source_sink_yaml, hgrid_fn, 
                         source_sink_in='source_sink.in'):
    """ Create source_sink.in based on hgrid and source_sink.yaml using the 
    create_source_sink_in function in schimpy preprocessor
    
    Parameters
    ----------
    source_sink_yaml : STR
        Filename for source_sink.yaml
    hgrid_fn : STR
        schism hgrid.gr3
    source_sink_in : STR, optional
        Filename for source_sink.in. The default is 'source_sink.in'.

    Returns
    -------
    None.

    """
    with open(source_sink_yaml, 'r') as file:
        source_sink = yaml.safe_load(file)
    s = create_schism_setup(hgrid_fn)
    s.create_source_sink_in(source_sink,source_sink_in)

def yaml2csv(source_yaml,source_csv): 
    """ Converting source yaml file to source csv file
    
    Parameters
    ----------
    source_yaml :YAML FILENAME 
        DESCRIPTION.
    source_csv : CSV FILENAME
        DESCRIPTION.

    Returns
    -------
    None.

    """
    with open(source_yaml, 'r') as file:
        source_sink = yaml.safe_load(file)
    potw = source_sink['sources']
    sites = potw.keys()
    sites = list(sites)
    a = [potw[k] for k in sites]
    a = np.array(a)
    x = a[:,0]
    y = a[:,1]
    df =  pd.DataFrame({'sites': sites,'x':x,'y':y})
    df.to_csv(source_csv)

def csv2yaml(source_csv,source_yaml):
    """ Converting from source csv to source yaml file
    
    Parameters
    ----------
    source_csv : CSV FILENAME
        DESCRIPTION.
        
    source_yaml : YAML FILENAME
        DESCRIPTION.

    Returns
    -------
    None

    """
    csv_data = pd.read_csv(source_csv)
    csv_data = csv_data[['site','utm_x','utm_y']]
    
    dict_file = {}
    for r in csv_data.iterrows():
        dict_file[r[1].site] = [float(r[1].utm_x),float(r[1].utm_y)]
    dict_file = {'sources':dict_file}
    with open(source_yaml, 'w') as file:
        documents = yaml.dump(dict_file, file)
if __name__ == "__main__":
    yaml_fn1 = 'source_sink.yaml'
    yaml_fn2 = 'potw_sources.yaml'
    new_yaml = 'source_sink_new.yaml'
    csv_fn1 = 'msource_formatted.csv'
    csv_fn2 = 'msource_potws.csv'
    csv_v_fn1 = 'vsource_formatted.csv'
    csv_v_fn2 = 'vsource_potws.csv'
    csv_merged = 'msource_new.csv'
    csv_merged_v = 'vsource_new.csv'
    hgrid = 'hgrid.gr3'
    concat_source_sink_yaml(yaml_fn1,yaml_fn2,new_yaml)
    concat_msource_csv(csv_fn1,csv_fn2,new_yaml,
                       csv_merged,freq='infer',how='left')
    concat_vsource_sink_csv(csv_v_fn1,csv_v_fn2,new_yaml,
                            'sources',csv_merged,freq='infer',how='left')
    write_source_sink_in(new_yaml,hgrid_fn)    