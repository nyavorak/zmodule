def alid_to_lc(input_path,output_path,save_file = True):
    """
    Transform ALID LC to readable LC photo-z readable
    """
    import pandas as pd
    #"/home/nyavorak/pyGRBz/notebooks/data/lc/ALID_20251111051015_cleaned.csv"
    
    #When the path from and for pipeline will be defined
    #alid or grb_name
    #input_path = 
    #output_path = 
    
    df = pd.read_csv(input_path)   
    df.columns = [col.lstrip() for col in df.columns]
    
    data=pd.DataFrame({})
    data ['time_since_burst']=[1]*df.shape[0]
    data ['Texp']=[1]*df.shape[0]
    data ['band']=[1]*df.shape[0]
    data ['zp']=["-"]*df.shape[0]
    data ['flux_err']=[1]*df.shape[0]
    data ['telescope']=["colibri"]*df.shape[0]
    data ['flux']=[1]*df.shape[0]
    data ['detection']=[1]*df.shape[0]
    data ['flux_unit']=["AB"]*df.shape[0]
    
    
    for a,b in enumerate(df["start_time_MJD"]):
        t = abs((df["Burst time (MJD)"][0]-((df["start_time_MJD"][a]+df["end_time_MJD"][a])/2))*(3600*24))
        data ['time_since_burst'][a] = t
        data ['Texp'][a] = df["N_images"][a]*60
        data ['band'][a] = df["filter"][a]    
        if np.isnan(df["mag"][a]):
            data['flux'][a] = df["3sigma upper limit"][a]
            data ['flux_err'][a] = 0.1
            data['detection'][a] = 0
        else:
            data['flux'][a] = df["mag"][a]
            data ['flux_err'][a] = df["magerr"][a]
            data['detection'][a] = 1
    data = data[['time_since_burst', 'band', 'flux', 'flux_err', 'Texp', 'zp', 'flux_unit', 'detection', 'telescope']]
    
    if save_file:
        data.to_csv(output_path,index=False,sep=" ")