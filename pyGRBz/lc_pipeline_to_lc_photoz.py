def alid_to_lc(input_path,output_path,ra,dec,mw_correction=1,save_file = True):
    """
    Transform ALID LC to readable LC photo-z readable
    """
    import pandas as pd    
    import numpy as np
    import os
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
        filename = os.path.basename(input_path)
        obj_name = os.path.splitext(filename)[0]
            
            #Open file to write header first, then the CSV data
        with open(output_path, 'w') as f:
            f.write(f"#name:{obj_name}\n")
            f.write(f"#type:lc\n")
            f.write(f"#RA_J2000:{ra}\n")
            f.write(f"#DEC_J2000:{dec}\n")
            f.write(f"#MW_corrected:{mw_correction}\n") #Default is not corrected
            f.write(f"#time_unit:s\n")
            f.write(f"#z:-99\n")
            f.write(f"#Av_host:-99\n")
            
            # Write the dataframe content
            data.to_csv(f, index=False, sep=" ")
