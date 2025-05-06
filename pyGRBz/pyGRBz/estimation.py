def stats(path_dir,ext_laws=['smc', 'lmc', 'mw', 'nodust','sne'],lim_bic = 2):
    """Compare detection if Multiple targets
       Best fit in chi-square
       Best fit in BIC within the threshold
    """    
    import numpy as np
    import os
    from astropy.io import ascii
    from astropy.table import Table, vstack
    import pandas as pd
    fit_results_smc = None
    fit_results_lmc = None
    fit_results_nodust = None
    fit_results_sne = None
    for ext in ext_laws:
        if ext == 'smc':
            fit_results_smc=ascii.read(path_dir+'best_fits_%s.dat' % ext)
        elif ext == 'mw':
            fit_results_mw=ascii.read(path_dir+'best_fits_%s.dat' % ext)
        elif ext == 'lmc':
            fit_results_lmc=ascii.read(path_dir+'best_fits_%s.dat' % ext)
        elif ext == 'nodust':
            fit_results_nodust=ascii.read(path_dir+'best_fits_%s.dat' % ext)
        elif ext == 'sne':
            fit_results_sne=ascii.read(path_dir+'best_fits_%s.dat' % ext)

    GRB_list = None
    if fit_results_smc is not None:
        GRB_list = fit_results_smc['name']
    elif fit_results_lmc is not None:
        GRB_list = fit_results_lmc['name']
    elif fit_results_nodust is not None:
        GRB_list = fit_results_nodust['name']
    elif fit_results_sne is not None:
        GRB_list = fit_results_sne['name']
    
    if False:
        print ('Not same number of detected GRBs')
    else:
        new_table=[]
        new_table_bic=[]
        new_table_conflict=[]
        nb_GRB_bic_det = 0
        nb_GRB_agreed = 0
    
        for GRB in GRB_list:
            if 'smc' in ext_laws:
                mask_smc = fit_results_smc['name'] == GRB
                sum_proba_smc = float(fit_results_smc['sum_proba'][mask_smc])
                bic_smc = float(fit_results_smc['bic'][mask_smc])
                Av_smc = float(fit_results_smc['best_Av'][mask_smc])
            else:
                mask_smc = []
                sum_proba_smc = np.nan
                bic_smc = 1e10
                Av_smc = np.nan
            if 'mw' in ext_laws:
                mask_mw = fit_results_mw['name'] == GRB
                sum_proba_mw = float(fit_results_mw['sum_proba'][mask_mw])
                bic_mw = float(fit_results_mw['bic'][mask_mw])
                Av_mw = float(fit_results_mw['best_Av'][mask_mw])
            else:
                mask_mw = []
                sum_proba_mw = np.nan
                bic_mw = 1e10
                Av_mw = np.nan
    
            if 'lmc' in ext_laws:
                mask_lmc = fit_results_lmc['name'] == GRB
                sum_proba_lmc = float(fit_results_lmc['sum_proba'][mask_lmc])
                bic_lmc = float(fit_results_lmc['bic'][mask_lmc])
                Av_lmc = float(fit_results_lmc['best_Av'][mask_lmc])
            else:
                mask_lmc = []
                sum_proba_lmc =np.nan
                bic_lmc = 1e10
                Av_lmc = np.nan
            if 'nodust' in ext_laws:
                mask_nodust = fit_results_nodust['name'] == GRB
                sum_proba_nodust = float(fit_results_nodust['sum_proba'][mask_nodust])
                bic_nodust = float(fit_results_nodust['bic'][mask_nodust])
                Av_nodust = float(fit_results_nodust['best_Av'][mask_nodust])
            else:
                mask_nodust = []
                sum_proba_nodust =np.nan
                bic_nodust = 1e10
                Av_nodust = np.nan
            if 'sne' in ext_laws:
                mask_sne = fit_results_sne['name'] == GRB
                sum_proba_sne = fit_results_sne['sum_proba'][mask_sne]
                bic_sne = fit_results_sne['bic'][mask_sne]
                Av_sne = fit_results_sne['best_Av'][mask_sne]
            else:
                mask_sne = []
                sum_proba_sne = np.nan
                bic_sne = 1e10
                Av_sne = np.nan                
            
            z_GRB_chi2 = 0
    
            if ('smc' in ext_laws) and (np.nanmax([sum_proba_smc,sum_proba_sne,sum_proba_lmc,sum_proba_mw,sum_proba_nodust]) == sum_proba_smc) :
                #print (GRB,'smc')
                choice_GRB_chi2 = 'smc'
                new_table.append(fit_results_smc[mask_smc])
                z_GRB_chi2 = fit_results_smc[mask_smc]['z_sim']
            elif ('mw' in ext_laws) and  (np.nanmax([sum_proba_smc,sum_proba_sne,sum_proba_lmc,sum_proba_mw,sum_proba_nodust]) == sum_proba_mw) :
                #print (GRB,'mw')
                choice_GRB_chi2 = 'mw'
                new_table.append(fit_results_mw[mask_mw])
                z_GRB_chi2 = fit_results_mw[mask_mw]['z_sim']
            elif ('lmc' in ext_laws) and  (np.nanmax([sum_proba_smc,sum_proba_sne,sum_proba_lmc,sum_proba_mw,sum_proba_nodust]) == sum_proba_lmc) :
                #print (GRB,'lmc')
                choice_GRB_chi2 = 'lmc'
                new_table.append(fit_results_lmc[mask_lmc])
                z_GRB_chi2 = fit_results_lmc[mask_lmc]['z_sim']
            elif ('nodust' in ext_laws) and  (np.nanmax([sum_proba_smc,sum_proba_sne,sum_proba_lmc,sum_proba_mw,sum_proba_nodust]) == sum_proba_nodust) :
                #print (GRB,'lmc')
                choice_GRB_chi2 = 'nodust'
                new_table.append(fit_results_nodust[mask_nodust])
                z_GRB_chi2 = fit_results_nodust[mask_nodust]['z_sim']
            elif ('sne' in ext_laws) and  (np.nanmax([sum_proba_smc,sum_proba_sne,sum_proba_lmc,sum_proba_mw,sum_proba_nodust]) == sum_proba_sne) :
                #print (GRB,'lmc')
                choice_GRB_chi2 = 'sne'
                new_table.append(fit_results_sne[mask_sne])
                z_GRB_chi2 = fit_results_sne[mask_sne]['z_sim']
                
    
    ################# DELTABIC limit        
            
            lim_bic = 2 # Threshold for bic comparison
    
            count_bic = 0
            list_bic = [bic_smc, bic_mw, bic_lmc, bic_nodust]
            choice_GRB_bic = ''
            z_GRB = 0
    
            if Av_mw == 0 and Av_smc == 0 and Av_lmc == 0 and Av_sne == 0 and Av_nodust == 0:
                new_table_bic.append(fit_results_nodust[mask_nodust])
                count_bic = "no_dust_condition"
            else:
                if ('smc' in ext_laws) and bic_sne-bic_smc>lim_bic and bic_lmc-bic_smc>lim_bic and bic_mw-bic_smc>lim_bic and bic_nodust-bic_smc>lim_bic :
                    new_table_bic.append(fit_results_smc[mask_smc])
                    nb_GRB_bic_det += 1
                    choice_GRB_bic = 'smc'
                    delta_bic_min = np.nanmin([bic_lmc-bic_smc, bic_sne-bic_smc,bic_mw-bic_smc, bic_nodust-bic_smc])
                    count_bic = 3
                    z_GRB = fit_results_smc[mask_smc]['z_sim']
                elif ('mw' in ext_laws) and bic_smc-bic_mw>lim_bic and bic_sne-bic_mw>lim_bic and bic_lmc-bic_mw>lim_bic and bic_nodust-bic_mw>lim_bic:
                    new_table_bic.append(fit_results_mw[mask_mw])
                    nb_GRB_bic_det += 1
                    choice_GRB_bic = 'mw'
                    delta_bic_min = np.nanmin([bic_smc-bic_mw,bic_sne-bic_mw,bic_lmc-bic_mw, bic_nodust-bic_mw])
                    count_bic = 3
                    z_GRB = fit_results_mw[mask_mw]['z_sim']
                elif ('sne' in ext_laws) and bic_smc-bic_sne>lim_bic and bic_lmc-bic_sne>lim_bic and bic_nodust-bic_sne>lim_bic and bic_mw-bic_sne>lim_bic:
                    new_table_bic.append(fit_results_sne[mask_sne])
                    nb_GRB_bic_det += 1
                    choice_GRB_bic = 'sne'
                    delta_bic_min = np.nanmin([bic_smc-bic_sne, bic_lmc-bic_sne, bic_nodust-bic_sne,bic_mw-bic_sne])
                    count_bic = 3
                    z_GRB = fit_results_sne[mask_sne]['z_sim']
                elif ('lmc' in ext_laws) and bic_sne-bic_lmc>lim_bic and bic_smc-bic_lmc>lim_bic and bic_mw-bic_lmc>lim_bic and bic_nodust-bic_lmc>lim_bic:
                    new_table_bic.append(fit_results_lmc[mask_lmc])
                    nb_GRB_bic_det += 1
                    choice_GRB_bic = 'lmc'
                    delta_bic_min = np.nanmin([bic_sne-bic_lmc,bic_smc-bic_lmc, bic_mw-bic_lmc, bic_nodust-bic_lmc])
                    count_bic = 3
                    z_GRB = fit_results_lmc[mask_lmc]['z_sim']
                elif ('nodust' in ext_laws) and bic_sne-bic_nodust and bic_smc-bic_nodust>lim_bic and bic_mw-bic_nodust>lim_bic and bic_lmc-bic_nodust>lim_bic:
                    new_table_bic.append(fit_results_nodust[mask_nodust])
                    nb_GRB_bic_det += 1
                    choice_GRB_bic = 'nodust'
                    delta_bic_min = np.nanmin([bic_smc-bic_nodust,bic_sne-bic_nodust, bic_mw-bic_nodust, bic_lmc-bic_nodust])
                    count_bic = 3
                    z_GRB = fit_results_nodust[mask_nodust]['z_sim']
                
                choice_GRB_conflict_1 = ''
                choice_GRB_conflict_2 = ''
                
                #### list_bic = [bic_smc, bic_mw, bic_lmc, bic_nodust]
                if count_bic == 0 and count_bic != "no_dust_condition":
                    min_bic_1 = np.nanmin(list_bic)
                    min_bic_idx_1 = list_bic.index(min_bic_1)
                    list_bic_reduced = list_bic[:min_bic_idx_1] + list_bic[min_bic_idx_1+1:]
                    min_bic_2 = np.nanmin(list_bic_reduced)
                    min_bic_idx_2 = list_bic.index(min_bic_2)
                    if ('smc' in ext_laws) and min_bic_idx_1 == 0:
                        new_table_conflict.append(fit_results_smc[mask_smc])
                        choice_GRB_conflict_1 = 'smc'
                    elif ('mw' in ext_laws) and min_bic_idx_1 == 1:
                        new_table_conflict.append(fit_results_mw[mask_mw])
                        choice_GRB_conflict_1 = 'mw'
                    elif ('sne' in ext_laws) and min_bic_idx_1 == 1:
                        new_table_conflict.append(fit_results_sne[mask_sne])
                        choice_GRB_conflict_1 = 'sne'
                    elif ('lmc' in ext_laws) and min_bic_idx_1 == 2:
                        new_table_conflict.append(fit_results_lmc[mask_lmc])
                        choice_GRB_conflict_1 = 'lmc'
                    elif ('nodust' in ext_laws) and min_bic_idx_1 == 3:
                        new_table_conflict.append(fit_results_nodust[mask_nodust])
                        choice_GRB_conflict_1 = 'nodust'
                    if ('smc' in ext_laws) and min_bic_idx_2 == 0:
                        new_table_conflict.append(fit_results_smc[mask_smc])
                        choice_GRB_conflict_2 = 'smc'
                    elif ('mw' in ext_laws) and min_bic_idx_2 == 1:
                        new_table_conflict.append(fit_results_mw[mask_mw])
                        choice_GRB_conflict_2 = 'mw'
                    elif ('sne' in ext_laws) and min_bic_idx_2 == 1:
                        new_table_conflict.append(fit_results_sne[mask_sne])
                        choice_GRB_conflict_2 = 'sne'
                    elif ('lmc' in ext_laws) and min_bic_idx_2 == 2:
                        new_table_conflict.append(fit_results_lmc[mask_lmc])
                        choice_GRB_conflict_2 = 'lmc'
                    elif ('nodust' in ext_laws) and min_bic_idx_2 == 3:
                        new_table_conflict.append(fit_results_nodust[mask_nodust])
                        choice_GRB_conflict_2 = 'nodust'
                        
            if choice_GRB_chi2 == choice_GRB_bic:
                nb_GRB_agreed += 1
    new_table=vstack(new_table)
    new_table.write(path_dir+'best_fits_combined.dat',format='ascii', overwrite=True)
    print("Best fit is given using the ",list(new_table['ext_law'][new_table['name']==list(set(GRB_list))])," extinction law")
    
    if new_table_bic != []:
        new_table_bic=vstack(new_table_bic)
    else:
        new_table_bic = Table(new_table_bic)
    new_table_bic.write(path_dir+'best_fits_combined_bic.dat',format='ascii', overwrite=True) 
    
    if new_table_conflict != []:
        new_table_conflict=vstack(new_table_conflict)
    else:
        new_table_conflict = Table(new_table_conflict)
    new_table_conflict.write(path_dir+'best_fits_combined_conflict.dat',format='ascii', overwrite=True)

    if nb_GRB_bic_det>1:
        percent_det_bic = nb_GRB_bic_det / len(GRB_list) * 100
        print('total GRBs resolved by BIC: ', percent_det_bic, '%')
        percent_agreed = nb_GRB_agreed / len(GRB_list) * 100
        print('nb GRB agreed chi2 - BIC: ', percent_agreed, '%')
    else:
        print("GRB not resolved with BIC threshold of ",lim_bic)
