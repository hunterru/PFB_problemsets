#!/usr/bin/env python3


# sort_by can be either 'Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', or 'MIT_specificity'

def range_overlap(user_bed_path, guide_loc_path, output_name, upstream , downstream , sort_by, de_novo, cloning_strategy):
   
    #note: if de_novo = True, then sort_by is automatically set to be "mismatch_score"
    import sys
    import pyranges as pr
    import pandas as pd
    from targetsite_to_primers import revcomp, startG, cloning_parameters, targetsite_to_primers

    #######################
    # read in the bed file from the user (right now this will be a file of TSS's for a specific gene list from Mina's function)
    #######################
    
    
    user_bed = pd.read_csv(user_bed_path, sep = '\t', header = None) # important: this assumes that the bed file has no column names

    
    # need to have first three columns be called 'Chromosome', 'Start', 'End' for a pyRanges object so here we will change the column names
    column_names_user = user_bed.columns.values
    column_names_user_list = list(column_names_user)
    column_names_user_list_str = [str(i) for i in column_names_user_list]

    column_names_user_list_str[0:6] = ['Chromosome', 'Start', 'End', 'Gene', '5', "Strand"]
    user_bed.columns = column_names_user_list_str
    
    # iterate over the rows and change the start and end positions in the user bed file based on the upstream and downstream arguments
    # also et the start column as a list to use later to determine the distance from TSS
    user_bed_start = [] 
    for index, row in user_bed.iterrows():
        if user_bed.at[index, 'Strand'] == '-':
            user_bed_start.append(user_bed.at[index, 'Start'])
            if user_bed.at[index, 'Start'] < downstream:
                user_bed.at[index, 'Start'] = 0
                user_bed.at[index, 'End'] += upstream
            else:   
                user_bed.at[index, 'End'] += upstream
                user_bed.at[index, 'Start'] -= downstream
        else:
            user_bed_start.append(user_bed.at[index, 'End'])
            if user_bed.at[index, 'Start'] < upstream:
                user_bed.at[index, 'Start'] = 0
                user_bed.at[index, 'End'] += downstream
            else:
                user_bed.at[index, 'Start'] -= upstream
                user_bed.at[index, 'End'] += downstream
                
    user_bed = user_bed.assign(Original_start = user_bed_start)
    
    user_bed_pyR = pr.PyRanges(user_bed) # convert the panda df to a pyRanges object which is required for the overlap function
    #user_bed_pyR_merge = user_bed_pyR.merge() # if we wanted to collapse overlapping ranges, we could use this, but removed it becuase we lose gene column when we do it
    

    # these are the guides already determined for human genome
    guide_locs = pd.read_csv(guide_loc_path, sep = '\t') # important: this assumes that the guide table has column names, need to change if this isnt true!
    column_names_gloc = guide_locs.columns.values
    column_names_gloc_list = list(column_names_gloc)
    column_names_gloc_list_str = [str(i) for i in column_names_gloc_list]

    column_names_gloc_list_str[0:3] = ['Chromosome', 'Start', 'End']
    guide_locs.columns = column_names_gloc_list_str
    guide_locs_noNaN = guide_locs.fillna(0)
    if de_novo == False:
        guide_locs_noNaN[['Doench2016_perc','Doench2016_score']] = guide_locs_noNaN.fusi.str.split('%', expand=True) 
        guide_locs_noNaN['Doench2016_score'] = guide_locs_noNaN['Doench2016_score'].str.replace(r"[\(\)]","") 
        guide_locs_noNaN[['Moreno_Matos_perc','Moreno_Matos_score']] = guide_locs_noNaN.crisprScan.str.split('%', expand=True)
        guide_locs_noNaN['Moreno_Matos_score'] = guide_locs_noNaN['Moreno_Matos_score'].str.replace(r"[\(\)]","")
        guide_locs_noNaN[['MIT_specificity']] = guide_locs_noNaN[['scoreDesc']]
        guide_locs_noNaN['MIT_specificity'] = guide_locs_noNaN['MIT_specificity'].str.replace(r"[A-Za-z\s\.\-]+","0")

        guide_locs_noNaN[["Doench2016_perc", "Doench2016_score", "Moreno_Matos_perc", "Moreno_Matos_score" ,"MIT_specificity"]] = guide_locs_noNaN[["Doench2016_perc", "Doench2016_score", "Moreno_Matos_perc", "Moreno_Matos_score" ,"MIT_specificity"]].apply(pd.to_numeric)
        
        guide_locs_noNaN_select = guide_locs_noNaN.iloc[ : ,[0,1,2,3,4,5,6,7,11,20,21,22,23,24]]
        
        if sort_by not in ['Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', 'MIT_specificity']:
            print("The sort_by argument must be 'Doench2016_perc', 'Doench2016_score', 'Moreno_Matos_perc', 'Moreno_Matos_score', 'MIT_specificity', because you failed to comply the default ('Doench2016_perc') will be used")
            sort_by = "Doench2016_perc"
        
        guide_locs_pyR = pr.PyRanges(guide_locs_noNaN_select)
        
    else: 
        sort_by = "mismatch_score"
        guide_locs_pyR = pr.PyRanges(guide_locs_noNaN)
    
    class NoOverlapError(Exception):
        pass
    
    if len(guide_locs_pyR.overlap(user_bed_pyR)) == 0:
        raise NoOverlapError("There are no overlaps between the user supplied ranges and the gRNAs!")
        
    guide_locs_pyR_overlap = guide_locs_pyR.overlap(user_bed_pyR)
    
    for k, guide_locs_pyR_overlap_df in guide_locs_pyR_overlap: # convert to pandas with this loop to more easily manipulate the df
        guide_locs_pyR_overlap_df   
    
    gene_list = []
    user_bed_start_ol = []
    for index_ol, row_ol in guide_locs_pyR_overlap_df.iterrows():
        for index_ub, row_ub in user_bed.iterrows():
            if row_ol[0] == row_ub[0]:
                if row_ol[1] in range(row_ub[1], (row_ub[2] + 1)) or row_ol[2] in range(row_ub[1], (row_ub[2] + 1)):
                    gene_list.append(row_ub[3])
                    user_bed_start_ol.append(row_ub[6])

    
    guide_locs_pyR_overlap_df = guide_locs_pyR_overlap_df.assign(Gene = gene_list)
    guide_locs_pyR_overlap_df = guide_locs_pyR_overlap_df.assign(Original_Start = user_bed_start_ol)
    guide_locs_pyR_overlap_df_sort = guide_locs_pyR_overlap_df.sort_values(["Gene", sort_by], ascending = [True, False])
    
    primer1_list = []
    primer2_list = []
    distance_from_tss = []
    for index, row in guide_locs_pyR_overlap_df_sort.iterrows():
        # get a list of distances from the end og each guide to the original TSS
        distance_temp = guide_locs_pyR_overlap_df_sort.at[index, 'End'] - guide_locs_pyR_overlap_df_sort.at[index, 'Original_Start']
        distance_from_tss.append(distance_temp)
        
        # use Minas functions to get lists of each primer
        gRNA = guide_locs_pyR_overlap_df_sort.at[index, 'guideSeq']
        (Grequired, p1specs, p2specs) = cloning_parameters(cloning_strategy)
        (primer1, primer2) = targetsite_to_primers(gRNA, cloning_strategy, Grequired, p1specs, p2specs)
        primer1_list.append(primer1)
        primer2_list.append(primer2)
    
    # add the various columns to the df that we populated in the above for loop
    guide_locs_pyR_overlap_df_sort = guide_locs_pyR_overlap_df_sort.assign(GuideEnd_to_TSS = distance_from_tss)
    guide_locs_pyR_overlap_df_sort = guide_locs_pyR_overlap_df_sort.drop(columns = ["Original_Start"])
    
    guide_locs_pyR_overlap_df_sort = guide_locs_pyR_overlap_df_sort.assign(Primer1 = primer1_list)
    guide_locs_pyR_overlap_df_sort = guide_locs_pyR_overlap_df_sort.assign(Primer2 = primer2_list)
    col_name_p1 = 'Primer1 - ' + cloning_strategy
    col_name_p2 = 'Primer2 - ' + cloning_strategy
    guide_locs_pyR_overlap_df_sort = guide_locs_pyR_overlap_df_sort.rename(columns = {"Primer1" : col_name_p1, "Primer2" : col_name_p2})
    
    guide_locs_pyR_overlap_df_sort.to_csv((output_name + '.txt'), sep = '\t', header = True)
    
    return(guide_locs_pyR_overlap_df_sort)

	
