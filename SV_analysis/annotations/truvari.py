#!/usr/bin/env python3

import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

chrs = ['chr1','chr2','chr3', 'chr4', 'chr5',
        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY']

genome='grch38'

file_path = str(sys.argv[1])
out_path= str(sys.argv[2])



if "alignment" in file_path:
    merged_tmp = joblib.load(file_path)
    merged_tmp = merged_tmp[merged_tmp.chrom.isin(chrs)]
    merged_tmp['filter'] = [','.join(map(str, l)) for l in merged_tmp['filter']]
    merged_tmp = merged_tmp[((merged_tmp['filter']=='PASS') | (merged_tmp["svtype"]=='INV') | (merged_tmp["svtype"]=='DUP')) & ~(merged_tmp["svtype"]=='UNK') & (merged_tmp["svlen"]>=50)]

    if "cutesv" or "sniffles2":
        sample_columns = [col for col in merged_tmp.columns if '_DV' in col]       
    elif "pbsv":
        sample_columns = [col for col in merged_tmp.columns if '_AD_ref' in col]
    elif "svdss":
        sample_columns = [col for col in merged_tmp.columns if '_GT' in col]
        merged_tmp[sample_columns] = merged_tmp[sample_columns].astype(str)
        merged_tmp = merged_tmp.replace("(None, None)", np.nan)
            
    c = merged_tmp[sample_columns].columns.to_numpy()
    c_tmp = []
    for DV in c:
        c_tmp.append(DV.split('.')[0])
    c_dict = dict(zip(c, c_tmp))

    caller=f"{DV.split('.')[2]}".split("_")[0]

    res = []
    sample_count=[]
    for x in merged_tmp[sample_columns].notna().to_numpy():
        res_temp1 = c[x].tolist()
        res_temp2 = []
        for i in res_temp1:
            res_temp2.append(c_dict[i])
        names=list(set(res_temp2))
        res.append(names)
        sample_count.append(len(names))

    tmp = pd.DataFrame(
        {
         'sample_count': sample_count,
         "sample_names":res
        },
        index=merged_tmp.index)

    tmp_list=[merged_tmp, tmp]
    merged = pd.concat(tmp_list, axis=1)
    merged['caller'] = caller

elif "assembly" in file_path:
    merged_tmp = joblib.load(file_path)
    merged_tmp = merged_tmp[merged_tmp.chrom.isin(chrs)]
    merged_tmp['filter'] = [','.join(map(str, l)) for l in merged_tmp['filter']]
    merged_tmp = merged_tmp[((merged_tmp['filter']=='PASS') | (merged_tmp["svtype"]=='INV') | (merged_tmp["svtype"]=='DUP')) & ~(merged_tmp["svtype"]=='UNK') & (merged_tmp["svlen"]>=50)]
    
    sample_columns = [col for col in merged_tmp.columns if '_GT' in col]
    merged_tmp[sample_columns] = merged_tmp[sample_columns].astype(str)
    merged_tmp2 = merged_tmp.replace("(None, None)", np.nan)
    
    c = merged_tmp2[sample_columns].columns.to_numpy()
    ids=[]
    for id in c:
        id_tmp=id.split("_GT")[0]
        if ":" in id_tmp:
           id_tmp=id_tmp.split(":")[1]
        ids.append(id_tmp)
    c_dict = dict(zip(c, np.array(ids)))
    res = []
    sample_count = []
    
    for x in merged_tmp2[sample_columns].notna().to_numpy():
        res_temp1 = c[x].tolist()
        res_temp2 = []
        for i in res_temp1:
            res_temp2.append(c_dict[i])
        names=list(set(res_temp2))
        res.append(names)
        sample_count.append(len(names))
    tmp = pd.DataFrame(
        {
         'sample_count':sample_count,
         "sample_names":res
        },
        index=merged_tmp2.index)




tmp_list=[merged_tmp, tmp]
merged = pd.concat(tmp_list, axis=1)


filename = file_path.split('/')[-1]
filename = filename.split('.')[0]

pheno_dict_tmp={
    'chrom':merged['chrom'],
    'start':merged['start'],
    'end':merged['end'],
    'id':merged.index,
    'svtype':merged['svtype'],
}
# print(pheno_dict_tmp)
pheno_sv_tmp = pd.DataFrame(pheno_dict_tmp)


pheno_sv_tmp = pheno_sv_tmp.reset_index(drop=True)

pheno_sv_tmp = pheno_sv_tmp.replace('DEL', 'deletion')
pheno_sv_tmp = pheno_sv_tmp.replace('INS', 'insertion')
pheno_sv_tmp = pheno_sv_tmp.replace('DUP', 'duplication')
pheno_sv_tmp = pheno_sv_tmp.replace('INV', 'inversion')

pheno_sv_tmp.to_csv(f'{out_path}/{filename}.bed', sep='\t', index=False, header=False)
