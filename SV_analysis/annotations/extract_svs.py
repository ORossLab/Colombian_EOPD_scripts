import pandas as pd
chrs = ['chr1','chr2','chr3', 'chr4', 'chr5',
        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY']

cols_to_keep=["chrom", "start", "end", "svtype", "svlen", "filter", "sample_count", "sample_names", "genotype", "caller"]

def cutesv_extraction(merged_tmp_cutesv):
    merged_tmp = merged_tmp_cutesv[merged_tmp_cutesv.chrom.isin(chrs)]
    merged_tmp['filter'] = [','.join(map(str, l)) for l in merged_tmp['filter']]
    merged_tmp = merged_tmp[((merged_tmp['filter']=='PASS') | (merged_tmp["svtype"]=='INV') | (merged_tmp["svtype"]=='DUP')) & ~(merged_tmp["svtype"]=='UNK') & (merged_tmp["svlen"]>=50)]
    merged_tmp

    sample_hom_het = [col for col in merged_tmp.columns if "_GT" in col]
    list_tmp=merged_tmp[sample_hom_het].astype(str).replace(["(0, 0)", "(0, 1)", "(1, 1)", "(None, None)"],
                                                   ['hom_ref', 'het', 'hom_alt', "nan"]).values

    hom_het=[]
    for row in list_tmp:
        remove_nan_row = [x for x in row if x != "nan"]
        hom_het.append(remove_nan_row)


    sample_columns = [col for col in merged_tmp.columns if '_DV' in col]
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
         "sample_names":res,
         "genotype":hom_het
        },
        index=merged_tmp.index)
    tmp_list=[merged_tmp, tmp]
    merged_cutesv = pd.concat(tmp_list, axis=1)
    merged_cutesv['caller'] = 'cutesv'

    # merged_cutesv['sample_names'] = merged_cutesv['sample_names'].apply(convert_to_list)
    # merged_cutesv['genotype'] = merged_cutesv['genotype'].apply(convert_to_list)
    merged_cutesv['sample_name_genotypes'] = merged_cutesv.apply(lambda row: list(zip(row['sample_names'], row['genotype'])), axis=1)

    merged_cutesv_cropped=merged_cutesv[cols_to_keep]
    return merged_cutesv
def sniffles2_extraction(merged_tmp_sniffles2):
    merged_tmp = merged_tmp_sniffles2[merged_tmp_sniffles2.chrom.isin(chrs)]
    merged_tmp['filter'] = [','.join(map(str, l)) for l in merged_tmp['filter']]
    merged_tmp = merged_tmp[((merged_tmp['filter']=='PASS') | (merged_tmp["svtype"]=='INV') | (merged_tmp["svtype"]=='DUP')) & ~(merged_tmp["svtype"]=='UNK') & (merged_tmp["svlen"]>=50)]

    sample_columns = [col for col in merged_tmp.columns if '_DV' in col]
    sample_hom_het = [col for col in merged_tmp.columns if "_GT" in col]
    list_tmp=merged_tmp[sample_hom_het].astype(str).replace(["(0, 0)", "(0, 1)", "(1, 1)", "(None, None)"],
                                                   ['hom_ref', 'het', 'hom_alt', "nan"]).values
    hom_het=[]
    for row in list_tmp:
        remove_nan_row = [x for x in row if x != "nan"]
        hom_het.append(remove_nan_row)

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
         "sample_names":res,
         "genotype":hom_het
        },
        index=merged_tmp.index)

    tmp_list=[merged_tmp, tmp]
    merged_sniffles2 = pd.concat(tmp_list, axis=1)
    merged_sniffles2['caller'] = 'sniffles2'

    # merged_sniffles2['sample_names'] = merged_sniffles2['sample_names'].apply(convert_to_list)
    # merged_sniffles2['genotype'] = merged_sniffles2['genotype'].apply(convert_to_list)
    merged_sniffles2['sample_name_genotypes'] = merged_sniffles2.apply(lambda row: list(zip(row['sample_names'], row['genotype'])), axis=1)

    merged_sniffles2_cropped=merged_sniffles2[cols_to_keep]
    return merged_sniffles2


def sniffles2_extraction(merged_tmp_sniffles2):
    merged_tmp = merged_tmp_sniffles2[merged_tmp_sniffles2.chrom.isin(chrs)]
    merged_tmp['filter'] = [','.join(map(str, l)) for l in merged_tmp['filter']]
    merged_tmp = merged_tmp[((merged_tmp['filter']=='PASS') | (merged_tmp["svtype"]=='INV') | (merged_tmp["svtype"]=='DUP')) & ~(merged_tmp["svtype"]=='UNK') & (merged_tmp["svlen"]>=50)]

    sample_columns = [col for col in merged_tmp.columns if '_DV' in col]
    sample_hom_het = [col for col in merged_tmp.columns if "_GT" in col]
    list_tmp=merged_tmp[sample_hom_het].astype(str).replace(["(0, 0)", "(0, 1)", "(1, 1)", "(None, None)"],
                                                   ['hom_ref', 'het', 'hom_alt', "nan"]).values
    hom_het=[]
    for row in list_tmp:
        remove_nan_row = [x for x in row if x != "nan"]
        hom_het.append(remove_nan_row)

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
         "sample_names":res,
         "genotype":hom_het
        },
        index=merged_tmp.index)

    tmp_list=[merged_tmp, tmp]
    merged_sniffles2 = pd.concat(tmp_list, axis=1)
    merged_sniffles2['caller'] = 'sniffles2'

    # merged_sniffles2['sample_names'] = merged_sniffles2['sample_names'].apply(convert_to_list)
    # merged_sniffles2['genotype'] = merged_sniffles2['genotype'].apply(convert_to_list)
    merged_sniffles2['sample_name_genotypes'] = merged_sniffles2.apply(lambda row: list(zip(row['sample_names'], row['genotype'])), axis=1)

    merged_sniffles2_cropped=merged_sniffles2[cols_to_keep]
    return merged_sniffles2


def pbsv_extraction(merged_tmp_pbsv):
    merged_tmp = merged_tmp_pbsv[merged_tmp_pbsv.chrom.isin(chrs)]
    merged_tmp['filter'] = [','.join(map(str, l)) for l in merged_tmp['filter']]
    merged_tmp = merged_tmp[((merged_tmp['filter']=='PASS') | (merged_tmp["svtype"]=='INV') | (merged_tmp["svtype"]=='DUP')) & ~(merged_tmp["svtype"]=='UNK') & (merged_tmp["svlen"]>=50)]


    sample_columns = [col for col in merged_tmp.columns if "_AD_ref" in col]

    sample_hom_het = [col for col in merged_tmp.columns if "_GT" in col]
    list_tmp=merged_tmp[sample_hom_het].astype(str).replace(["(0, 0)", "(0, 1)", "(1, 1)", "(None, None)"],
                                                   ['hom_ref', 'het', 'hom_alt', "nan"]).values
    hom_het=[]
    for row in list_tmp:
        remove_nan_row = [x for x in row if x != "nan"]
        hom_het.append(remove_nan_row)

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
         "sample_names":res,
         "genotype":hom_het
        },
        index=merged_tmp.index)

    tmp_list=[merged_tmp, tmp]
    merged_pbsv = pd.concat(tmp_list, axis=1)
    merged_pbsv['caller'] = 'pbsv'

    # merged_pbsv['sample_names'] = merged_pbsv['sample_names'].apply(convert_to_list)
    # merged_pbsv['genotype'] = merged_pbsv['genotype'].apply(convert_to_list)
    merged_pbsv['sample_name_genotypes'] = merged_pbsv.apply(lambda row: list(zip(row['sample_names'], row['genotype'])), axis=1)

    merged_pbsv_cropped=merged_pbsv[cols_to_keep]
    return merged_pbsv


def svdss_extraction(merged_tmp_svdss):
    merged_tmp = merged_tmp[merged_tmp.chrom.isin(chrs)]
    merged_tmp['filter'] = [','.join(map(str, l)) for l in merged_tmp['filter']]
    merged_tmp = merged_tmp[((merged_tmp['filter']=='PASS') | (merged_tmp["svtype"]=='INV') | (merged_tmp["svtype"]=='DUP')) & ~(merged_tmp["svtype"]=='UNK') & (merged_tmp["svlen"]>=50)]
    merged_tmp

    sample_columns = [col for col in merged_tmp.columns if "_GT" in col]

    sample_hom_het = [col for col in merged_tmp.columns if "_GT" in col]
    list_tmp=merged_tmp[sample_hom_het].astype(str).replace(["(0, 0)", "(0, 1)", "(1, 1)", "(None, None)"],
                                                   ['hom_ref', 'het', 'hom_alt', "nan"]).values
    hom_het=[]
    for row in list_tmp:
        remove_nan_row = [x for x in row if x != "nan"]
        hom_het.append(remove_nan_row)

    c = merged_tmp[sample_columns].columns.to_numpy()
    c_tmp = []
    for DV in c:
        c_tmp.append(DV.split('.')[0])
    c_dict = dict(zip(c, c_tmp))

    caller=f"{DV.split('.')[2]}".split("_")[0]

    res = []
    sample_count=[]
    for x in merged_tmp[sample_columns].to_numpy():
        mask=['None' not in str(sample) for sample in x]
        res_temp1 = c[mask].tolist()
        res_temp2 = []
        for i in res_temp1:
            res_temp2.append(c_dict[i])
        names=list(set(res_temp2))
        res.append(names)
        sample_count.append(len(names))

    tmp = pd.DataFrame(
        {
         'sample_count': sample_count,
         "sample_names":res,
         'genotype':hom_het
        },
        index=merged_tmp.index)
    tmp_list=[merged_tmp, tmp]
    merged_svdss = pd.concat(tmp_list, axis=1)
    merged_svdss['caller'] = 'svdss'

    # merged_svdss['sample_names'] = merged_svdss['sample_names'].apply(convert_to_list)
    # merged_svdss['genotype'] = merged_svdss['genotype'].apply(convert_to_list)
    merged_svdss['sample_name_genotypes'] = merged_svdss.apply(lambda row: list(zip(row['sample_names'], row['genotype'])), axis=1)


    merged_svdss_cropped=merged_svdss[cols_to_keep]

    return merged_svdss