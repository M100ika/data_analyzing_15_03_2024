import numpy as np
import pandas as pd

def inmap_readaction(meta_map_init):
    meta_map_df = meta_map_init.copy()
    meta_map_df = meta_map_df.drop(columns=['SequencingRun'])
    for column in meta_map_df.select_dtypes(include=['object']):  
        meta_map_df[column] = meta_map_df[column].str.replace('.extendedFrags.fastq.gz', '', regex=False)

    meta_map_df['fastqFile'] = meta_map_df['fastqFile'].apply(lambda x: x.split('_')[0])
    meta_map_df = meta_map_df.rename(columns={'fastqFile':'ID'})
    meta_map_df['ID'] = pd.to_numeric(meta_map_df['ID'], errors='coerce')
    return meta_map_df


def write_uniq_data(meta_map_df, meta, flag = False):
    matching_values = meta_map_df['ID'].isin(meta['ID'])
    matched_values = meta_map_df['ID'][matching_values]

    unique_in_meta_map_df = meta_map_df[~meta_map_df['ID'].isin(meta['ID'])]['ID']

    unique_in_meta = meta[~meta['ID'].isin(meta_map_df['ID'])]['ID']

    if flag == True:
        file = 'data/output_data/uniq_data.csv'
        with open(file, 'w') as uniq_data_file:
            uniq_data_file.write(f"Совпадающие значения в обоих DataFrame ({len(matched_values.unique())}): {matched_values.unique()}\n\n"
                                f"Уникальные значения в in.map ({len(unique_in_meta_map_df.unique())}): {unique_in_meta_map_df.unique()}\n\n"
                                f"Уникальные значения в metadata ({len(unique_in_meta.unique())}): {unique_in_meta.unique()}\n\n"
                                f"Общее кол-во пропущенных данных: {len(unique_in_meta_map_df.unique()) + len(unique_in_meta.unique())}\n\n\n")
    else: 
        print(f'Write uniq data function. Вывод уникальный значений:')
        print((f"Совпадающие значения в обоих DataFrame ({len(matched_values.unique())}): {matched_values.unique()}\n\n"
                                f"Уникальные значения в in.map ({len(unique_in_meta_map_df.unique())}): {unique_in_meta_map_df.unique()}\n\n"
                                f"Уникальные значения в metadata ({len(unique_in_meta.unique())}): {unique_in_meta.unique()}\n\n"
                                f"Общее кол-во пропущенных данных: {len(unique_in_meta_map_df.unique()) + len(unique_in_meta.unique())}"))
        

def merge_meta_inmap(meta_map_df, meta):
    merged_df = pd.merge(meta_map_df, meta, on='ID')
    return merged_df[['#SampleID', 'ID','GROUP']]


def merge_genus_meta(final_merged_df, genus_df):
    genus_df_transposed = genus_df.T.reset_index()
    genus_df_transposed.columns = ['#SampleID'] + list(genus_df['Genus'])
    genus_df_transposed = genus_df_transposed.drop(genus_df_transposed.index[0])
    merged_genus_df = pd.merge(final_merged_df, genus_df_transposed, on='#SampleID', how='left')
    return merged_genus_df.drop(columns=['ID'])

def replace_group(merged_genus_df):
    replace_dict = {
        1: "AstNIS",
        2: "Astreg",
        3: "KokNis",
        4: "Kokreg",
        5: "SemNis",
        6: "SemReg",
        7: "UstNis",
        8: "Ustreg",
        9: "KarNis",
        10: "KarReg"
    }

    return merged_genus_df['GROUP'].replace(replace_dict)


def main():
    otu_df = pd.read_csv('/home/esp/data_analyze/15.03.2024/data/init_data/Res1/OTU.txt', sep='\t')
    hiera_blast = pd.read_csv('/home/esp/data_analyze/15.03.2024/data/init_data/Res1/hiera_BLAST.txt', sep='\t')
    meta_map_init = pd.read_csv('/home/esp/data_analyze/15.03.2024/data/init_data/Res1/primary/in.map', sep='\t')
    meta = pd.read_csv('/home/esp/data_analyze/15.03.2024/data/init_data/Res1/metadata_geo.csv',  sep='\t')
    genus_df = pd.read_csv('/home/esp/data_analyze/15.03.2024/data/init_data/Res1/higherLvl/Genus.txt', sep='\t')

    meta_map_df = inmap_readaction(meta_map_init)

    write_uniq_data(meta_map_df, meta, False)
    final_merged_df = merge_meta_inmap(meta_map_df, meta)
    merged_genus_df = merge_genus_meta(final_merged_df, genus_df)
    merged_genus_df = replace_group(merged_genus_df)
    

main()