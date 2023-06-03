import sys
import re
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 200)
  
VCF_COLUMNS = ['CHROM', 'POS',	'ID',	'REF',	'ALT', 	'QUAL', 	'FILTER', 	'INFO', 'FORMAT', 'custom_somatico']

INFO_FIELDS_PATTERNS = { 'FractionInformativeReads': r'FractionInformativeReads=(\d+\.?\d+)',
                                                  'MQ': r'MQ=(\d+\.?\d+)',
                                                  'MQRankSum': r'MQRankSum=(\d+\.?\d+)', 
                                                  'ReadPosRankSum': r'ReadPosRankSum=(\d+\.?\d+)', 
                                                  'hotspot': r'\b(hotspot)\b', 
                                                  'GermlineStatus': r'GermlineStatus=(Germline_DB|Germline_Proxi|Somatic|Somatic_Putative_DB)'}


def parse_filters(df):
    '''Parse the VCF FILTER field.'''
    filters = []
    df_filter = df['FILTER'].str.split(';', expand=True)
    for col in df_filter.columns:
        filters.extend(list(df_filter[col].value_counts().index))
    for filter in set(filters):
        df.loc[df['FILTER'].str.contains(filter, na=False), filter] = True


def parse_infos(df):
    '''Parse the VCF INFO field.'''
    df_infos = df['INFO'].str.split(';',expand=True)
    for col in df_infos.columns:
        col = df_infos[col].value_counts().index
        for item in col:
            if '=' in item:
                item = item.split('=')[0]
            print(item)
            if  item in INFO_FIELDS_PATTERNS.keys():
                df[item] = df['INFO'].str.extract(INFO_FIELDS_PATTERNS[item])


if __name__ == '__main__':
    vcf_file = sys.argv[1]
    df = pd.read_csv(vcf_file, sep='\t', comment='#', names=VCF_COLUMNS)
    df = df[['CHROM', 'POS', 'FILTER', 'INFO']]
    parse_filters(df)
    parse_infos(df)
    df.fillna('', inplace=True)
    df.drop(['FILTER', 'INFO'], axis=1, inplace=True)
    df.to_csv('teste.csv', index=False)
    print(df)
