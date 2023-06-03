import sys
import re
import pandas as pd


INDEX_FILTER = 6

INDEX_INFO = 7
  
INFO_FIELDS = { 'FractionInformativeReads': r'FractionInformativeReads=(\d+\.\d+)',
                              'MQ': r'MQ=(\d+\.\d+)',
                              'MQRankSum': r'MQRankSum=(\d+\.\d+)', 
                              'ReadPosRankSum': r'ReadPosRankSum=(\d+\.\d+)',
                              'hotspot': r'\b(hotspot)\b', 
                              'GermlineStatus': r'.*;(?:GermlineStatus=(Germline_DB|Germline_Proxi|Somatic|Somatic_Putative_DB))?'}

FILTER_FIELDS = ['base_quality', 'filtered_reads', 'fragment_length', 
                                'low_depth', 'low_frac_info_reads', 'low_normal_depth', 
                                'long_indel', 'mapping_quality', 'non_homref_normal', 
                                'no_reliable_supporting_read', 'panel_of_normals',
                                'read_position', 'RMxNRepeatRegion', 'str_contraction',
                                'too_few_supporting_reads', 'Somatic_Quality', 'weak_evidence']


def parse_filters(filters):
  '''Parse the VCF filters field.
  VCF header: #  CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	custom_somatico
  '''
  filters_status = {}
  for filter in FILTER_FIELDS:
    filters_status[filter] = str(True) if filter in filters else str(False)
  return filters_status


def parse_infos(infos):
  '''Parse the VCF info field.
  VCF header: #  CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	custom_somatico
  '''
  info_fields = {}
  info_fields['MQ'] = re.search(INFO_FIELDS['MQ'], infos).group(1) 
  info_fields['MQRankSum'] = re.search(INFO_FIELDS['MQRankSum'], infos).group(1) if re.search(INFO_FIELDS['MQRankSum'], infos) else None
  info_fields['FractionInformativeReads'] = re.search(INFO_FIELDS['FractionInformativeReads'], infos).group(1) 
  info_fields['ReadPosRankSum'] = re.search(INFO_FIELDS['ReadPosRankSum'], infos).group(1) if re.search(INFO_FIELDS['ReadPosRankSum'], infos) else None
  info_fields['hotspot'] =re.search(INFO_FIELDS['hotspot'], infos).group(1) if re.search(INFO_FIELDS['hotspot'], infos) else None
  info_fields['GermlineStatus'] = re.search(INFO_FIELDS['GermlineStatus'], infos).group(1) 
  return info_fields


def print_as_table(filename, data, filter_fields, info_fields):
  '''Print each data in a dataframe and export it to csv.
  '''
  pd.set_option('display.max_columns', None)
  pd.set_option('display.max_rows', None)
  pd.set_option('display.width', 200)
  header = [*['CHROM', 'POS'], *filter_fields, *list(info_fields.keys())]
  df = pd.DataFrame(data, columns=header)
  df.to_csv(vcf_file + '.table.csv', index=False)
  print(df)


if __name__ == '__main__':
    vcf_file = sys.argv[1]
    with open(vcf_file, 'r') as f:
    	lines = f.readlines()
    	data = []
    	for line in lines:
    		if '#' not in line:
    			cols = line.split('\t')
    			chrom = cols[0]
    			pos = cols[1]
    			filter_fields = parse_filters(cols[INDEX_FILTER].split(';'))
    			info_fields = parse_infos(cols[INDEX_INFO])
    			data.append([*[chrom, pos], *filter_fields.values(), *info_fields.values()] )
    	print_as_table(vcf_file, data, filter_fields, info_fields)
    