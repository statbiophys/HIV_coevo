import pandas as pd
import numpy as np
import sys

# Arg1: file_name of the input table
#print (sys.argv[1])
# Arg1: file_name of the output table
#print (sys.argv[1])


def select_first_V(row):
    row['v_call'] = row['v_call'].split(',')[0]
    return row
def select_first_J(row):
    row['j_call'] = row['j_call'].split(',')[0]
    return row
def set_junct_start_find(row):
    row['junction_start'] = row['sequence'].find(row['junction'])
    return row


f = pd.read_csv(sys.argv[1], sep='\t')

print('Filtering out bad clonotypes')
f = f[f['productive']=='T']
f = f[f['vj_in_frame']=='T']
f = f[f['locus']=='IGH']

print('Considering only the first V J match')
f = f.apply(select_first_V, axis=1)
f = f.apply(select_first_J, axis=1)

print('Identifying the beginning of the CDR3 seq')
f['junction_start'] = np.zeros(len(f))
f = f.apply(set_junct_start_find, axis=1)

f = f[['sequence_id', 'duplicate_count', 'sequence', 'sequence_alignment', 'v_call',
       'j_call', 'junction', 'junction_start', 'junction_length', 'v_sequence_start',
       'v_sequence_end', 'v_germline_start', 'v_germline_end', 'j_sequence_start',
       'j_sequence_end', 'j_germline_start', 'j_germline_end']]

f.to_csv(sys.argv[2],sep='\t')
