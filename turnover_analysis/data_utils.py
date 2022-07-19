import numpy as np
import pandas as pd
import os
from Bio import SeqIO, AlignIO, Seq


# Functions for import-export of data


def read_times(path):
    times = []
    time_file = open(path)
    for l in time_file:
        times.append([int(t) for t in l.split('\t')[:-1]])
    time_file.close()
    return times


def read_fasta(path):
    """
    Read a fasta file and return it as a dictionary
    """
    fasta = SeqIO.parse(path, "fasta")
    seqs = dict()
    for record in fasta:
        seqs[record.id] = str(record.seq)
    return seqs


def write_fasta(path, seqs):
    """
    Write the seqs dictionary as a fasta file at the given path
    """
    f = open(path, 'w')
    for k, v in seqs.items():
        f.write('>'+k+'\n')
        f.write(v+'\n')
    f.close()
    
    
def import_hiv_fasta(directory, samples=[]):
    """
    Import all the hiv fasta files in the directory. Alternatively
    a list of samples to import can be specified
    """
    seqs = dict()
    for file_name in os.listdir(directory):
        split_name = file_name.split('_')
        if len(split_name) == 3 and split_name[2] == 'hiv.fasta':
            sample_id = split_name[0][3:]+'_'+split_name[1][3:]
            if len(samples) > 0 and sample_id not in samples:
                continue
            s = read_fasta(directory+"/"+file_name)
            seqs[sample_id] = s
    return seqs


def read_hiv_msa_fasta(path):
    """
    Read a msa list of hiv sequences saved in a fasta file.
    The fasta keys needs to be formatted as id_count_time.
    It returns a dataFrame containing this information
    """
    hiv_msa = read_fasta(path)
    framed_seqs = pd.DataFrame()
    framed_seqs['id'] = hiv_msa.keys()
    framed_seqs.index = np.arange(len(framed_seqs))
    aux_index = framed_seqs.id.str.split('_', expand=True)
    framed_seqs['counts'] = aux_index[1].astype(int)
    framed_seqs['time'] = aux_index[2].astype(int)
    framed_seqs['msa_seq'] = framed_seqs['id'].map(hiv_msa)
    framed_seqs.index = framed_seqs['id']
    return framed_seqs.drop(columns='id')


def import_fam_frames(directory, samples=[], keep_seq=False):
    """
    Import all the abr family dataframes in a directory and return the 
    frame having the family identifier as index (grouping together all the
    clonotypes of each family)
    """
    fams = dict()
    for file_name in os.listdir(directory):
        
        if file_name.split('_')[1][:3] == 'day': # Loop over the clone files
            pat = file_name.split('_')[0][3:]
            time = file_name.split('_')[1][3:]
            id_ = pat+'_'+time
            if len(samples) > 0 and id_ not in samples:
                continue
            fam_frame = pd.read_csv(directory+"/"+file_name,sep='\t')
            fam_frame['frequency'] = np.ones(len(fam_frame))/len(fam_frame)
            if keep_seq:
                fam_frame['seq_counts'] = fam_frame['counts'].copy()
                fams[pat+'_'+time] = fam_frame.groupby('family').agg({'counts':sum, 'frequency':sum, 'seq':list, 'seq_counts':list})
            else:
                fam_frame['seq_counts'] = fam_frame['counts'].copy()
                fams[pat+'_'+time] = fam_frame.groupby('family').agg({'counts':sum, 'frequency':sum, 'seq_counts':list})
            
    return fams


def import_fam_frame_pat(directory, pat, keep_seq=False):
    """
    Import all the abr family dataframes in a dir of a patient
    """
    samples = []
    for file_name in os.listdir(directory):
        if file_name.split('_')[1][:3] == 'day':
            p = int(file_name.split('_')[0][3:])
            if p == int(pat):
                samples.append(str(p)+'_'+file_name.split('_')[1][3:])
    return import_fam_frames(directory, samples, keep_seq)


def build_fam_timeframe(sing_fam_frames, groupbyfam=True):
    """
    Given the dictionary of abr family frames of a given patient, it builds
    a unique dataframe of families across time points. It also returns
    the list of times.
    The family time frame has the fields:
        - time_occ: number of time points in which the family is present
        - frequencies: list of frequencies of the seqs in the family at each 
            time point in which the family is present
        - counts: list of counts of the seqs in the family at each 
            time point in which the family is present
        - time: time indexes in which the family is present
        - seq_counts: list of the list of counts of all the sequences in a family
            across time points
        - seqs: list of lists of the nucl sequences in the family across time points
    """
    
    sort_sing_fams = []
    fam_times = []
    for id_, c in sing_fam_frames.items():
        p, time = id_.split('_')[0], id_.split('_')[1]
        sort_sing_fams.append(c)
        fam_times.append(int(time))
    sort_fam_times = np.argsort(fam_times)
    fam_times = np.array(fam_times)[sort_fam_times]
    sort_sing_fams = [sort_sing_fams[t] for t in sort_fam_times]
    
    all_fam_dict = pd.DataFrame(columns=['counts', 'frequency', 'seqs', 'seq_counts'])
    for i,s in enumerate(sort_sing_fams):
        sers = pd.DataFrame(s)
        sers.index = sers.index + ['_'+str(i)]*len(sers)
        all_fam_dict = all_fam_dict.append(sers)
        
    fam_frame = pd.DataFrame()
    fam_frame['id'] = all_fam_dict.index
    fam_frame.index = np.arange(len(fam_frame))
    aux_index = fam_frame.id.str.split('_', expand=True)
    fam_frame['counts'] = fam_frame['id'].map(all_fam_dict['counts']).astype(int)
    fam_frame['log+1_counts'] = np.log(fam_frame['counts'])+1
    fam_frame['frequencies'] = fam_frame['id'].map(all_fam_dict['frequency']).astype(float)
    fam_frame['len'] = aux_index[2].astype(int)
    if 'seqs' in all_fam_dict:
        fam_frame['seqs'] = fam_frame['id'].map(all_fam_dict['seq'])
    fam_frame['seq_counts'] = fam_frame['id'].map(all_fam_dict['seq_counts'])
    fam_frame['time'] = aux_index[4].astype(int)
    fam_frame['family'] = aux_index[0]+'_'+aux_index[1]+'_'+aux_index[2]+'_'+aux_index[3]
    fam_frame.index = fam_frame['id']
    
    if groupbyfam:
        fam_frame['time_occ'] = np.ones(len(fam_frame))
        fam_frame = fam_frame.groupby('family').agg({'time_occ':sum, 'frequencies':list, 'counts':list, 
                                            'log+1_counts':list, 'time':list, 'seq_counts':list, 'seqs':list})
    return fam_frame, fam_times


def find_adj_times(times):
    """
    Given a list of integers it finds the sequence2 of adjacent numbers (slow!).
    It is applied to the time field of row.
    It returns the list of indexes from which the adjacent integers start, and the
    list of numbers of consecutive numbers.
    """
    first_time_index, n_step = 0, 0
    first_time_indexes, n_steps = [], []
    for i in range(1,len(times)):
        if times[i]-times[i-1] != 1:
            if n_step != 0:
                first_time_indexes.append(first_time_index)
                n_steps.append(n_step+1)
                n_step = 0
            first_time_index = i
        else: n_step += 1
            
    if n_step != 0:  
        first_time_indexes.append(first_time_index)
        n_steps.append(n_step+1)
        
        
def find_adj_times(row):
    """
    Given a list of integers it finds the sequence2 of adjacent numbers (slow!).
    It is applied to the time field of row.
    It creates the list of indexes from which the adjacent integers start, and the
    list of numbers of consecutive numbers.
    """
    times = np.sort(row['time'])
    first_time_index, n_step = 0, 0
    first_time_indexes, n_steps = [], []
    for i in range(1,len(times)):
        if times[i]-times[i-1] != 1:
            if n_step != 0:
                first_time_indexes.append(first_time_index)
                n_steps.append(n_step+1)
                n_step = 0
            first_time_index = i
        else: n_step += 1
            
    if n_step != 0:  
        first_time_indexes.append(first_time_index)
        n_steps.append(n_step+1)
            
    row['first_time_index'] = first_time_indexes
    row['n_adj_times'] = n_steps
    row['n_adj_seqs'] = len(n_steps)
    return row


def _discard_times(times_to_discard, frame):
    max_t = max(frame['time'].values)
    for t in times_to_discard[::-1]:
        frame = frame[frame['time'] != t]
        new_t = t
        for t1 in range(t+1, max_t+1):
            if t1 in frame['time'].values:
                frame.loc[frame['time'] == t1,'time'] = new_t
                new_t += 1
    return frame

def align_times(hiv_seqs, fams, hiv_times, fam_times):
    """
    Discard the non overlapping times between abr and hiv
    """
    time_and = np.array(list(set(fam_times).intersection(set(hiv_times))), dtype=int)
    time_and = np.sort(time_and)
    time_diff_hiv = np.array(list(set(hiv_times).difference(set(fam_times))), dtype=int)
    bad_hiv_time_indexes = [np.where(hiv_times == t)[0][0] for t in time_diff_hiv]
    time_diff_fam = np.array(list(set(fam_times).difference(set(hiv_times))), dtype=int)
    bad_fam_time_indexes = [np.where(fam_times == t)[0][0] for t in time_diff_fam]
    hiv_seqs = _discard_times(bad_hiv_time_indexes, hiv_seqs)
    fams = _discard_times(bad_fam_time_indexes, fams)
    return hiv_seqs, fams, time_and