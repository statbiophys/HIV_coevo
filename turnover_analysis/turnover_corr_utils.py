import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO, Seq
from itertools import product



def build_consensus_sequence(seqs, weights=[]):
    """
    Builds the consensus sequence among seqs (of same length) of nucleotides+gap-symbol. 
    The sequences can be weighted. The consensus has the symbols with larger weight.
    """
    
    letters = ['A', 'T', 'C', 'G', '-']
    letter_index = {letters[i]:i for i in range(len(letters))}

    if len(weights) == 0:
        weights = np.ones(len(seqs[0]))
    
    count_letter_mat = np.zeros((len(seqs[0]), 5), dtype=float)
    for s, w in zip(seqs, weights):
        count_indexes = np.array([letter_index[l] for l in s])
        count_letter_mat[np.arange(len(s)), count_indexes] += w

    consensus = ''
    for c in np.argmax(count_letter_mat, axis=1):
        consensus+=letters[c]
    return consensus



# BUILDING SNP DATAFRAMES

def build_snp_frame(seqs, reference_seq, seq_ids=[], counts=[]):
    """
    It builds the Single Nucleotide Polimorfism composition of a list of sequences, seqs,
    with respect to a reference sequence. Meaning that only the mutations that are different
    from the reference are tracked.
    Sequences must be framed and multiple aligned, gaps "-" are treated as letters.
    It returns a dataFrame with indexes as [mutation poistion]_[letter], and columns:
        - position: the bp position from the beginning of the seq on the mutation
        - letter: nucleotide of the mutation
        - frequency: fraction of sequences that contain that snp
        - syn: whether the mutation is synonimous with respect the reference seq
        - counts: abundance of that mutation, computed if counts arg is a list of ints which
            specifies the abundance of each sequence
        - seq_id: list of ids of the sequences containing that mutation, computed only if
            seq_ids arg is a list containing the id for each sequence.
    """
    
    # The sequences are converted in a matrix of letters
    seqs = np.array([list(s) for s in seqs])
    N = len(seqs[:,0])
    
    # Reference in aminoacid
    aa_reference = Seq.translate(reference_seq, gap='-')
    
    snp_frame = pd.DataFrame(columns=['id', 'position', 'letter', 'frequency', 'syn', 'counts', 'seq_ids'])
    
    # Iteration over all the positions of the sequences.
    for i, site_composisiton in enumerate(seqs.T):
        
        # Finding reference codon
        codon_index = int(i/3)
        reference_codon = reference_seq[codon_index*3:codon_index*3+3]
        # Locating the mutations
        is_mutation = site_composisiton != reference_seq[i]
        mutations = site_composisiton[is_mutation]
        mut_letter, mut_occ = np.unique(mutations, return_counts=True)
        
        # Iteration over the mutations
        for l, c in zip(mut_letter, mut_occ):
            
            if aa_reference[codon_index] == '-' or l == '-':
                is_syn = False # A gap is always considered as a non-syn
            else:
                mut_codon = list(reference_codon)
                mut_codon[i%3] = l
                mut_aa = Seq.translate(''.join(mut_codon))
                is_syn = mut_aa==aa_reference[codon_index]
            
            snp_dict = {'id': str(i)+'_'+str(l), 'position': i, 'letter': l, \
                        'syn': is_syn, 'frequency': c/N}
            mut_indexes = np.arange(N)[site_composisiton == l]
            if len(counts) == N:
                snp_dict['counts'] = sum([counts[i] for i in mut_indexes])
            if len(seq_ids) == N:
                snp_dict['seq_ids'] = [seq_ids[i] for i in mut_indexes]
            snp_frame = snp_frame.append(snp_dict, ignore_index=True)
    
    if len(counts) != N: snp_frame = snp_frame.drop(columns={'counts'})
    if len(seq_ids) != N: snp_frame = snp_frame.drop(columns={'seq_ids'}) 
        
    return snp_frame.set_index('id')


def build_time_snp_frame(seq_time_frame, consensus):
    """
    Given the msa_frame with columns: time, counts, msa_seq, it builds a new dataframe containing
    the information about SNPs in time.
    The returned frame has the same fields position, letter, syn as the build_snp_frame method,
    but now there is a time field which is a ordered list of time indexes in which the mutation 
    is non-zero. The frequencies and counts fields are list containing frequenceise and counts 
    at the times specified in the time filed (in a sort of sparse representation).
    """
    
    time_snp_frame = pd.DataFrame(columns=['id', 'position', 'letter', 'syn', 'frequency', 'counts'])
    T = max(seq_time_frame['time'])

    for t in range(T+1):
        time_frame = seq_time_frame[seq_time_frame['time'] == t]
        snp_frame = build_snp_frame(time_frame['msa_seq'].values, consensus, counts=time_frame.counts.values)
        snp_frame['time'] = np.ones(len(snp_frame))*t
        snp_frame['id'] = snp_frame.index
        snp_frame = snp_frame.set_index(np.arange(len(snp_frame)))
        time_snp_frame = time_snp_frame.merge(snp_frame, how='outer')

    time_snp_frame = time_snp_frame.groupby('id').agg({'position':'first', 'letter':'first', 'syn': 'first', 
                                                       'frequency':list, 'counts':list, 'time':list})
    return time_snp_frame.sort_values('position').rename(columns={'frequency':'frequencies'})


def _build_abr_time_snp_frame_singlefam(row, T):
    """
    Same procedure of build_time_snp_frame but for a single family row of the abr_time_frame.
    """
    time_snp_frame = pd.DataFrame(columns=['id', 'position', 'letter', 'frequency', 'counts'])
    seqs_dict = dict()
    consensus = build_consensus_sequence(row['seqs'][row['first_time_index'][0]])
    for it_start, nt in zip(row['first_time_index'], row['n_adj_times']):
        for i in range(it_start, it_start+nt):
            t = row['time'][i]
            snp_frame = build_snp_frame(row['seqs'][i], consensus, counts=row['seq_counts'][i])
            snp_frame['time'] = np.ones(len(snp_frame))*t
            snp_frame['id'] = snp_frame.index
            snp_frame = snp_frame.set_index(np.arange(len(snp_frame)))
            time_snp_frame = time_snp_frame.merge(snp_frame, how='outer')

    time_snp_frame = time_snp_frame.groupby('id').agg({'position':'first', 'letter':'first', 'syn':'first',
                                                       'frequency':list, 'counts':list, 'time':list})
    return time_snp_frame


def build_abr_time_snp_frame(fam_time_frame, T):
    """
    Same procedure of build_time_snp_frame but for an abr_time_frame.
    """
    snps_dict = dict()
    for family, row in fam_time_frame.iterrows():
        snps_dict[family] = _build_abr_time_snp_frame_singlefam(row, T)
    return  snps_dict



# COMPUTING TURNOVER OF SNPS OR FAMILIES OF AB


def _get_snp_trajectories(time_snp_frame, nT, label='frequencies'):
    """
    Given the time snp frame, it returns the trajectories of the 
    frequencies in an array form
    """
    times = time_snp_frame.time.values
    freqs = time_snp_frame[label].values
    freq_trajs = np.zeros((len(time_snp_frame), nT))
    for i in range(len(times)):
        np.put(freq_trajs[i], times[i], freqs[i])
    return freq_trajs


def get_snp_abs_turnover(time_snp_frame, nT, label='frequencies'):
    """
    Get the absolute turnover measure of a given time snp frame
    """
    freq_trajs = _get_snp_trajectories(time_snp_frame, nT, label)
    return np.sum(np.abs(freq_trajs[:,1:] - freq_trajs[:,:-1]),axis=0)


def compute_abr_snp_turnover_singlesnp(snp_traj, fam_info, is_abs=True):
    """
    Compute the turnover of a snps of a family.
    """
    if is_abs:
        tr = np.abs(snp_traj[1:]-snp_traj[:-1])
    else:
        tr = snp_traj[1:]-snp_traj[:-1]
    # Removing the extremes of the trajectory, where the family was not present
    #for i_start, n in zip(fam_info.first_time_index, fam_info.n_adj_times):
    #    snps_times = fam_info.time[i_start:i_start+n]
    #    if snps_times[0]>0:
    #        tr[snps_times[0]-1]=0
    #    if snps_times[-1] < len(tr):
    #        tr[snps_times[-1]]=0
    return tr


def compute_abr_snp_turnover_singlefam(snps, fam_info, T, traj_label='frequency', is_abs=True):
    """
    Compute the turnover of a snps of a family
    """
    turnovers = np.zeros(T-1)
    for index, row in snps.iterrows():
        f = np.zeros(T)
        np.put(f, row['time'], row[traj_label])
        turnovers += compute_abr_snp_turnover_singlesnp(f, fam_info, is_abs)
    return np.array(turnovers)


def compute_abr_snp_turnover(snps_dict, sub_fam_time_frame, T, traj_label='frequency', is_abs=True):
    turnovers = np.zeros(T-1)
    for family, snps in snps_dict.items():
        fam_info = sub_fam_time_frame.loc[family]
        turnovers += compute_abr_snp_turnover_singlefam(snps, fam_info, T, traj_label, is_abs)
    return np.array(turnovers)


def compute_abr_snp_turnover_synratio(snps_dict, sub_fam_time_frame, T, traj_label='frequency', is_abs=True):
    turnover_syn, turnover_nsyn = np.zeros(T-1), np.zeros(T-1)
    syn_prob_rand_codon = build_syn_prob_rand()
    for family, snps in snps_dict.items():
        fam_info = sub_fam_time_frame.loc[family]
        consensus = build_consensus_sequence(fam_info['seqs'][0])
        syn_prob = compute_syn_prob_in_seq(consensus, syn_prob_rand_codon)
        c = syn_prob/(1-syn_prob)
        turnover_syn += compute_abr_snp_turnover_singlefam(snps[snps['syn']==True], fam_info, T, traj_label, is_abs)
        turnover_nsyn += compute_abr_snp_turnover_singlefam(snps[snps['syn']==False], fam_info, T, traj_label, is_abs)*c
    return np.array(turnover_nsyn)/np.array(turnover_syn)



def build_syn_prob_rand():
    """
    It builds the dictionary which associates at each existing codon the probability
    of observing a synonimous mutation arising in a complete random scenario
    """
    nucleotides = ['A','C','T','G']
    expect_syn = dict()

    for l1,l2,l3 in product(nucleotides,nucleotides,nucleotides):
        codon = ''.join([l1,l2,l3])
        aa = Seq.translate(codon)
        n_syn = 0
        for lett in nucleotides:
            for i in range(3):
                if lett != codon[i]:
                    codon_l = list(codon)
                    codon_l[i] = lett
                    n_syn += Seq.translate(''.join(codon_l)) == aa
        expect_syn[codon] = n_syn/9.0
        expect_syn['---'] = 0
    return expect_syn


def compute_syn_prob_in_seq(seq, expect_syn):
    """
    It computes the probability of synonimous mutations that can appear in
    the sequence given the probability of appearence at each codon specified
    in expect_syn
    """
    fract_syn = 0
    for i in range(int(len(seq)/3.0)):
        codon = seq[3*i:3*(i+1)]
        fract_syn += expect_syn[codon]
    return fract_syn/float(int(len(seq)/3.0))


