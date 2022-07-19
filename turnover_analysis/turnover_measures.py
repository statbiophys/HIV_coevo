import numpy as np
import turnover_corr_utils as ut


class hiv_info_container:
    def __init__(self, msa_frame, time_snp_frame, syn_prob_rand, times):
        self.msa_frame = msa_frame
        self.time_snp_frame = time_snp_frame
        self.times = times
        self.nT = len(times)
        self.syn_prob_rand = syn_prob_rand
        
class abr_info_container:
    def __init__(self, time_frame, sub_time_frame, snp_dict, n_seq, syn_prob_rand, times):
        self.time_frame = time_frame
        self.sub_frame = sub_time_frame
        self.snp_dict = snp_dict
        self.n_seq = n_seq
        self.times = times
        self.nT = len(times)
        self.syn_prob_rand = syn_prob_rand
        
        
# GENERIC MEASURES

class delta_times():
    def run(self, info):
        return info.times[1:]-info.times[:-1]
    def label(self):
        return 'Time point difference'
    def short_label(self):
        return 'delta_time'
    
    
# HIV TURNOVER MEASURES 

class hiv_nseq():
    def run(self, info):
        n_seqs = [len(info.msa_frame[info.msa_frame['time']==t]) for t in range(max(info.msa_frame['time']+1))]
        return (np.array(n_seqs)[1:]+np.array(n_seqs)[:-1])/2
    def label(self):
        return 'Hiv unique sequence number'
    def short_label(self):
        return 'hiv_nseq'

class hiv_turnover_freq():
    def run(self, info):
        return ut.get_snp_abs_turnover(info.time_snp_frame, info.nT)
    def label(self):
        return 'Hiv snp turnover frequency'
    def short_label(self):
        return 'hiv_turn_freq'
    
class hiv_turnover_synratio_freq():
    def run(self, info):
        consensus_t0 = ut.build_consensus_sequence(info.msa_frame[info.msa_frame['time'] == 0]['msa_seq'].values)
        syn_prob = ut.compute_syn_prob_in_seq(consensus_t0, info.syn_prob_rand)
        frame_syn = info.time_snp_frame[info.time_snp_frame['syn']==True]
        frame_nsyn = info.time_snp_frame[info.time_snp_frame['syn']==False]
        turn_freq_syn = ut.get_snp_abs_turnover(frame_syn, info.nT)
        turn_freq_nsyn = ut.get_snp_abs_turnover(frame_nsyn, info.nT)
        return  turn_freq_nsyn / turn_freq_syn * syn_prob / (1-syn_prob)
    def label(self):
        return 'Hiv snp turnover nsyn-syn ratio frequency'
    def short_label(self):
        return 'hiv_turn_synratio_freq'
    
class hiv_turnover_count():
    def run(self, info):
        return ut.get_snp_abs_turnover(info.time_snp_frame, info.nT, 'counts')
    def label(self):
        return 'Hiv snp turnover count'
    def short_label(self):
        return 'hiv_turn_count'
    
class hiv_turnover_synratio_count():
    def run(self, info):
        consensus_t0 = ut.build_consensus_sequence(info.msa_frame[info.msa_frame['time'] == 0]['msa_seq'].values)
        syn_prob = ut.compute_syn_prob_in_seq(consensus_t0, info.syn_prob_rand)
        frame_syn = info.time_snp_frame[info.time_snp_frame['syn']==True]
        frame_nsyn = info.time_snp_frame[info.time_snp_frame['syn']==False]
        turn_freq_syn = ut.get_snp_abs_turnover(frame_syn, info.nT, 'counts')
        turn_freq_nsyn = ut.get_snp_abs_turnover(frame_nsyn, info.nT, 'counts')
        return  turn_freq_nsyn / turn_freq_syn * syn_prob / (1-syn_prob)
    def label(self):
        return 'Hiv snp turnover nsyn-syn ratio count'
    def short_label(self):
        return 'hiv_turn_synratio_count'
    

# ABR TURNOVER MEASURES
    
class abr_nseq():
    def run(self, info):
        return (np.array(info.n_seq)[1:]+np.array(info.n_seq)[:-1])/2
    def label(self):
        return 'Abr unique sequences number'
    def short_label(self):
        return 'abr_nseq'
    
class abr_lin_turnover_freq():
    def run(self, info):
        return ut.get_snp_abs_turnover(info.time_frame, info.nT)
    def label(self):
        return 'Abr lineage turnover frequency'
    def short_label(self):
        return 'abr_lin_turn_freq'

class abr_lin_turnover_count():
    def run(self, info):
        return ut.get_snp_abs_turnover(info.time_frame, info.nT, 'counts')
    def label(self):
        return 'Abr lineage turnover count'
    def short_label(self):
        return 'abr_lin_turn_count'
    
class abr_turnover_count_large():
    def run(self, info):
        freq_trajs = ut._get_snp_trajectories(info.time_frame, info.nT, 'counts')
        turnovers = np.abs(freq_trajs[:,1:] - freq_trajs[:,:-1])
        large_turn_average = []
        for t in range(info.nT-1):
            vals = turnovers[:,t]
            vals = vals[vals != 0]
            sort_vals = np.sort(vals)
            large_turn_average.append(sort_vals[int(0.90*len(vals)):].sum())
        return large_turn_average
    def label(self):
        return 'Abr lineage turnover >90 count'
    def short_label(self):
        return 'abr_lin_turn_large'
    
class abr_turnover_count_small():
    def run(self, info):
        freq_trajs = ut._get_snp_trajectories(info.time_frame, info.nT, 'counts')
        turnovers = np.abs(freq_trajs[:,1:] - freq_trajs[:,:-1])
        large_turn_average = []
        for t in range(info.nT-1):
            vals = turnovers[:,t]
            vals = vals[vals != 0]
            sort_vals = np.sort(vals)
            large_turn_average.append(sort_vals[:int(0.50*len(vals))].sum())
        return large_turn_average
    def label(self):
        return 'Abr lineage turnover <50 count'
    def short_label(self):
        return 'abr_lin_turn_small'
    
class abr_snp_turnover_freq():
    def run(self, info):
        return ut.compute_abr_snp_turnover(info.snp_dict, info.sub_frame, info.nT)
    def label(self):
        return 'Abr snp turnover frequency'
    def short_label(self):
        return 'abr_turn_freq'
    
class abr_snp_turnover_count():
    def run(self, info):
        return ut.compute_abr_snp_turnover(info.snp_dict, info.sub_frame, info.nT, 'counts')
    def label(self):
        return 'Abr snp turnover count'
    def short_label(self):
        return 'abr_turn_count'
        

    