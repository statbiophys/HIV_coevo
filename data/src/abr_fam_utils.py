import numpy as np
import pandas as pd
import atriegc
from matplotlib import cm
import matplotlib.pyplot as plt
from itertools import combinations



class abr_families_gen(object):
    
    def __init__(self, frames, count_labels, seq_labels, Vgene_labels, Jgene_labels, ignoreJ=False, verbose=True):
        """
        Class that builds a list the abr-family dataframes from a 
        list of clonotype dataframes. The frames need to contain a count filed, a nucelotide 
        seqeunce field, a Vgene field and a Jgene field. 
        The corresponding labels of those fields for each frame need to be specified as arguments. 
        The abr families are build using single linkage clustering (using hamming distance
        between sequencies) where the threshold (over the hamming distance rescaled by the 
        sequence length, in [0,1]) is specified as argument.
        The clustering is performed within each class of clonotypes defined by a given V, 
        J genes and a given length (VJL class).
        """
        
        # Arguments check
        self._args_check(frames, count_labels, seq_labels, Vgene_labels, Jgene_labels)
        self.ignoreJ = ignoreJ
        self.verbose=verbose
        
        # Frame pre-processing
        self.frames = []
        for i, frame in enumerate(frames):
            f = frame[[self.count_l[i], self.seq_l[i], self.Vgene_l[i], self.Jgene_l[i]]].copy()
            self.frames.append(self._preprocess_frame(f, i))
    
        # Merging all the dataframes
        self.merged_frame = pd.DataFrame()
        for frame in self.frames:
            if len(self.merged_frame) == 0:
                self.merged_frame = frame.copy()
            else:
                self.merged_frame = self._merge_frames(self.merged_frame, frame.copy())
                
        if self.ignoreJ:
            self.merged_frame = self.merged_frame.groupby('VL-seq').agg({
                'counts':sum, 'Vgene':'first', 'seq':'first',
                'VL':'first', 'lengths': 'first', 'VL-seq':'first'
            })
        else:
            self.merged_frame = self.merged_frame.groupby('VJL-seq').agg({
                'counts':sum, 'Vgene':'first', 'Jgene':'first', 'seq':'first',
                'VJL':'first', 'lengths': 'first', 'VJL-seq':'first'
            })
                
            
    def run(self, threshold, min_len=20, compute_translator=False):
        """
        It returns the list of dataframes with counts for each abr family. The indexes of those
        families are in the form: "Vgene_Jgene_sequence-length_cluster-id".
        The cluster ids are integers that can be translated into the list of nucleotide 
        sequences using class variable cluster_translator.
        """
        # Short sequences are filtered away
        self.min_len = min_len
        filt_frame = self.merged_frame[self.merged_frame['lengths'] > min_len]
        
        # Computing the single linkage clusters of the sequences in each VJL class
        self.abr_fam_clusters = self._compute_clusters(filt_frame, threshold)
        if compute_translator:
            self.cluster_translator = self._compute_cluster_translator()
        
        # Building the count dataframe where the indexes are the abr families
        if self.verbose: print('Identifying clusters in samples...')
        fams = []
        for i in range(len(self.frames)):
            if self.verbose: print(str(i+1)+'/'+str(len(self.frames)))
            fams.append(self._make_families_frame(i))
            
        return fams

    
    def plot_seq_dist_per_length(self, ax, n_largest_len_classes, max_class_size=1000):
        ax.set_yscale('log')
        ax.set_xlabel('Hamm dist / CDR3 length', fontsize=14)
        ax.set_ylabel('pdf', fontsize=14)
        
        dists_per_len = self._compute_dist_per_len(n_largest_len_classes, max_class_size)
        ls = {int(l) for l in dists_per_len.keys()}
        for l, dists in dists_per_len.items():
            l = int(l)
            bins=(np.arange(0,max(dists)*l+1)+0.5)/l
            h=plt.hist(dists, bins=bins, density=True, color='b', alpha=0)

            l_norm = (l-min(ls))/(max(ls)-min(ls))
            im = ax.plot(h[1][1:]-1/l/2, h[0], c=cm.viridis(l_norm), lw=2)

        sm = cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=min(ls), vmax=max(ls)))
        sm._A = []
        cbar = plt.colorbar(sm)
        cbar.set_label('CDR3 length', fontsize=14)
        return ax

        
    def _args_check(self, frames, count_labels, seq_labels, Vgene_labels, Jgene_labels):
        
        self.count_l = count_labels
        self.seq_l = seq_labels
        self.Vgene_l = Vgene_labels
        self.Jgene_l = Jgene_labels
        
        if type(count_labels) == str: 
            self.count_l = [count_labels for _ in range(len(frames))]
        if type(seq_labels) == str: 
            self.seq_l = [seq_labels for _ in range(len(frames))]
        if type(Vgene_labels) == str: 
            self.Vgene_l = [Vgene_labels for _ in range(len(frames))]
        if type(Jgene_labels) == str: 
            self.Jgene_l = [Jgene_labels for _ in range(len(frames))]

        if len(frames) != len(self.count_l) or len(frames) != len(self.seq_l) or \
           len(frames) != len(self.Vgene_l) or len(frames) != len(self.Jgene_l):
            raise Exception('Inconsistent lengths of the arguments')
            

    def _preprocess_frame(self, frame, i):
        """
        Adding auxiliary columns to the dataframes containing length, VJL class and 
        unique identifier.
        Columns are renamed in the default way.
        """
        #f = frame.copy()
        f = frame
        # Length of the nucleotide sequence
        f['lengths'] = [len(ncl) for ncl in f[self.seq_l[i]]]
        if not self.ignoreJ:
            # Vgene-Jgene-length- class of the abr
            f['VJL'] = f[self.Vgene_l[i]] + '_' + f[self.Jgene_l[i]] + '_' + f['lengths'].astype(str)
            # Vgene-Jgene-length-sequence of the abr (uniquely identify the clone)
            f['VJL-seq'] = f['VJL'] + '_' + f[self.seq_l[i]]
            f.rename(columns={self.count_l[i]:'counts', self.seq_l[i]:'seq', 
                              self.Vgene_l[i]:'Vgene', self.Jgene_l[i]:'Jgene'}, inplace=True)
            f = f.groupby('VJL-seq').agg({
                'counts':sum,'Vgene':'first','Jgene':'first','seq':'first',
                'VJL':'first','lengths':'first', 'VJL-seq':'first'
            })
        else:
            # Vgene-Jgene-length- class of the abr
            f['VL'] = f[self.Vgene_l[i]] + '_' + f['lengths'].astype(str)
            # Vgene-Jgene-length-sequence of the abr (uniquely identify the clone)
            f['VL-seq'] = f['VL'] + '_' + f[self.seq_l[i]]
            f.rename(columns={self.count_l[i]:'counts', self.seq_l[i]:'seq', 
                              self.Vgene_l[i]:'Vgene'}, inplace=True)
            f = f.groupby('VL-seq').agg({
                'counts':sum,'Vgene':'first','seq':'first',
                'VL':'first','lengths':'first', 'VL-seq':'first'
            })
        f = f.set_index(np.arange(len(f)))
        return f


    def _merge_frames(self, f1, f2):
        """
        Merging two abr dataframes on the unique identifier VJL-seq.
        """
        lab = 'VJL'
        if self.ignoreJ: lab = 'VL'
        merged_frame = pd.merge(f1, f2, on=lab+'-seq', suffixes=('', '1'), how='outer')
        merged_frame['counts'].fillna(int(0),inplace=True)
        merged_frame['counts1'].fillna(int(0),inplace=True)
        merged_frame['counts'] = merged_frame['counts']+merged_frame['counts1']
        merged_frame['seq'].fillna(merged_frame['seq1'], inplace=True)
        merged_frame['Vgene'].fillna(merged_frame['Vgene1'], inplace=True)
        if not self.ignoreJ:
            merged_frame['Jgene'].fillna(merged_frame['Jgene1'], inplace=True)
        merged_frame[lab].fillna(merged_frame[lab+'1'], inplace=True)
        merged_frame['lengths'].fillna(merged_frame['lengths1'], inplace=True)
        if not self.ignoreJ:
            merged_frame = merged_frame[['counts','Vgene','Jgene','seq',lab,'lengths',lab+'-seq']]
        else:
            merged_frame = merged_frame[['counts','Vgene','seq',lab,'lengths',lab+'-seq']]
        return merged_frame


    def _compute_class_clusters(self, class_seqs, threshold):
        """
        Computing the abr-family single-linkage cluster of a given VJL class.
        It returns a dictionary where the sequence (key) is associated with the
        identifier of its cluster (an integer value).
        """
        trie = atriegc.TrieNucl()
        trie.insert_list(class_seqs)
        l_threshold = int(len(class_seqs[0])*threshold)
        t_cluster = trie.clusters(l_threshold)
        return t_cluster


    def _compute_clusters(self, frame, threshold):
        """
        Computing the abr-family single-linkage cluster of a sample.
        It returns a dictionary of a dictionary. The first key is the VJL class,
        the associated dictionary has the sequence (key) associated with the
        identifier of its cluster (an integer value).
        """
        if self.verbose: print('Computing clusters...')
        lab = 'VJL'
        if self.ignoreJ: lab = 'VL'
            
        clusters = dict()
        VJLs_uniq = set(frame[lab])
        for VJL in VJLs_uniq:
            sub_frame = frame[frame[lab] == VJL]
            clusters[VJL] = self._compute_class_clusters(sub_frame['seq'].values, threshold)
        return clusters


    def _compute_cluster_translator(self):
        """
        It returns a dictionary of a dictionary. The first key is the VJL class,
        the associated dictionary has cluster identifier (key) associated with 
        list of sequences belonging to the cluster.
        """
        if self.verbose: print('Computing cluster translator...')
            
        clust_transl = dict()
        for k, d in self.abr_fam_clusters.items():
            clust_transl[k] = pd.DataFrame({'seq' : list(d.keys()), 'cluster' : list(d.values()),
                                            'n_seq' : np.ones(len(d))})
            clust_transl[k] = clust_transl[k].groupby('cluster').agg({'seq':lambda x: list(x),
                                                                      'n_seq':np.sum})
        return clust_transl


    def _make_families_frame(self, i):
        
        #frame = frame_clone.copy()
        #frame = frame[frame['lengths'] > self.min_len]
        #frame['occ'] = np.ones(len(frame))
        lab = 'VJL'
        if self.ignoreJ: lab = 'VL'
        for id_, row in self.frames[i].iterrows():
            if row['lengths'] > self.min_len:
                fams_dict = self.abr_fam_clusters[row[lab]]
                cluster = int(fams_dict[row['seq']])
                self.frames[i].loc[id_, 'family'] = row[lab] + '_' + str(cluster)
        #frame['family'] = frame['VJL'] + '_' + frame['cluster'].astype(int).astype(str)
        #frame = frame.groupby('VJL-cluster').agg({'counts':sum, 'occ':sum})
        #frame['family'] = frame.index
        return self.frames[i].drop([lab+'-seq', 'lengths', lab], axis=1)
    
    
    def _compute_dist_per_len(self, n_largest_len_classes, max_class_size=1000):
        lens = {int(l) for l in set(self.merged_frame.lengths)}
        ls, count_ls = np.unique(self.merged_frame.lengths, return_counts=True)
        dists_len = dict()
        lens_list = ls[np.argsort(count_ls)][-n_largest_len_classes:]
        for l in lens_list:
            sub_frame = self.merged_frame[self.merged_frame['lengths'].astype(int)==int(l)]
            if len(sub_frame)<=1:
                continue
            max_n_clones = min(len(sub_frame),max_class_size)
            dists_len[l] = self._compute_h_dists(list(sub_frame['seq'][:max_n_clones]))
        return dists_len
            
            
    def _hamm_dist(self, seq1, seq2):
        l1, l2 = np.array(list(seq1)), np.array(list(seq2))
        return np.sum(l1 != l2)

    
    def _compute_h_dists(self, seqs):
        dists = []
        l = len(seqs[0])
        for pair in combinations(seqs, 2):
            dists.append(self._hamm_dist(pair[0], pair[1])/l)
        return np.array(dists)