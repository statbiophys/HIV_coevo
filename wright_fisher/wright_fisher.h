#ifndef WRIGHT_FISHER_H
#define WRIGHT_FISHER_H


#include "utils.h"


const bool debug=false;
const int seq_len = 1;
const int max_unique_gnt = 1000;

class Population;
class Genotype;


using i_i_d = std::tuple<int, int, double>;
using wf_pop_traj = std::vector<std::vector<i_i_d>>;
using wf_traj = std::vector<wf_pop_traj>;
using wf_info = std::vector<std::map<int, str>>;
using gntp = Genotype*;
using popp = Population*;
using poppv = std::vector<popp>;


void generate_rand_seq(gsl_rng* rng, int* seq);


/* Abstract Genotype defined by a binary string which interact with a population
   and receives a fitness poportional to a function of some average fitness weight
   (like activation probabilities) to specify in the derived class.
   The constant of proportinality is the ftn_factor.
   Mutations flip a loci. It is assumed that maximum one mutation can appear
   in a string (mut_rate*seq_len << 1). */
class Genotype {

    protected:

        // Index within the population vectors
        int index;
        // Id of the parent genotype. Equal to same id if it is a founder.
        long parent; 
        // Unique identifier within the population
        long id;
        // Number of individuals
        int abundance;
        // Average exponential fitness
        double av_exp_fitness;
        // Mutation rate per genotype (per sequence)
        double mut_rate;
        // Constant of proportinality of the final fitness
        double ftn_factor;
        // List of fitness weights toward each genotype of the interacting population
        // The final fitness is ftn_function of this average times ftn_factor.
        // fw_d ftn_weights;
        double ftn_weights[max_unique_gnt];
        // Effective length of the ftn_weight array
        int ftn_weights_size;

        // Overlap with another strain
        int affinity(int other_sequence[]);

        // Explicit computation of fitness weight given affinity of the two seqs
        virtual double ftn_weight(int affinity){
            std::runtime_error("ftn_weight to implement");
            return -1;
        }
        // Function of the average fitness weights
        virtual double ftn_function(double average_ftn_weight) { 
            return average_ftn_weight; 
        }
        virtual double affinity_from_weight(double fw) {
            std::runtime_error("affinity_from_weight to implement");
            return -1;
        }
        virtual Genotype* construct_new_mutant(int index, int id, int seq[]){
            std::runtime_error("new_mutant to implement");
            return NULL;
        }
 
    public:

        // Index of the interacting population
        int my_inter_pop;
        // Sequence characterizing the genotype
        int sequence[seq_len];

        Genotype();
        Genotype(
            int index, long id, int abundance, long parent, int seq[], 
            double mut_rate, double ftn_fact
        );

        long get_id() const { return id; }
        void set_id(long new_id) { id = new_id; }
        int get_index() const { return index; }
        void set_index(int index) { this->index = index; }
        int get_abundance() const { return abundance; }
        void set_abundance(int ab) { abundance = ab; }
        double get_av_exp_fitness() const { return av_exp_fitness; }
        void set_av_exp_fitness(const Population& interacting_pop);
        double get_mut_rate() const { return mut_rate; }
        bool operator==(const gntp g); 
        str get_info();
        gntp new_mutant(int index, long id, const poppv& pop, gsl_rng* rng);

        // Methods to update the fitness weights
        void remove_adv_gnt_weight(int index);
        void add_adv_gnt_weight(const gntp gnt);
        // Direct computation of weights across all adversary population 
        void init_ftn_weights(const popp adv_pop);
        // Computation of weights by modifying the ones of the parent genotype
        void init_ftn_weights(double parent_fw[], int mut_site, const popp pop);
};


/* Genotype defined by a binary string which interact with a population
   and receives a fitness poportional to the log of the average activation
   probability (1+exp[-beta(affinity/seq_len - affinity_th)])^-1 
   over that population. */
class Genotype_Sigm : public Genotype {

    protected:
        // Affinity th normalized by the sequence length, in [-1,1]
        double affinity_th;
        double beta;

        double ftn_function(double average_ftn_weight);
        double ftn_weight(int affinity);
        double affinity_from_weight(double fw);
        Genotype* construct_new_mutant(int index, int id, int seq[]);

    public:
        Genotype_Sigm(
            int index, long id, int ab, long parent, int seq[], double mut_rate, double ftn_fact,
            double aff_th, double beta):
            Genotype{index, id, ab, parent, seq, mut_rate, ftn_fact}, affinity_th{aff_th}, beta{beta} {};
};


/* Genotype defined by a binary string which interact with a population
   and receives a fitness poportional to the log of the average activation
   probability (1+exp[-beta(affinity/seq_len - affinity_th)])^-1 
   over that population. */
class Genotype_Lin : public Genotype {

    protected:

        double ftn_weight(int affinity);
        double affinity_from_weight(double fw);
        Genotype* construct_new_mutant(int index, int id, int seq[]);

    public:
        Genotype_Lin(int index, long id, int ab, long parent, int seq[], double mut_rate, double ftn_fact):
            Genotype{index, id, ab, parent, seq, mut_rate, ftn_fact} {};
        double set_exp_fitness(const poppv& populations);
};





/* It contains information of a set of evolving genotypes */
class Population {

    protected:

        int pop_id;
        long max_gnt_id;
        int pop_size;
        int n_unique_gnt;
        gntp genotypes[max_unique_gnt];

        void remove_adv_pop_gnt(int gnt_index);
        void add_adv_pop_gnt(const gntp gnt);

    public:
        int interacting_pop_id; 
        // Whether the genotypes with same sequences are collapsed together. It slows performance
        bool collapse_unique;

        Population(int pop_id, int inter_pop_id, gntp genotypes[], int n_unique_gnt, bool collapse_unique=false);
        ~Population() { 
            for (int i=0; i<n_gnt(); i++){ 
                delete genotypes[i]; 
                genotypes[i] = NULL;
            }
        }

        long gnt_id(int gnt_index) const { return (*genotypes[gnt_index]).get_id(); }
        int gnt_ab(int gnt_index) const { return (*genotypes[gnt_index]).get_abundance(); }
        double gnt_exp_ftn(int gnt_index) const { return (*genotypes[gnt_index]).get_av_exp_fitness(); }
        str gnt_info(int gnt_index) const {return (*genotypes[gnt_index]).get_info(); }
        gntp get_gnt(int gnt_index) { return genotypes[gnt_index]; } 
        gntp get_gnt(int gnt_index) const { return genotypes[gnt_index]; } 
        int size() const { return pop_size; }
        int n_gnt() const { return n_unique_gnt; }
        int id() const {return pop_id; }

        void init_ftn_weights(const poppv& populations);
        void set_gnt_ab(int gnt_index, int ab, poppv& populations);
        void mutate(wf_info& gnt_info, poppv& populations, gsl_rng* rng);
        void set_av_exp_fitness(const poppv& populations);  
};


/* Wright Fisher model for coevolving populations */
class Wright_Fisher {
    
    private:

        /* Random number generator */
        gsl_rng* rng;
        /* Abundance trajectory for every population at each time step.
           The type uint_uint_d gives the abundance as value for every 
           genotype id as key. */
        wf_traj traj;
        /* Information of each created genotype durning the process 
           for each population */
        wf_info gnt_info;
        int n_generations, traj_step;
        bool print_info;

    public:
        Wright_Fisher(gsl_rng* rng, bool print_info) : rng(rng), print_info{print_info} {}
        /* Run the WF model for a given number of t_steps. The init_pop sets which populations
           are interacting and from which starting abundance. traj_step is the number of steps
           between each save of the abundance trajectory. */
        void run(int t_steps, poppv& init_pop, int traj_step);
        /* Print the abundance trajectories for each population in a given directory */
        void print_trajectory(str out_dir, str name="") const;
        /* It compute the average fitness trajectory */
        void compute_av_fitness(doublev2& av_fitness, str out_path="");
        /* It compute the absolute turnover trajectory at delta_t steps. n_points says how many 
           time points of the turnover are generated considering the last part of the trajectory.
           n_points=0 for computing all the possible points. */
        void compute_turnover(doublev2& turnover, int delta_t=1, int n_points=0, str out_path="");
};


#endif