#include "wright_fisher.h"



Genotype::Genotype(
    int index, long id, int abundance, long parent, int seq[], 
    double mut_rate, double ftn_fact
):
    index{index}, id{id}, abundance{abundance}, mut_rate{mut_rate}, av_exp_fitness{1}, 
    parent{parent}, ftn_factor{ftn_fact}, ftn_weights_size{0} {
    for (int i=0; i<seq_len; i++) sequence[i] = seq[i];
}


Genotype::Genotype():
    index{0}, id{0}, abundance{0}, mut_rate{0}, av_exp_fitness{1}, parent{0}, ftn_factor{0}, ftn_weights_size{0} { }


int Genotype::affinity(int other_seq[]){
    int affinity = 0;
    for (int i=0; i<seq_len; i++) affinity += other_seq[i]*sequence[i];
    return affinity;
}


void Genotype::init_ftn_weights(const popp adv_pop) {
    for(int i=0; i<(*adv_pop).n_gnt(); i++) {
        ftn_weights[i] = ftn_weight( affinity(adv_pop->get_gnt(i)->sequence) );
        ftn_weights_size++;
        // for(int j=0; j<seq_len; j++)
        //     std::cout << sequence[j] << (adv_pop->get_gnt(i)->sequence)[j] << " ";
        // std::cout << affinity(adv_pop->get_gnt(i)->sequence) << " " << ftn_weights[i] << "\t";
    }
    // std::cout << "\n";
}


void Genotype::init_ftn_weights(double parent_fw[], int mut_site, const popp adv_pop){
    int my_site = sequence[mut_site];
    for (int i=0; i<(*adv_pop).n_gnt(); i++){
        Genotype& gnt = *(*adv_pop).get_gnt(i);
        int other_site = gnt.sequence[mut_site];
        double aff = affinity_from_weight(parent_fw[i]);
        aff += 2*other_site*my_site;
        ftn_weights[i] = ftn_weight(aff);
        ftn_weights_size++;
    }
}


void Genotype::remove_adv_gnt_weight(int index) {
    for(int i=index; i<ftn_weights_size-1; i++)
        ftn_weights[i] = ftn_weights[i+1];     
    ftn_weights_size--;
}


void Genotype::add_adv_gnt_weight(const gntp gnt) {
    if (ftn_weights_size != (*gnt).index)
        throw std::runtime_error("Inconsitent gnt index of " + std::to_string(id));
    ftn_weights[ftn_weights_size] = ftn_weight(affinity(gnt->sequence));
    ftn_weights_size++;
}


void Genotype::set_av_exp_fitness(const Population& pop) {
    // To avoid computing affinities in fake populations
    if (ftn_factor==0) return;
    
    double av_fw = 0;
    for (int i=0; i<pop.n_gnt(); i++){
        Genotype& gnt = *pop.get_gnt(i);
        av_fw += ftn_weights[i] * gnt.get_abundance();
    }
    av_exp_fitness = exp( ftn_function(av_fw / (double)pop.size()) * ftn_factor );
}


gntp Genotype::new_mutant(int index, long id, const poppv& pops, gsl_rng* rng) {
    int new_seq[seq_len];
    for (int i=0; i<seq_len; i++) new_seq[i] = sequence[i];
    unsigned long int mut_loci = gsl_rng_uniform_int(rng, seq_len);
    new_seq[mut_loci] *= -1;
    Genotype* gnt = construct_new_mutant(index, id, new_seq);
    gnt->init_ftn_weights(ftn_weights, mut_loci, pops[my_inter_pop]);
    return gnt;
}


bool Genotype::operator==(const gntp g) {
    for (int l=0; l<seq_len; l++)
        if (sequence[l] != g->sequence[l]) 
            return false;
    return true;
}


str Genotype::get_info() { 
    str seq_str = std::to_string(parent)+"_";
    for (int i=0; i<seq_len; i++){
        if (sequence[i]==1) seq_str += "+";
        if (sequence[i]==-1) seq_str += "-";
    }
    return seq_str; 
}


double Genotype_Sigm::ftn_weight(int affinity) {
    return 1/( 1+exp( -beta*( affinity/(float)seq_len - affinity_th ) ) );
}


double Genotype_Sigm::affinity_from_weight(double fw){
    return seq_len*(affinity_th - log(1/fw - 1)/beta);
}


double Genotype_Sigm::ftn_function(double average_ftn_weight) {
    return log(average_ftn_weight);
}


Genotype* Genotype_Sigm::construct_new_mutant(int index, int new_id, int seq[]){
    Genotype* gnt = new Genotype_Sigm(
        index, new_id, 1, get_id(), seq, mut_rate, ftn_factor, affinity_th, beta
    );
    gnt->my_inter_pop = my_inter_pop;
    return gnt;
}


double Genotype_Lin::ftn_weight(int affinity) {
    return affinity/2.0;
}


double Genotype_Lin::affinity_from_weight(double fw){
    return fw * 2;
}


Genotype* Genotype_Lin::construct_new_mutant(int index, int new_id, int seq[]){
    Genotype* gnt = new Genotype_Lin(
        index, new_id, 1, get_id(), seq, mut_rate, ftn_factor
    );
    (*gnt).my_inter_pop = my_inter_pop;
    return gnt;
}


void generate_rand_seq(gsl_rng* rng, int* seq){
    for (int i=0; i<seq_len; i++)
        seq[i] = gsl_rng_uniform_int(rng, 2)*2 - 1;
}


Population::Population(int id, int inter_pop_id, gntp genotypes[], int n_unique_gnt, bool collapse_unique) : 
pop_id{id}, interacting_pop_id{inter_pop_id}, n_unique_gnt{n_unique_gnt}, collapse_unique{collapse_unique} {
    // Setting the population size and the maximum index
    pop_size = 0;
    max_gnt_id = 0;
    std::set<int> ids = std::set<int>();
    for (int i=0; i<n_unique_gnt; i++){
        Genotype& gnt = (*genotypes[i]);
        this->genotypes[i] = genotypes[i];
        ids.insert(gnt.get_id());
        pop_size += gnt.get_abundance();
        if (gnt.get_id() > max_gnt_id) max_gnt_id = gnt.get_id();
        gnt.my_inter_pop = inter_pop_id;
    }
    // Check on repeated genotype ids
    if (ids.size() != n_unique_gnt)
        throw std::runtime_error("Repeated genotype ids in population "+std::to_string(id));
} 


void Population::set_gnt_ab(int gnt_index, int ab, poppv& populations) {
    // If 0 abundance the genotype is removed from the population
    if (ab <= 0){
        if (debug) std::cout << "Erased " << (*genotypes[gnt_index]).get_id() << "\n";
        delete genotypes[gnt_index];
        genotypes[gnt_index] = NULL;

        for(int i=gnt_index; i<n_gnt()-1; i++){
            (*genotypes[i+1]).set_index(i);
            genotypes[i] = genotypes[i+1];            
        }
        n_unique_gnt--;
        // The adversary population is updated
        (*populations[interacting_pop_id]).remove_adv_pop_gnt(gnt_index);
    }
    else
        // Normally set abundance if poistive
        (*genotypes[gnt_index]).set_abundance(ab);
}


void Population::init_ftn_weights(const poppv& pops) {
    for (int i=0; i<n_unique_gnt; i++)
        genotypes[i]->init_ftn_weights(pops[interacting_pop_id]);
}


void Population::remove_adv_pop_gnt(int gnt_index){
    for(int i=0; i<n_gnt(); i++)
        genotypes[i]->remove_adv_gnt_weight(gnt_index);
}


void Population::add_adv_pop_gnt(const gntp gnt){
    for(int i=0; i<n_gnt(); i++)
        genotypes[i]->add_adv_gnt_weight(gnt);
}


void Population::set_av_exp_fitness(const poppv& populations){
    for (int i=0; i<n_unique_gnt; i++) {
        (*genotypes[i]).set_av_exp_fitness(*populations[interacting_pop_id]);
    }
}


void Population::mutate(wf_info& gnt_info, poppv& populations, gsl_rng* rng){

    // Iteration over all the genotypes before mutations
    int old_n_gnt = n_unique_gnt;
    for (int i=0; i<old_n_gnt; i++){

        Genotype& gnt = (*genotypes[i]);

        // Creating a binom distributed number of mutations per genotype
        int n_mut = gsl_ran_binomial(rng, gnt.get_mut_rate(), gnt.get_abundance());
        for (int j=0; j<n_mut; j++){
            // Creation of the new mutant
            gntp new_gen = gnt.new_mutant(n_unique_gnt, max_gnt_id+1, populations, rng);

            int g_found_i = -1;
            if (collapse_unique) {
                // Check among the existing genotypes
                for (int k=0; k<n_gnt(); k++) {
                    if ((*genotypes[k]) == new_gen) { 
                        g_found_i = k;
                        break;
                    }
                }                
                if (g_found_i >= 0) {
                    //delete new_gen;
                    //new_gen = NULL;
                    Genotype& equal_gnt = *genotypes[g_found_i];
                    equal_gnt.set_abundance(equal_gnt.get_abundance() + 1);
                    if (debug) 
                        std::cout << "New mut from " << gnt.get_id() << " collapsed in " << equal_gnt.get_id() << "\n";
                } 
                // Check among the extinct genotypes and rename the id if found
                str seq = new_gen->get_info();
                seq = seq.substr(seq.find('_')+1);
                std::map<int, str>::iterator it;
                for (it=gnt_info[id()].begin(); it!=gnt_info[id()].end(); it++){
                    str old_seq = it->second;
                    old_seq = old_seq.substr(old_seq.find('_')+1);
                    if (seq == old_seq){
                        new_gen->set_id(it->first);
                        break;
                    }
                }

                if (g_found_i >= 0) {
                    delete new_gen;
                    new_gen = NULL;
                } 
            }

            // Mutant added to the population
            if (g_found_i < 0) {
                genotypes[n_unique_gnt] = new_gen;
                max_gnt_id++;
                n_unique_gnt++;
                if (n_unique_gnt >= max_unique_gnt)
                    throw std::runtime_error("Max number of unique genotypes reached");
                if (debug) std::cout << "New mut " << max_gnt_id << " from " << gnt.get_id() << "\n";

                (*populations[interacting_pop_id]).add_adv_pop_gnt(new_gen);
                
                if (gnt_info.size() > 0)
                    gnt_info[id()][(*new_gen).get_id()] = (*new_gen).get_info();
            }
        }

        if (n_mut > 0) {
            // Updating mother gnt abundance
            int new_ab = gnt.get_abundance() - n_mut;
            set_gnt_ab(i, new_ab, populations);
            if (new_ab <= 0){ 
                i--; // all the other indexes have been decr by one because of the extinction
                old_n_gnt--; 
            }
        }
    }
}


void Wright_Fisher::run(int t_steps, poppv& populations, int traj_step){

    if (debug) std::cout << "Realization starts\n";

    for (popp pop : populations) {
        if(pop->collapse_unique && !print_info){
            print_info = true;
            std::cout << "Print info set true, required for collapse unique\n";
        }
    }

    int t_traj = 0;
    this->traj_step = traj_step;
    int n_traj_steps = t_steps/(float)traj_step+1;
    traj = wf_traj( populations.size(), wf_pop_traj(n_traj_steps) );

    gnt_info = std::vector<std::map<int, str>>();
    if (print_info){
        gnt_info = std::vector<std::map<int, str>>(populations.size());
        for(popp pop : populations)
            for  (int i=0; i<(*pop).n_gnt(); i++)
                gnt_info[(*pop).id()].insert(std::make_pair((*pop).gnt_id(i), (*pop).gnt_info(i)));
    }

    for(popp pop_ptr : populations)
        pop_ptr->init_ftn_weights(populations);

    n_generations = t_steps;
    
    for (int t=0; t<t_steps; t++){

        if (debug) std::cout << "Time " << t << "\n";

        for(popp pop_ptr : populations){
            Population& pop = (*pop_ptr);

            // Update the fitness for every genotype and computing the multinomial weights
            doublev samp_weights = doublev(pop.n_gnt());
            pop.set_av_exp_fitness(populations);
            for (int i=0; i<pop.n_gnt(); i++)
                samp_weights[i] = pop.gnt_exp_ftn(i) * pop.gnt_ab(i);

            if (debug) {
                for (int i=0; i<pop.n_gnt(); i++)
                    std::cout << pop.gnt_id(i) << ":" << pop.gnt_exp_ftn(i) << "\t";
                std::cout << "\n";
            }

            // Multinomial sampling with GSL libraries
            std::vector<unsigned int> new_abundances = std::vector<unsigned int>(pop.n_gnt());
            gsl_ran_multinomial(rng, pop.n_gnt(), pop.size(), &samp_weights[0], &new_abundances[0]);
            for (int i=pop.n_gnt()-1; i>=0; i--) pop.set_gnt_ab(i, new_abundances[i], populations);

            if (debug) {
                for (int i=0; i<pop.n_gnt(); i++)
                    std::cout << pop.gnt_id(i) << ":" << pop.gnt_ab(i) << "\t";
                std::cout << "\n";
            }

            // Mutation step
            pop.mutate(gnt_info, populations, rng);
                
            if (debug) std::cout << "\n";
        }

        // Saving the trajectory
        if (t%traj_step == 0){
            for (int p=0; p<populations.size(); p++){
                const Population& pop = (*populations[p]);
                std::vector<i_i_d> info_list(pop.n_gnt());
                for (int i=0; i<pop.n_gnt(); i++){
                    i_i_d info {pop.gnt_id(i), pop.gnt_ab(i), pop.gnt_exp_ftn(i)};
                    info_list[i] = info;
                }
                traj[p][t_traj] = info_list;
            }
            t_traj++;
        }
    }
}


void Wright_Fisher::print_trajectory(str out_dir, str name) const {
    if (name != "")
        name = name+"_";
    for (int p=0; p<traj.size(); p++){
        std::ofstream out;
        out.open(out_dir+"/trajectory_"+name+std::to_string(p)+".txt");
        for (const std::vector<i_i_d>& info_list : traj[p]){
            for (const i_i_d& it: info_list)
                out << std::get<0>(it) << ":" << std::get<1>(it) << "," << std::get<2>(it) << "\t";
            out << "\n";
        }
        out.close();

        if (print_info) {
            out.open(out_dir+"/gnt_info_"+name+std::to_string(p)+".txt");
            std::map<int, str>::iterator it;
            std::map<int, str> m = gnt_info[p];
            for (it = m.begin(); it != m.end(); it++)
                out << it->first << ":" << it->second << "\n";
            out.close();
        }
    }
}


void Wright_Fisher::compute_av_fitness(doublev2& av_fitness, str out_path) {
    av_fitness = doublev2(traj.size(), doublev(traj[0].size()));
    for (int p=0; p<traj.size(); p++){
        for(int t=0; t<traj[p].size(); t++){
            int n = 0;
            for (const i_i_d& it: traj[p][t]) {
                av_fitness[p][t] += std::get<1>(it) * log(std::get<2>(it));
                n += std::get<1>(it);
            }
            av_fitness[p][t] /= (float)n;
        }
    }

    if (out_path != "")
        print_2dtraj(av_fitness, out_path);
}


void Wright_Fisher::compute_turnover(doublev2& turnover, int delta_t, int n_points, str out_path) {
    if (delta_t % traj_step != 0)
        throw std::runtime_error("Turnover dt incompatible with trajectory step");
    if ((n_points+1)*delta_t+1 > traj_step*traj[0].size())
        throw std::runtime_error("Turnover points incompatible with n generations");

    int traj_dt = delta_t/(float)traj_step;

    if (n_points == 0)
        n_points = traj[0].size()/(float)traj_dt - 1; 

    turnover = doublev2(2, doublev(n_points));
    //std::cout << n_points << " "  << turnover.size() << " "  << traj[0].size() << "\n";

    for (int p=0; p<traj.size(); p++) {
        //std::cout << "pop " << p << " turn " << delta_t << "\n";
        int t_start = traj[p].size() - (n_points+1)*traj_dt + traj_dt;
        dict1000<int,int> old_ids = dict1000<int,int>();
        for (int t=t_start; t<traj[p].size(); t+=traj_dt) {
            int i = (t-t_start)/(float)traj_dt;
            //std::cout << i << " t1 " << t  << " t0 " << t-traj_dt << "\n";
            int turn = 0;
            old_ids.clear();
            for (const i_i_d& it: traj[p][t-traj_dt])
                old_ids.add(std::get<0>(it), std::get<1>(it));
            for (const i_i_d& it: traj[p][t]) {
                int id = std::get<0>(it);
                // If the genotype is not among the ones at previous step, 
                // the turnover contribution is its abundance.
                if (!old_ids.is_present(id)){
                    turn += std::get<1>(it);
                    //std::cout << "new " << id << ":" << std::get<1>(it) << "\t";
                }
                // Turnover as the abs difference between abundances
                else{
                    int delta = std::get<1>(it) - old_ids[id];
                    turn += std::abs(delta);
                    //std::cout << " old " << id << ":" << std::abs(delta) << "\t";
                    old_ids.remove(id);
                }
            }
            // For all the not-found old genotypes the turnover is the old abundance
            for(int j=0; j<old_ids.size(); j++){
                turn += old_ids.get_value(j);
                //std::cout  << " extinct " << old_ids.get_key(j) << ":" << old_ids.get_value(j) << "\t";
            }
            
            turnover[p][i] = turn;
            //std::cout << i <<  "turn " << turn << "\n";
        }
    }

    if (out_path != "")
        print_2dtraj(turnover, out_path);
}