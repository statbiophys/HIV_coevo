#include "../wright_fisher.h"


// Population of binary spins that coevolve

// Compilation string:
// g++ -o spin50.exe spin50.cpp ../wright_fisher.cpp ../utils.cpp -lgsl -lgslcblas -std=c++17

// Execution string:
// ./spin50.exe [directory where the parameter file is saved]

// Parameter file must be coomposed of lines like:
// prameter_id \t parameter_value



int main(int argc, char *argv[]){

    if (seq_len != 50)
        throw std::runtime_error("Set sequence length to 50 in wright_fisher.h");

    // Importing the parameters
    Param params = import_param_from_main(argc, argv);
    bool print_info = params.get_double("print_info");
    int n_realizations = params.get_double("n_realizations");
    int n_gen = params.get_double("n_generations");
    int n_init_seq = params.get_double("n_init_seq");
    int pop_size1 = params.get_double("pop_size1");
    int pop_size2 = params.get_double("pop_size2");
    doublev turn_long_dts = params.get_doublev("turn_long_dts");
    int traj_step = params.get_double("traj_step");
    double mut_rate1 = params.get_double("mut_rate1");
    double mut_rate2 = params.get_double("mut_rate2");
    double sel_coef1 = params.get_double("sel_coef1");
    double sel_coef2 = params.get_double("sel_coef2");

    // Creating and seeding the random number generator
    gsl_rng* rng = get_rand_num_generator(params.get_double("seed"));

    doublev2 av_fitness(2, doublev(n_gen / (float)traj_step));
    std::vector<doublev2> av_turnover_long(turn_long_dts.size());
    for (int i=0; i<av_turnover_long.size(); i++) {
        av_turnover_long[i] = doublev2(2, doublev(n_gen / (float)turn_long_dts[i]-1));
    }

    Timer timer = Timer();
    Perc perc(10, n_realizations);
    // Iteration over the realizations
    for (int r=0; r<n_realizations; r++){

        // Building the populations
        int seq1[seq_len], seq2[seq_len];
        gntp gnts1[max_unique_gnt], gnts2[max_unique_gnt];
        for(int g=0; g<n_init_seq; g++){
            generate_rand_seq(rng, seq1);
            gnts1[g] = new Genotype_Lin(g, g, pop_size1/n_init_seq, 0, seq1, mut_rate1, sel_coef1);
            generate_rand_seq(rng, seq2);
            gnts2[g] = new Genotype_Lin(g, g, pop_size2/n_init_seq, 0, seq1, mut_rate2, sel_coef2);
        }

        // Initializing the populations
        Population pop1 = Population(0, 1, gnts1, n_init_seq, false);
        Population pop2 = Population(1, 0, gnts2, n_init_seq, false);
        poppv populations = poppv { &pop1, &pop2 };

        // Initializing and running the WF
        Wright_Fisher wf = Wright_Fisher(rng, print_info);
        wf.run(n_gen, populations, params.get_double("traj_step"));
        if (r == 0)
            wf.print_trajectory(argv[1]);

        // Computing the observables
        doublev2 traj;
        wf.compute_av_fitness(traj);
        for(int i=0; i<n_gen/(float)traj_step; i++){
            av_fitness[0][i] += traj[0][i] / (float)n_realizations;
            av_fitness[1][i] += traj[1][i] / (float)n_realizations;
        }
        
        for (int k=0; k<av_turnover_long.size(); k++){
            doublev2 turn = doublev2();
            wf.compute_turnover(turn, turn_long_dts[k]);
            for(int i=0; i<av_turnover_long[k][0].size(); i++){
                av_turnover_long[k][0][i] += turn[0][i] / (float)n_realizations; 
                av_turnover_long[k][1][i] += turn[1][i] / (float)n_realizations; 
            }
        }
        
        perc.step(r); 
    }

    // Printing the observables
    print_2dtraj(av_fitness, (str)argv[1]+"/av_fitness.txt");
    for (int k=0; k<av_turnover_long.size(); k++){
        print_2dtraj(av_turnover_long[k], (str)argv[1]+"/turnover_"+std::to_string((int)turn_long_dts[k])+".txt");
    }

    std::cout << "elapsed time: " << timer.elapsed() << "\n";

    return 0;
}