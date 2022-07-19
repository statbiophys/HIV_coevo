#include "../wright_fisher.h"


// Population of binary spins that coevolve

// Compilation string:
// g++ -o single_spin_scan.exe single_spin_scan.cpp ../wright_fisher.cpp ../utils.cpp -lgsl -lgslcblas -std=c++17

// Execution string:
// ./single_spin_scan.exe [directory where the parameter file is saved]

// Parameter file must be coomposed of lines like:
// prameter_id \t parameter_value



int main(int argc, char *argv[]){

    if (seq_len != 1)
        throw std::runtime_error("Set sequence length to 1 in wright_fisher.h");

    // Creating and seeding the random number generator
    gsl_rng* rng = get_rand_num_generator();

    // Importing the parameters
    Param params = import_param_from_main(argc, argv);
    int n_realizations = params.get_double("n_realizations");
    int relax_time_base = params.get_double("relax_time_base");
    int relax_time_rand = params.get_double("relax_time_rand");
    int n_turn_points = params.get_double("n_turn_points");
    int pop_size1 = params.get_double("pop_size1");
    int pop_size2 = params.get_double("pop_size2");
    int traj_step = params.get_double("traj_step");
    int traj_len = params.get_double("traj_len");
    int traj_step2 = params.get_double("traj_step2");
    double mut_rate1 = params.get_double("mut_rate1");
    double mut_rate2 = params.get_double("mut_rate2");
    doublev fitn_coefs1 = params.get_doublev("sel_coefs1");
    doublev fitn_coefs2 = params.get_doublev("sel_coefs2");
    doublev turn_long_dts = params.get_doublev("turn_long_dt");

    Timer timer = Timer();
    // Iteration over the fitness coefficients of the two populations
    for (int i=0; i<fitn_coefs1.size(); i++) {

        // Only the turnovers are stored and printed.
        // For each value of dt considered, the turnovers of all the realizations
        // are appended and printed in a single file.
        std::vector<doublev2> turnovers( 
            turn_long_dts.size(), doublev2( 2, doublev(n_turn_points*n_realizations) ) 
        );

        // First run for the trajectory
        Wright_Fisher wf = Wright_Fisher(rng, true);
        int seq[1] {1};
        gntp gnts1[max_unique_gnt];
        int start_ab1 = gsl_rng_uniform(rng) * pop_size1;
        gnts1[0] = new Genotype_Lin(0, 0, pop_size1, start_ab1, seq, mut_rate1, fitn_coefs1[i]);
        gntp gnts2[max_unique_gnt];
        int start_ab2 = gsl_rng_uniform(rng) * pop_size2;
        gnts2[0] = new Genotype_Lin(0, 0, pop_size2, start_ab2, seq, mut_rate2, fitn_coefs2[i]);
        Population pop1 = Population(0, 1, gnts1, 1, true);
        Population pop2 = Population(1, 0, gnts2, 1, true);
        poppv populations = poppv { &pop1, &pop2 };
        wf.run(traj_len, populations, traj_step2);
        wf.print_trajectory((str)argv[1], std::to_string(i));

        // Iteration over realizations
        for (int r=0; r<n_realizations; r++){

            // Building the populations
            int seq[1] {1};
            gntp gnts1[max_unique_gnt];
            int start_ab1 = gsl_rng_uniform(rng) * pop_size1;
            gnts1[0] = new Genotype_Lin(0, 0, pop_size1, start_ab1, seq, mut_rate1, fitn_coefs1[i]);
            gntp gnts2[max_unique_gnt];
            int start_ab2 = gsl_rng_uniform(rng) * pop_size2;
            gnts2[0] = new Genotype_Lin(0, 0, pop_size2, start_ab2, seq, mut_rate2, fitn_coefs2[i]);

            // Initializing the populations
            Population pop1 = Population(0, 1, gnts1, 1, true);
            Population pop2 = Population(1, 0, gnts2, 1, true);
            poppv populations = poppv { &pop1, &pop2 };

            // Running the WF for n_gen that allow to compute the max dt turnover
            int n_gen = relax_time_base + (n_turn_points+1)*turn_long_dts[turn_long_dts.size()-1];
            // A random initial time is added to remove trivial dependences on init cond
            n_gen += gsl_rng_uniform(rng) * relax_time_rand;
            wf.run(n_gen, populations, traj_step);

            doublev2 traj;
            for (int dti=0; dti<turn_long_dts.size(); dti++){
                int dt = turn_long_dts[dti];
                wf.compute_turnover(traj, dt, n_turn_points);
                for(int i=0; i<n_turn_points; i++){
                    int i0 = n_turn_points * r;
                    turnovers[dti][0][i0+i] = traj[0][i]; 
                    turnovers[dti][1][i0+i] = traj[1][i]; 
                }
            }
        }

        // Printing the observables at given fitness coef
        str is = std::to_string(i);
        for (int dti=0; dti<turn_long_dts.size(); dti++){
            str name = "/turn"+is+"_dt"+std::to_string(int(turn_long_dts[dti]))+".txt";
            print_2dtraj(turnovers[dti], (str)argv[1]+name);
        }

        std::cout << "Iteration " << i+1 << "/" << fitn_coefs1.size();
        std::cout << ", s1: " << fitn_coefs1[i] << ", s2: " << fitn_coefs2[i] << " completed\n";
    }
    
    std::cout << "elapsed time: " << timer.elapsed() << "\n";

    return 0;
}