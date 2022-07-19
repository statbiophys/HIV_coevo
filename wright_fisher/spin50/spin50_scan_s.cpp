#include "../wright_fisher.h"





// Population of binary spins that coevolve

// Compilation string:
// g++ -o spin50_scan_s.exe spin50_scan_s.cpp ../wright_fisher.cpp ../utils.cpp -lgsl -lgslcblas -std=c++17

// Execution string:
// ./spin50_scan.exe [directory where the parameter file is saved]

// Parameter file must be coomposed of lines like:
// prameter_id \t parameter_value


template <typename T>
void delete_pointed_to(T* const ptr)
{
    delete ptr;
}


int main(int argc, char *argv[]){

    if (seq_len != 50)
        throw std::runtime_error("Set sequence length to 50 in wright_fisher.h");

    // Creating and seeding the random number generator
    gsl_rng* rng = get_rand_num_generator();

    // Importing the parameters
    //str path = "scan_s2/";
    //Param params = import_param_from_main(2, path);
    Param params = import_param_from_main(argc, argv);
    int n_realizations = params.get_double("n_realizations");
    int relax_time_base = params.get_double("relax_time_base");
    int relax_time_rand = params.get_double("relax_time_rand");
    int n_turn_points = params.get_double("n_turn_points");
    int n_init_seq = params.get_double("n_init_seq");
    int pop_size1 = params.get_double("pop_size1");
    int pop_size2 = params.get_double("pop_size2");
    int traj_step = params.get_double("traj_step");
    double mut_rate1 = params.get_double("mut_rate1");
    double mut_rate2 = params.get_double("mut_rate2");
    doublev fitn_coefs = params.get_doublev("sel_coefs");
    doublev turn_long_dts = params.get_doublev("turn_long_dt");

    Timer timer = Timer();
    // Iteration over the fitness coefficients of the two populations
    for (int i=0; i<fitn_coefs.size(); i++) {
        for (int j=0; j<fitn_coefs.size(); j++) {

            // Only the turnovers are stored and printed.
            // For each value of dt considered, the turnovers of all the realizations
            // are appended and printed in a single file.
            std::vector<doublev2> turnovers( 
                turn_long_dts.size(), doublev2( 2, doublev(n_turn_points*n_realizations) ) 
            );

            for (int r=0; r<n_realizations; r++){

                // Building the populations
                int seq1[seq_len], seq2[seq_len];
                gntp gnts1[max_unique_gnt], gnts2[max_unique_gnt];
                for(int g=0; g<n_init_seq; g++){
                    generate_rand_seq(rng, seq1);
                    gnts1[g] = new Genotype_Lin(g, g, pop_size1/n_init_seq, 0, seq1, mut_rate1, fitn_coefs[i]);
                    generate_rand_seq(rng, seq2);
                    gnts2[g] = new Genotype_Lin(g, g, pop_size2/n_init_seq, 0, seq1, mut_rate2, -fitn_coefs[j]);
                }

                // Initializing the populations
                Population* pop1 = new Population(0, 1, gnts1, n_init_seq, false);
                Population* pop2 = new Population(1, 0, gnts2, n_init_seq, false);
                poppv populations = poppv { pop1, pop2 };

                Wright_Fisher wf = Wright_Fisher(rng, false);
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
                
                std::for_each(populations.begin(), populations.end(), delete_pointed_to<Population>);
                //std::cout << gnts1[0]->get_id() << "\n";
            }

            // Printing the observables at given fitness coef
            for (int dti=0; dti<turn_long_dts.size(); dti++){
                str is = std::to_string(i);
                str js = std::to_string(j);
                str name = "/turn"+is+"_"+js+"_dt"+std::to_string(int(turn_long_dts[dti]))+".txt";
                print_2dtraj(turnovers[dti], (str)argv[1]+name);
            }

            std::cout << "Iteration " << j+fitn_coefs.size()*i+1 << "/" << fitn_coefs.size()*fitn_coefs.size();
            std::cout << ", s1: " << fitn_coefs[i] << ", s2: " << fitn_coefs[j] << " completed\n";
        }
    }
    
    std::cout << "elapsed time: " << timer.elapsed() << "\n";

    return 0;
}