#include "utils.h"


gsl_rng* get_rand_num_generator(int fix_seed){
    int seed = fix_seed;
    if (fix_seed == 0)
        seed = std::chrono::system_clock::now().time_since_epoch().count();
    const gsl_rng_type* T;
    gsl_rng* r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, seed);
    return r;
}


Param import_param_from_main(int argc, char *argv[]) {
    if (argc != 2)
        throw std::runtime_error("The path containing the parameter file must be passed during execution");

    str path(argv[1]);
    return Param(path + "/param.txt");
}

Param import_param_from_main(int argc, str path) {
    if (argc != 2)
        throw std::runtime_error("The path containing the parameter file must be passed during execution");

    return Param(path + "/param.txt");
}


Param::Param(const str& file_path){

    d = doubled();
    dv = doublevd();
    s = strd();
    
    std::ifstream param_file (file_path);
    if (!param_file.is_open())
        throw std::runtime_error("Error in opening the parameter file at "+file_path);

    str line;
    while ( getline (param_file, line) ) {
        std::size_t tab_pos = line.find("\t");
        str key = line.substr(0,tab_pos);
        str value = line.substr(tab_pos+1, str::npos);

        if (value.find(",") != str::npos){
            dv[key] = str2vecd(value, ",", false); // Parse a vector
        }
        else{
            try {
                double vald = std::stod(value); // Parse a double
                d[key] = vald;
            } catch (std::invalid_argument const&){
                s[key] = value; // Parse a string if stod gives exception
            }
        }
    }
    param_file.close();

    if (d.size() == 0 && dv.size() == 0 && s.size() == 0)
        throw std::runtime_error("Empty parameter file");
}


double Param::get_double(str param_id) const {
    if (d.find(param_id) == d.end())
        throw std::runtime_error("Parameter " + param_id + " not found.");
    else return d.at(param_id);
}

doublev Param::get_doublev(str param_id) const {
    if (dv.find(param_id) == dv.end())
        throw std::runtime_error("Parameter " + param_id + " not found.");
    else return dv.at(param_id);
}

str Param::get_str(str param_id) const {
    if (s.find(param_id) == s.end())
        throw std::runtime_error("Parameter " + param_id + " not found.");
    else return s.at(param_id);
}


doublev Param::str2vecd(const str& line, str separator, bool sep_at_end) {
    std::size_t sep_pos = line.find(separator);
    if (sep_pos == str::npos) {
        if (sep_at_end)
            throw std::runtime_error(separator + " separator not found in " + line);
        else {
            doublev v = {std::stod(line.substr(0, sep_pos))};
            return v;
        }
    }

    str elem = line.substr(0, sep_pos);
    doublev v = doublev(0);
    try {
        v.push_back(std::stod(elem));
    }
    catch (std::exception& e){
        v.push_back(0);
    }

    while (true){
        std::size_t next_sep_pos = line.find(separator, sep_pos+1);
        if (sep_at_end && next_sep_pos == str::npos) break;
        str elem = line.substr(sep_pos+1, next_sep_pos-sep_pos);
        try{
            v.push_back(std::stod(elem));
        }
        catch (std::exception& e){
            v.push_back(0);
        }
        if (!sep_at_end && next_sep_pos == str::npos) break;
        sep_pos = next_sep_pos;
    }

    return v;
}


Perc::Perc(int perc_step, int max_steps) : m_max_steps{max_steps}, m_last_perc{0.0} {
    if (perc_step<1 || perc_step>100) {
        std::cout << "Invalid percentage step. Set by default to 10%\n";
        m_perc_step = 10;
    }
    else m_perc_step = perc_step;
};


void Perc::step(int curr_step) {
    double perc = (double)curr_step/(double)m_max_steps*100;
    if (perc >= m_last_perc){
        std::cout << round(perc) << "%\n";
        m_last_perc = round(perc) + m_perc_step;
    }
}


void print_2dtraj(const doublev2& traj, str file_path) {
    std::ofstream out;
    out.open(file_path);

    for (int t=0; t<traj[0].size(); t++){
        for (int p=0; p<traj.size(); p++)
            out << traj[p][t] << "\t";
        out << "\n";
    }

    out.close();
}