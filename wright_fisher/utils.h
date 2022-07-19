#ifndef UTILS_H
#define UTILS_H


#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <math.h> 
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <stdio.h>
#include <string>

// Including libraries from GNU Scientific Library (to be 
// installed https://www.gnu.org/software/gsl/doc/html/index.html)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// alias for string
using str = std::string;
// aliases for vectors
using intv = std::vector<int>;
using intv2 = std::vector<intv>;
using doublev = std::vector<double>;
using doublev2 = std::vector<doublev>;
// aliases for dictionaries
using doubled = std::map<str, double>;
using doublevd = std::map<str, doublev>;
using strd = std::map<str, str>;



/* Parameters. They can be doubles, vector of doubles or strings */
class Param {

    private:
	    doubled d;
	    doublevd dv;
	    strd s;

        doublev str2vecd(const str& line, str separator, bool sep_at_end);

    public:
        Param(const str& path);

        double get_double(str param_id) const;
        doublev get_doublev(str param_id) const;
        str get_str(str param_id) const;
};


/* Class for displaying percentage in the standard output. */
class Perc {
	private:
		int m_perc_step;
		int m_max_steps;
		double m_last_perc;

	public:
		Perc(int perc_step, int max_steps);
		void step(int curr_step);
};


/* Class for measuring the time between the reset and en enlapsed call */
class Timer {
	private:
		using clock_t = std::chrono::high_resolution_clock;
		using second_t = std::chrono::duration<double, std::ratio<1> >;
		std::chrono::time_point<clock_t> m_beg;
	public:
		Timer() : m_beg(clock_t::now()) { };
		/* Set the onset for the time measure */
		void reset() { m_beg = clock_t::now(); }
		/* Get the time in seconds enlapsed from reset */
		double elapsed() const { return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count(); }
};


const int max_size=300;

template <typename T, typename U> class dict1000 {
	private:
		int current_size;
		T keys[max_size];
		U values[max_size];

		int find_key(T k, bool rise_err=true){
			for (int i=0; i<current_size; i++)
				if (keys[i] == k)
					return i;
			if (rise_err) std::runtime_error("Key not present in dict1000");
			return current_size;
		}

	public:
		dict1000(){ current_size = 0; }

		U operator[](const T& key) const { return values[find_key(key)]; }

		U& operator[](const T& key) { return values[find_key(key)]; }

		int size() const { return current_size; }

		T get_key(int i) const { 
			if (i>current_size) std::runtime_error("Size of dict1000 exceded");
			return keys[i];
		}

		U get_value(int i) const { 
			if (i>current_size) std::runtime_error("Size of dict1000 exceded");
			return values[i];
		}

		bool is_present (const T& k) {
			int i = find_key(k, false);
			if (i == current_size) return false;
			return true;
		}

		void add(const T& key, const U& value) {
			if (current_size == max_size) 
				std::runtime_error("Max size of dict1000 exceded");
			if (find_key(key, false) != current_size) 
				std::runtime_error("Key already present in dict1000");
			keys[current_size] = key;
			values[current_size] = value;
			current_size++;
		}

		void remove(const T& key){
			int i = find_key(key);
			for (int j=i; j<current_size-1; j++){
				keys[j] = keys[j+1];
				values[j] = values[j+1];
			}
			current_size--;
		}

		void clear() { current_size = 0; }
};



/* Building the random number generator of GSL. If fix_seed=0 the seed is random by default
   otherwise it correspond to fix_seed. */
gsl_rng* get_rand_num_generator(int fix_seed=0);

/* Importing the parameters from the main arguments. It checks if the 
   paramenter file path is given at execution and translates the text
   file in the param struct */
Param import_param_from_main(int argc, char *argv[]);
Param import_param_from_main(int argc, str path);

void print_2dtraj(const doublev2& traj, str file_path);


#endif