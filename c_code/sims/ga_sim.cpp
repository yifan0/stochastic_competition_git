#include <cstdio>
#include <random>	// rand()
#include <string.h> // memcpy()
#include <queue>	// queue
#include <time.h>
#include <chrono>
#include <set>
#include <map>
#include <stack>
#include <numeric>
#include <algorithm>
#include "sfmt.h"
#include "sfmt.cpp"
#include "userintf.cpp"
#include "macdecls.h"
#include "stocc.h"
#include "stoc1.cpp"
#include <mpi.h>
#include <iostream>	   // std::cout
#include <fstream>	   // std::ofstream
#include <cxxopts.hpp> // to handle cmdline args
#include "ga++.h"
#include "tree.h"
#include "measure1D.cpp"
#include "width2D.h"
#include "diff2D.h"
#include "print_stats.h"
#include "print_msg.h"
#include <mxx/reduction.hpp>
using namespace std;

#define GHOSTS 2
#define GA_DIMS 2 // use 2 dims for array, just set one dim to size 1 when using 1D

typedef std::tuple<int, int, cell_type> cell_update;

// Select locations of events based on probability p events in area 0 to n
// return: vector of event locations
set<int> event_list(CRandomSFMT1& rng, StochasticLib1& srng, size_t n, double p) {
	size_t count = srng.Binomial(n, p);
	set<int> events;
	while (events.size() < count) {
		int index = rng.IRandomX(0,n-1);
		if(!events.count(index)) {
			events.insert(index);
		}
	}
	return events;
}


int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);
	int mpi_num_procs = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
	GA_Initialize();
	char init_err_msg[] = "MA_init failed";
	int heap = 3000000, stack = 3000000;
	if (!MA_init(C_DBL, stack, heap))
		GA_Error(init_err_msg, stack + heap);
	int me = GA_Nodeid();
	int nprocs = GA_Nnodes();

	std::chrono::time_point<std::chrono::system_clock> start_time;
	std::chrono::time_point<std::chrono::system_clock> end_time, rep_start_time, rep_end_time, out_start_time, out_end_time;

	cxxopts::Options options("mpirun -n N ./ga_sim", "Stocastic competition simulation on a grid. Each cell contains a single individual with a trait value. Each step, individuals reproduce with mutation and speciation, and compete with neighbors. Output is the trait value across the grid, and optionally a tree file of speciation events.");

	options.add_options()
		("s,size", "grid size (size^dims)", cxxopts::value<int>()->default_value("500"))
		("t,tree", "generate tree file", cxxopts::value<bool>()->default_value("false"))
		("c,specrate", "speciation rate", cxxopts::value<double>()->default_value("0.0001"))
		("d,dims", "dimension of grid (size^dims)", cxxopts::value<int>()->default_value("2"))
		("p,patch_count", "number of individuals per patch (inverse of p)", cxxopts::value<double>()->default_value("1"))
		("r,reps", "number of repetitions", cxxopts::value<int>()->default_value("1"))
		("m,mutsize", "maximum change in mutation event", cxxopts::value<double>()->default_value("0.1"))
		("f,file_out", "output file on/off", cxxopts::value<bool>()->default_value("false"))
		("o,outfile", "output file prefix (for csv and/or tree files)", cxxopts::value<std::string>()->default_value("out"))
		("n,nsteps", "number of steps to run", cxxopts::value<int>()->default_value("100"))
		("u,diff", "run difference statistic", cxxopts::value<bool>()->default_value("false"))
		("v,div", "run diversity statistic", cxxopts::value<bool>()->default_value("false"))
		("w,width", "run width statistic", cxxopts::value<bool>()->default_value("false"))
        ("g,ghosts", "frequency to update ghosts (synchronize every g steps)", cxxopts::value<int>()->default_value("1"))
		("h,help", "Print usage")
		;

	auto result = options.parse(argc, argv);

	if (result.count("help")) {
		if (me == 0)
			std::cout << options.help() << std::endl;
		MPI_Finalize();
		exit(0);
	}
	int size = result["size"].as<int>();
	int nrep = result["reps"].as<int>();
	double specrate = result["specrate"].as<double>();
	double mutsize = result["mutsize"].as<double>();
	std::string outfile = result["outfile"].as<std::string>();
	int nsteps = result["nsteps"].as<int>();
	double p = 1 / result["patch_count"].as<double>();
	int timescale = 100 * (size * 1.0 / p);
	int endtime = timescale / nsteps;
    int ghost_sync_period = result["ghosts"].as<int>();
	int ndims = result["dims"].as<int>();
	if (ndims != 1 && ndims != 2) {
		if (me == 0)
			println("Only 1 and 2 dimensions are supported.");
		MPI_Finalize();
		exit(1);
	}
    // always use a "2D" array, but for 1D have the second dimension be 1
	int dims[GA_DIMS];
	int grid_ld[GA_DIMS];
	int lo[GA_DIMS], hi[GA_DIMS];
	int ghost_grid_ld[GA_DIMS - 1];
	int ghost_dims[GA_DIMS];
	int ghost_width[GA_DIMS];
	double *land_grid_ptr, *ghost_grid_ptr;
	FILE *fp;

	// boolean options
	bool generate_tree = result["tree"].as<bool>();
	bool file_out = result["file_out"].as<bool>();
	bool run_div_stat = result["div"].as<bool>();
	bool run_width_stat = result["width"].as<bool>();
	bool run_diff_stat = result["diff"].as<bool>();
    if (ndims == 1 && (run_div_stat || run_width_stat || run_diff_stat)) {
        run_div_stat = true;
        run_width_stat = true;
        run_diff_stat = true;
    }
    if (ndims == 2 && (run_div_stat || run_diff_stat)) {
        run_div_stat = true;
        run_diff_stat = true;
    }
    if (ndims == 2 && (run_div_stat || run_diff_stat || run_width_stat)) {
        file_out = true;
    }

	start_time = std::chrono::system_clock::now();

	// grid for average across reps
	GA_Mask_sync(0, 0); // turns off sync when updating ghosts
    dims[1] = size;
    ghost_width[1] = GHOSTS;
    int total_cells = size;
    if (ndims == 2) {
        dims[0] = size;
        ghost_width[0] = GHOSTS;
        total_cells = size*size;
    }
    else {
        dims[0] = 1;
        ghost_width[0] = 0;
    }
	char land_grid_name[] = "land grid";
	int ga_land_grid = NGA_Create_ghosts(C_DBL, GA_DIMS, dims, ghost_width, land_grid_name, NULL);
	if (ga_land_grid == 0) {
		char create_err[] = "Failure for NGA_Create_ghosts()";
		GA_Error(create_err, 1);
	}

	NGA_Distribution(ga_land_grid, me, lo, hi);
	if (lo[0] < 0) {
		debug("Owns no elements");
		return 1;
	}
	NGA_Access(ga_land_grid, lo, hi, &land_grid_ptr, grid_ld);

	println("Inputs:");
	println("\trepetitions = %d", nrep);
	println("\tsize = %d", size);
	println("\tdimensions = %d", ndims);
	println("\tindividuals per patch = %f", 1 / p);
	println("\tmutation size = %f", mutsize);
	println("\tspeciation rate = %.2e", specrate);
	println("\ttimescale = %d", timescale);
	println("\tnsteps = %d", nsteps);
	println("\tend time = %d", endtime);
	println("\tprocesses = %d", mpi_num_procs);
	println("\tgenerate tree = %s", generate_tree ? "true" : "false");
	println("\tfile output = %s", file_out ? "true" : "false");
	println("\tstatistics run: diversity = %s, width = %s, difference = %s", run_div_stat ? "true" : "false", run_width_stat ? "true" : "false", run_diff_stat ? "true" : "false");
    if (ndims == 1 && run_div_stat) {
        println("NOTE: for 1D, either all or no statistics are run.");
    } else if (ndims == 2 && run_div_stat) {
        println("NOTE: for 2D, either diversity and difference are run, or neither is run. Width is turned on/off separately.");
    }
    if (ndims == 2 && (run_div_stat || run_width_stat || run_diff_stat)) {
        println("NOTE: for 2D, to calculate summary statistics, output file must be used.");
    }
	println("");
	fflush(stdout);

	// Random number generation
	CRandomSFMT1 RanGen(time(0) + me * 10); // Agner Combined generator
	StochasticLib1 sto(time(0) + me * 7);	// Stochastic RNG
	int local_rows = (hi[0] - lo[0] + 1);
	int local_cols = (hi[1] - lo[1] + 1);
    //println("local_rows = %d, local_cols = %d", local_rows, local_cols);
	int local_area = local_cols * local_rows;
	double land_mask_data[local_area];
	double* land_mask[local_rows];
	for (size_t i = 0; i < local_rows; i++)
		land_mask[i] = land_mask_data + i * local_cols;

	for (int rep = 0; rep < nrep; rep++) {
		rep_start_time = std::chrono::system_clock::now();

		double one = 1;
		GA_Fill(ga_land_grid, &one);
		GA_Update_ghosts(ga_land_grid);
		NGA_Access_ghosts(ga_land_grid, ghost_dims, &ghost_grid_ptr, ghost_grid_ld);
		if (!ghost_grid_ptr) {
			char grid_err[] = "NULL pointer for ghost grid.";
			GA_Error(grid_err, 1);
		}

		// zero out mask except local edges
		memset(land_mask_data, 0, sizeof(*land_mask_data));
		for (size_t i = 0; i < local_rows; i++) {
			land_mask[i][0] = 1;
			land_mask[i][local_cols - 1] = 1;
		}
		for (size_t i = 0; i < local_cols; i++) {
			land_mask[0][i] = 1;
			land_mask[local_rows - 1][i] = 1;
		}

		// invasion rule variables
		cell_type neighborhood[8];
		cell_type inv[8];
		cell_type inv_sum = 0;
		int inv_index = 0;

		// track speciation events
		vector<tuple<size_t, cell_type, cell_type>> speciation_events; // only used when generate_tree == true

		int row, col;
		for (int step = 0; step < timescale; step++) {

			// speciation rule
			set<int> spec_events = event_list(RanGen, sto, local_area, specrate);
			for(int index : spec_events) {
				size_t i = index / local_cols; // row
				size_t j = index % local_cols; // col
				float ratio = (1 + RanGen.Random() * mutsize);
				if (RanGen.Random() < 0.5) ratio = 1 / ratio;
				float probsuccess = p * ratio / (p * (ratio - 1) + 1);
				if (RanGen.Random() <= probsuccess) {
					cell_type old_val = ghost_grid_ptr[(i + ghost_width[0]) * ghost_grid_ld[0] + j + ghost_width[1]];
					cell_type new_val = old_val * ratio;
					ghost_grid_ptr[(i + ghost_width[0]) * ghost_grid_ld[0] + j + ghost_width[1]] = new_val;
					if (generate_tree)
						speciation_events.push_back(make_tuple(step, old_val, new_val));
				}

				// set mask for [i][j] and surrounding cells unless on edge
				for (int x = -1; x <= 1; x++) {
					for (int y = -1; y <= 1; y++) {
						row = i + x;
						col = j + y;
						if (row >= local_rows - 1 || row <= 0)
							continue;
						if (col >= local_cols - 1 || col <= 0)
							continue;
						cell_type local_max = ghost_grid_ptr[(row + ghost_width[0]) * ghost_grid_ld[0] + col + ghost_width[1]];
						for (int xx = -1; xx <= 1; xx++) {
							for (int yy = -1; yy <= 1; yy++) {
								int neighbor_row = row + xx;
								int neighbor_col = col + yy;
								local_max = std::max(local_max, ghost_grid_ptr[(neighbor_row + ghost_width[0]) * ghost_grid_ld[0] + neighbor_col + ghost_width[1]]);
							}
						}
						if(row >= 0 && col >= 0 && row < local_rows && col < local_cols)
							land_mask[row][col] = local_max * p / (local_max * p + ghost_grid_ptr[(row + ghost_width[0]) * ghost_grid_ld[0] + col + ghost_width[1]] * (1 - p));
					}
				}
			}

			if (step % ghost_sync_period == 0)
				GA_Update_ghosts(ga_land_grid);

			// invasion rule
			std::vector<cell_update> updates;
			for (size_t i = 0; i < local_rows; i++) {
				for (size_t j = 0; j < local_cols; j++) {
					if (land_mask[i][j] != 0) {
						double randval = RanGen.Random(); //random_float();
						if (randval < land_mask[i][j]) {
							inv_sum = 0;
							inv_index = 0;
							for (int x = -1; x <= 1; x++) {
								for (int y = -1; y <= 1; y++) {
									if ((x != 0 || y != 0)) {
										neighborhood[inv_index] = ghost_grid_ptr[(i + x + ghost_width[0]) * ghost_grid_ld[0] + j + y + ghost_width[1]];
										inv[inv_index] = p * neighborhood[inv_index] / (p * neighborhood[inv_index] + ghost_grid_ptr[(i + ghost_width[0]) * ghost_grid_ld[0] + j + ghost_width[1]] * (1 - p));
										inv_sum += inv[inv_index];
										inv_index++;
									}
								}
							}
							if (randval <= inv_sum / 8) {
								// Get random element with weighted probabilities
								double weighted_rand = RanGen.Random() * inv_sum;
								inv_index = 0;
								while (weighted_rand > inv[inv_index]) {
									weighted_rand -= inv[inv_index];
									inv_index++;
								}
								if (neighborhood[inv_index] != ghost_grid_ptr[(i + ghost_width[0]) * ghost_grid_ld[0] + j + ghost_width[1]]) {
									updates.push_back({i, j, neighborhood[inv_index]});
								}
							}
						}
					}
				}
			}

			for (const auto &[i, j, val] : updates) {
				ghost_grid_ptr[(i + ghost_width[0]) * ghost_grid_ld[0] + j + ghost_width[1]] = val;
				for (int x = -1; x <= 1; x++) {
					for (int y = -1; y <= 1; y++) {
						row = i + x;
						col = j + y;
						if (row >= local_rows - 1 || row <= 0)
							continue;
						if (col >= local_cols - 1 || col <= 0)
							continue;
						cell_type local_max = 0;
						for (int xx = -1; xx <= 1; xx++) {
							for (int yy = -1; yy <= 1; yy++) {
								row = i + x + xx;
								col = j + y + yy;
								if (row >= local_rows - 1 || row <= 0)
									continue;
								if (col >= local_cols - 1 || col <= 0)
									continue;
								if (ghost_grid_ptr[(row + ghost_width[0]) * ghost_grid_ld[0] + col + ghost_width[1]] != val && ghost_grid_ptr[(row + ghost_width[0]) * ghost_grid_ld[0] + col + ghost_width[1]] > local_max) {
									local_max = ghost_grid_ptr[(row + ghost_width[0]) * ghost_grid_ld[0] + col + ghost_width[1]];
								}
							}
						}
						row = i + x;
						col = j + y;
						if (local_max == 0) {
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols)
								land_mask[row][col] = 0;
						}
						else {
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols)
								land_mask[row][col] = local_max * p / (local_max * p + ghost_grid_ptr[(row + ghost_width[0]) * ghost_grid_ld[0] + col + ghost_width[1]] * (1 - p));
						}
					}
				}
			}

			// renormalize every nstep steps
			if (step % endtime == endtime - 1) {
				cell_type land_grid_mean = 0;
				for (size_t i = 0; i < local_rows; i++)
					for (size_t j = 0; j < local_cols; j++)
						land_grid_mean += ghost_grid_ptr[(i + ghost_width[0]) * ghost_grid_ld[0] + j + ghost_width[1]];
				land_grid_mean = land_grid_mean / (total_cells);
				MPI_Allreduce(MPI_IN_PLACE, &land_grid_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				println("Global average at step %d = %f", step, land_grid_mean);
				if (land_grid_mean > 100) {
					if (me == 0)
						println("Normalizing by %f", land_grid_mean);

					if (generate_tree) {
						for (tuple<size_t, cell_type, cell_type> &event : speciation_events) {
							get<1>(event) = get<1>(event) / land_grid_mean;
							get<2>(event) = get<2>(event) / land_grid_mean;
						}
					}

					for (size_t i = 0; i < ghost_dims[0]; i++)
						for (size_t j = 0; j < ghost_dims[1]; j++)
							ghost_grid_ptr[i * ghost_dims[1] + j] /= land_grid_mean;
				}
			}
		}

		rep_end_time = std::chrono::system_clock::now();
		std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
		println("Rep %d run time = %fs", rep, rep_time.count());
		fflush(stdout);

		out_start_time = std::chrono::system_clock::now();
		if (file_out) {
			outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".csv";
			fp = fopen(outfile.c_str(), "w");
			GA_Print_csv_file(fp, ga_land_grid);
			fclose(fp);
            println("Wrote results to file %s", outfile.c_str());
		}

		// gather speciation events
		if (generate_tree) {
			vector<tuple<size_t, cell_type, cell_type>> global_speciation_events = mxx::gatherv(speciation_events, 0);
			// gather extant species set
			debug("Gathering extant species list");
			set<cell_type> extant_species(ghost_grid_ptr+ghost_width[0]*ghost_grid_ld[0]+ghost_width[1], ghost_grid_ptr+(local_rows + ghost_width[0]) * ghost_grid_ld[0] + local_cols + ghost_width[1]);
			vector<cell_type> extant_species_vec(extant_species.begin(), extant_species.end());
			debug("Local extant species: %lu", extant_species.size());
			vector<cell_type> global_extant_species_vec = mxx::gatherv(extant_species_vec, 0);
			if (me == 0) {
				set<cell_type> global_extant_species(global_extant_species_vec.begin(), global_extant_species_vec.end());
				debug("Global extant species: %lu", global_extant_species.size());
				map<cell_type, speciation_tree_node*> leaves;
				// sort global_speciation_events by timestep
				sort(global_speciation_events.begin(), global_speciation_events.end());
				speciation_tree_node *tree = new speciation_tree_node(get<1>(global_speciation_events[0]), 0, nullptr);
				leaves[get<1>(global_speciation_events[0])] = tree;
				debug("Creating tree");
				for (tuple<size_t, cell_type, cell_type> event : global_speciation_events) {
					size_t timestep = get<0>(event);
					cell_type old_val = get<1>(event);
					cell_type new_val = get<2>(event);
					// instead of finding the node by traversing the tree, just look it up in a hash map (for performance)
					// speciation_tree_node *parent_node = find_node(tree, old_val);
					speciation_tree_node *parent_node = leaves[old_val];
					if (!parent_node) {
						char find_err[] = "Did not find node with value";
						debug("Error: did not find node with value %f", old_val);
						debug("%s", toString_final(tree, timescale).c_str());
						GA_Error(find_err, 1);
						return 1;
					}
					speciation_tree_node *old_child = new speciation_tree_node(old_val, timestep, parent_node);
					speciation_tree_node *new_child = new speciation_tree_node(new_val, timestep, parent_node);
					parent_node->left_child = old_child;
					parent_node->right_child = new_child;
					leaves[old_val] = old_child;
					leaves[new_val] = new_child;
				}

				// prune extinct species
				debug("Pruning tree");
				tree = prune(tree, global_extant_species);
				debug("Pruned, root time = %lu", tree->time);

				if (file_out) {
					outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".tree";
					ofstream fout(outfile.c_str());
					fout << toString_final(tree, timescale) << endl;
					delete_tree(tree);
					fout.close();
                    println("Wrote results to file %s", outfile.c_str());
				}
			}
		}

		if (run_div_stat || run_width_stat || run_diff_stat) {
			std::string div_file = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + "_div.csv";
			std::string diff_file = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + "_diff.csv";;
			std::string width_file = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + "_width.csv";;
			if (ndims == 1) {
				vector<double> local_landscape_vec(land_grid_ptr, land_grid_ptr+local_cols);
				std::string landscape_str(local_landscape_vec.begin(), local_landscape_vec.end());
				vector<double> landscape = mxx::gatherv(local_landscape_vec, 0);
				if (me == 0) {
					vector<double> diversity;
					vector<double> difference;
					vector<double> width;
					measure1D(landscape, diversity, difference, width);
					print_double_stat_to_csv(diversity, div_file);
					println("Saved diversity to %s", div_file.c_str());
					print_double_stat_to_csv(difference, diff_file);
					println("Saved difference to %s", diff_file.c_str());
					print_double_stat_to_csv(width, width_file);
					println("Saved width to %s", width_file.c_str());
				}
			} else if (ndims == 2) {
                // 2D stats are run from CSV file
                if (run_width_stat) {
                    if (me == 0) {
                        std::string species_count_file = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + "_species_count.csv";
                        std::string fitness_file = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + "_fitness.csv";
                        std::string slope_file = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + "_slope.csv";
                        tuple<array<double,10>,vector<int>,vector<double>,double> result = width2D(outfile);
                        array<double,10> width_arr = get<0>(result);
                        vector<int> species_count = get<1>(result);
                        vector<double> species_fitness = get<2>(result);
                        double coeff = get<3>(result);

                        print_int_stat_to_csv(species_count, species_count_file.c_str());
                        println("Saved species count to %s", species_count_file.c_str());
                        print_double_to_csv(coeff, slope_file.c_str());
                        println("Saved slope to %s", slope_file.c_str());
                        print_double_stat_to_csv(species_fitness, fitness_file.c_str());
                        println("Saved species fitness to %s", fitness_file.c_str());
                        print_double_arr_to_csv(width_arr.data(), width_arr.size(), width_file.c_str());
					    println("Saved width to %s", width_file.c_str());
                    }
                }
                if (run_div_stat) { // div and diff are either both run or neither are run
                    if (me == 0) {
                        tuple<array<double,7>,array<double,7>> result = diff2D(outfile);
                        array<double,7> diff_arr = get<0>(result);
                        array<double,7> div_arr = get<1>(result);
                        print_double_arr_to_csv(diff_arr.data(), diff_arr.size(), diff_file);
                        println("Saved difference to %s", diff_file.c_str());
                        print_double_arr_to_csv(div_arr.data(), div_arr.size(), div_file);
                        println("Saved diversity to %s", div_file.c_str());
                    }
                }
			} else {
				println("Summary statistics unavailable for %d dimensions", ndims);
			}
		}

		out_end_time = std::chrono::system_clock::now();
		std::chrono::duration<double> out_time = out_end_time - out_start_time;
		println("Output run time = %fs", out_time.count());
	}

	end_time = std::chrono::system_clock::now();
	std::chrono::duration<double> total_time = end_time - start_time;
	println("Total run time = %fs", total_time.count());
	fflush(stdout);

	GA_Destroy(ga_land_grid);

	GA_Terminate();
	MPI_Finalize();
}
