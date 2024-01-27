#include <cstdio>
#include <string.h> // memcpy()
#include <queue>	// queue
#include <omp.h>
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
#include <mxx/reduction.hpp>
using namespace std;

#define RANDOM_FLOAT RanGen.Random()

// #define println(...) { if(me == 0) { printf(__VA_ARGS__); printf("\n"); } }
// #define debug(...) { if(true) { printf("Rank %d: ", me); printf(__VA_ARGS__); printf("\n"); fflush(stdout); } }
// #define print(...) { printf(__VA_ARGS__); }

#define NDIM 2
#define GHOSTS 2

typedef std::tuple<int, int, cell_type> cell_update;

// Select locations of events based on probability p events in area 0 to n
// return: vector of event locations
set<int> event_list(CRandomSFMT1& rng, StochasticLib1& stoc_rng, size_t n, double p) {
	size_t count = stoc_rng.Binomial(n, p);
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
	int nrep = 10;
	int size = 100;
	double p = 0.1;
	double mutsize = 0.1;
	double specrate = 0.0001; // 1.0e-4
	int timescale = 100 * size / p;
	int nsteps = 100;
	int endtime = timescale / nsteps;
	std::string outfile;
	std::chrono::time_point<std::chrono::system_clock> start_time, end_time, rep_start_time, rep_end_time, out_start_time, out_end_time;
	int me, nprocs;
	int dims[NDIM];
	int grid_ld[NDIM];
	int lo[NDIM], hi[NDIM];
	int ghost_grid_ld[NDIM - 1];
	int ghost_dims[NDIM];
	int ghost_width[NDIM];
	double *land_grid_ptr, *ghost_grid_ptr;
	int heap = 3000000, stack = 3000000;
	FILE *fp;

	MPI_Init(&argc, &argv);
	int num_procs = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	GA_Initialize();
	char init_err_msg[] = "MA_init failed";
	if (!MA_init(C_DBL, stack, heap))
		GA_Error(init_err_msg, stack + heap);

	me = GA_Nodeid();
	nprocs = GA_Nnodes();

	cxxopts::Options options("mpirun -n N ./ga_sim", "Stocastic competition simulation with GlobalArray parallelism");

	options.add_options()
		("s,size", "grid size", cxxopts::value<int>()->default_value("500"))
		("c,specrate", "speciation rate", cxxopts::value<double>()->default_value("0.0001"))
		("r,reps", "number of repetitions", cxxopts::value<int>()->default_value("1"))
		("m,mutsize", "maximum change in mutation event", cxxopts::value<double>()->default_value("0.1"))
		("o,outfile", "output file prefix", cxxopts::value<std::string>()->default_value("out"))
		("h,help", "Print usage")
		;

	auto result = options.parse(argc, argv);

	if (result.count("help")) {
		if (me == 0)
			std::cout << options.help() << std::endl;
		MPI_Finalize();
		exit(0);
	}
	size = result["size"].as<int>();
	nrep = result["reps"].as<int>();
	specrate = result["specrate"].as<double>();
	mutsize = result["mutsize"].as<double>();

	timescale = 100 * (size * 1.0 / p);
	endtime = timescale / nsteps;

	start_time = std::chrono::system_clock::now();

	// grid for average across reps
	GA_Mask_sync(0, 0); // turns off sync when updating ghosts
	for (size_t i = 0; i < NDIM; i++) {
		dims[i] = size;
		ghost_width[i] = GHOSTS;
	}
	char land_grid_name[] = "land grid";
	int ga_land_grid = NGA_Create_ghosts(C_DBL, NDIM, dims, ghost_width, land_grid_name, NULL);
	if (ga_land_grid == 0) {
		char create_err[] = "Failure for NGA_Create_ghosts()";
		GA_Error(create_err, 1);
	}

	NGA_Distribution(ga_land_grid, me, lo, hi);
	if (lo[0] < 0) {
		debug("Owns no elements");
	}
	NGA_Access(ga_land_grid, lo, hi, &land_grid_ptr, grid_ld);

	println("Inputs:");
	println("\trepetitions = %d", nrep);
	println("\tsize = %d%s%d", size, "x", size);
	println("\tindividuals per patch = %f", 1 / p);
	println("\tmutation size = %f", mutsize);
	println("\tspeciation rate = %.2e", specrate);
	println("\ttimescale = %d", timescale);
	println("\tnsteps = %d", nsteps);
	println("\tend time = %d", endtime);
	println("\tprocesses = %d", num_procs);
	fflush(stdout);

	cell_type land_grid_mean;

	int local_rows = (hi[0] - lo[0] + 1);
	int local_cols = (hi[1] - lo[1] + 1);
	int local_area = local_cols * local_rows;
	double land_mask_data[local_area];
	double* land_mask[local_rows];
	vector<tuple<size_t, cell_type, cell_type>> process_speciation_events;
	vector<tuple<size_t, cell_type, cell_type>> global_speciation_events;

#pragma omp parallel
{
	#pragma omp single
	println("threads  = %d", omp_get_num_threads());

	CRandomSFMT1 RanGen(time(0) + num_procs * 10 + omp_get_thread_num()); // Agner Combined generator
	StochasticLib1 sto(time(0) + me * 7);	// Stochastic RNG

	#pragma omp for
	for (size_t i = 0; i < local_rows; i++)
		land_mask[i] = land_mask_data + i * local_cols;

	for (int rep = 0; rep < nrep; rep++) {
		rep_start_time = std::chrono::system_clock::now();

		double one = 1;
		#pragma omp single
		{
		GA_Fill(ga_land_grid, &one);
		GA_Update_ghosts(ga_land_grid);
		NGA_Access_ghosts(ga_land_grid, ghost_dims, &ghost_grid_ptr, ghost_grid_ld);
		}
		if (!ghost_grid_ptr) {
			char grid_err[] = "NULL pointer for ghost grid.";
			GA_Error(grid_err, 1);
		}

		// zero out mask except local edges
		#pragma omp single
		memset(land_mask_data, 0, sizeof(*land_mask_data));
		#pragma omp for
		for (size_t i = 0; i < local_rows; i++) {
			land_mask[i][0] = 1;
			land_mask[i][local_cols - 1] = 1;
		}
		#pragma omp for
		for (size_t i = 0; i < local_cols; i++) {
			land_mask[0][i] = 1;
			land_mask[local_rows - 1][i] = 1;
		}
		#pragma omp barrier

		// invasion rule variables
		cell_type neighborhood[8];
		cell_type inv[8];
		cell_type inv_sum = 0;
		int inv_index = 0;

		// track speciation events, speciation_events is local to each thread
		vector<tuple<size_t, cell_type, cell_type>> speciation_events;
		#pragma omp single
		process_speciation_events.clear();
		global_speciation_events.clear();

		int row, col;
		for (int step = 0; step < timescale; step++) {

			// speciation rule
			#pragma omp for
			for(size_t i = 0; i < local_rows; i++) { // i as row
				set<int> spec_events = event_list(RanGen, sto, local_cols, specrate);
				for(int index : spec_events) {
					size_t j = index; // col
					float ratio = (1 + RANDOM_FLOAT * mutsize);
					if (RANDOM_FLOAT < 0.5) ratio = 1 / ratio;
					float probsuccess = p * ratio / (p * (ratio - 1) + 1);
					if (RANDOM_FLOAT <= probsuccess) {
						cell_type old_val = ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS];
						cell_type new_val = old_val * ratio;
						ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS] = new_val;
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
							cell_type local_max = ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS];
							for (int xx = -1; xx <= 1; xx++) {
								for (int yy = -1; yy <= 1; yy++) {
									int neighbor_row = row + xx;
									int neighbor_col = col + yy;
									local_max = std::max(local_max, ghost_grid_ptr[(neighbor_row + GHOSTS) * ghost_grid_ld[0] + neighbor_col + GHOSTS]);
								}
							}
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols) {
								#pragma omp atomic write
								land_mask[row][col] = local_max * p / (local_max * p + ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] * (1 - p));;
							}
						}
					}
				}
			}

			#pragma omp single
			if (step % 10 == 0) {
				GA_Update_ghosts(ga_land_grid);
			}

			// invasion rule
			std::vector<cell_update> updates;
			#pragma omp for
			for (size_t i = 0; i < local_rows; i++) {
				for (size_t j = 0; j < local_cols; j++) {
					if (land_mask[i][j] != 0) {
						double randval = RANDOM_FLOAT;
						if (randval < land_mask[i][j]) {
							inv_sum = 0;
							inv_index = 0;
							for (int x = -1; x <= 1; x++) {
								for (int y = -1; y <= 1; y++) {
									if ((x != 0 || y != 0)) {
										neighborhood[inv_index] = ghost_grid_ptr[(i + x + GHOSTS) * ghost_grid_ld[0] + j + y + GHOSTS];
										inv[inv_index] = p * neighborhood[inv_index] / (p * neighborhood[inv_index] + ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS] * (1 - p));
										inv_sum += inv[inv_index];
										inv_index++;
									}
								}
							}
							if (randval <= inv_sum / 8) {
								// Get random element with weighted probabilities
								double weighted_rand = RANDOM_FLOAT * inv_sum;
								inv_index = 0;
								while (weighted_rand > inv[inv_index]) {
									weighted_rand -= inv[inv_index];
									inv_index++;
								}
								if (neighborhood[inv_index] != ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS]) {
									updates.push_back({i, j, neighborhood[inv_index]});
								}
							}
						}
					}
				}
			}
		
			// every thread has a private updates
			cell_type local_max;
			for (const auto &[i, j, val] : updates) {
				ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS] = val;
				for (int x = -1; x <= 1; x++) {
					for (int y = -1; y <= 1; y++) {
						row = i + x;
						col = j + y;
						if (row >= local_rows - 1 || row <= 0)
							continue;
						if (col >= local_cols - 1 || col <= 0)
							continue;
						local_max = 0;
						for (int xx = -1; xx <= 1; xx++) {
							for (int yy = -1; yy <= 1; yy++) {
								row = i + x + xx;
								col = j + y + yy;
								if (row >= local_rows - 1 || row <= 0)
									continue;
								if (col >= local_cols - 1 || col <= 0)
									continue;
								if (ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] != val && ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] > local_max) {
									local_max = ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS];
								}
							}
						}
						row = i + x;
						col = j + y;
						if (local_max == 0) {
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols) {
								#pragma omp atomic write
								land_mask[row][col] = 0;
							}
						}
						else {
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols) {
								#pragma omp atomic write
								land_mask[row][col] = local_max * p / (local_max * p + ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] * (1 - p));
							}
						}
					}
				}
			}

			// renormalize every nstep steps
			if (step % endtime == endtime - 1) {
				land_grid_mean = 0;
				#pragma omp for reduction(+:land_grid_mean)
				for (size_t i = 0; i < local_rows; i++)
					for (size_t j = 0; j < local_cols; j++)
						land_grid_mean += ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS];
				#pragma omp single
				{
				land_grid_mean = land_grid_mean / (size * size);
				MPI_Allreduce(MPI_IN_PLACE, &land_grid_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				}
				if (land_grid_mean > 100) {
					if (me == 0) {
						#pragma omp single
						println("Normalizing by %f", land_grid_mean);
					}

					for (tuple<size_t, cell_type, cell_type> &event : speciation_events) {
						get<1>(event) = get<1>(event) / land_grid_mean;
						get<2>(event) = get<2>(event) / land_grid_mean;
					}

					#pragma omp for
					for (size_t i = 0; i < ghost_dims[0]; i++)
						for (size_t j = 0; j < ghost_dims[1]; j++)
							ghost_grid_ptr[i * ghost_dims[1] + j] /= land_grid_mean;
				}
			}
		}

		#pragma omp single
		{
		rep_end_time = std::chrono::system_clock::now();
		std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
		println("Rep %d run time = %fs", rep, rep_time.count());
		fflush(stdout);
		
		out_start_time = std::chrono::system_clock::now();
		outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".csv";
		fp = fopen(outfile.c_str(), "w");
		GA_Print_csv_file(fp, ga_land_grid);
		fclose(fp);
		}

		// gather speciation events
		#pragma omp critical
		{
			process_speciation_events.insert(process_speciation_events.end(), speciation_events.begin(), speciation_events.end());
		}
		#pragma omp barrier
		#pragma omp single
		{
		global_speciation_events = mxx::gatherv(process_speciation_events, 0);
		// gather extant species set
		set<cell_type> extant_species(ghost_grid_ptr+GHOSTS*ghost_grid_ld[0]+GHOSTS, ghost_grid_ptr+(local_rows + GHOSTS) * ghost_grid_ld[0] + local_cols + GHOSTS);
		vector<cell_type> extant_species_vec(extant_species.begin(), extant_species.end());
		debug("Local extant species: %d", extant_species.size());
		vector<cell_type> global_extant_species_vec = mxx::gatherv(extant_species_vec, 0);
		if (me == 0) {
			set<cell_type> global_extant_species(global_extant_species_vec.begin(), global_extant_species_vec.end());
			debug("Global extant species: %d", global_extant_species.size());
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
			debug("Pruned, root time = %d", tree->time);

			outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".tree";
			ofstream fout(outfile.c_str());
			fout << toString_final(tree, timescale) << endl;
			delete_tree(tree);
			fout.close();
		}
		out_end_time = std::chrono::system_clock::now();
		std::chrono::duration<double> out_time = out_end_time - out_start_time;
		println("Output run time = %fs", out_time.count());
		} // end #pragma omp single
	}

} /* #pragma omp parallel */
	end_time = std::chrono::system_clock::now();
	std::chrono::duration<double> total_time = end_time - start_time;
	println("Total run time = %fs", total_time.count());
	println("Wrote results to file %s", outfile.c_str());
	fflush(stdout);

	GA_Destroy(ga_land_grid);

	GA_Terminate();
	MPI_Finalize();
}
