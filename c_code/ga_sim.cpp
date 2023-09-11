#include <cstdio>
#include <random>	// rand()
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
#include <mpi.h>
#include <iostream>	   // std::cout
#include <fstream>	   // std::ofstream
#include <cxxopts.hpp> // to handle cmdline args
#include "ga++.h"
#include "tree.h"
#include <mxx/reduction.hpp>
using namespace std;

#define println(...) { if(me == 0) { printf(__VA_ARGS__); printf("\n"); } }
#define debug(...) { if(true) { printf("Rank %d: ", me); printf(__VA_ARGS__); printf("\n"); fflush(stdout); } }
#define print(...) { printf(__VA_ARGS__); }

#define NDIM 2
#define GHOSTS 2

typedef std::tuple<int, int, cell_type> cell_update;

int main(int argc, char *argv[])
{
	int nrep = 10;
	int size = 100;
	double p = 0.1;
	double mutsize = 0.1;
	double specrate = 0.0001; // 1.0e-4
	double invrate = 0.2;
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

	if (result.count("help"))
	{
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
	for (size_t i = 0; i < NDIM; i++)
	{
		dims[i] = size;
		ghost_width[i] = GHOSTS;
	}
	char land_grid_name[] = "land grid";
	int ga_land_grid = NGA_Create_ghosts(C_DBL, NDIM, dims, ghost_width, land_grid_name, NULL);
	if (ga_land_grid == 0)
	{
		char create_err[] = "Failure for NGA_Create_ghosts()";
		GA_Error(create_err, 1);
	}

	NGA_Distribution(ga_land_grid, me, lo, hi);
	if (lo[0] < 0)
	{
		debug("Owns no elements");
		return 1;
	}
	NGA_Access(ga_land_grid, lo, hi, &land_grid_ptr, grid_ld);

	println("Inputs:");
	println("\trepetitions = %d", nrep);
	println("\tsize = %d%s%d", size, "x", size);
	println("\tindividuals per patch = %f", 1 / p);
	println("\tmutation size = %f", mutsize);
	println("\tspeciation rate = %.2e", specrate);
	println("\tinvasion rate = %f", invrate);
	println("\ttimescale = %d", timescale);
	println("\tnsteps = %d", nsteps);
	println("\tend time = %d", endtime);
	println("\tprocesses = %d", num_procs);
	println("");
	fflush(stdout);

	// Random number generation
	CRandomSFMT1 RanGen(time(0) + me * 10); // Agner Combined generator
	int local_rows = (hi[0] - lo[0] + 1);
	int local_cols = (hi[1] - lo[1] + 1);
	int local_area = local_cols * local_rows;
	double land_mask_data[local_area];
	double* land_mask[local_rows];
	for (size_t i = 0; i < local_rows; i++)
		land_mask[i] = land_mask_data + i * local_cols;
	// distribution for speciation events
	std::default_random_engine generator;
	std::binomial_distribution<int> speciation_distribution(local_area, specrate);

	for (int rep = 0; rep < nrep; rep++)
	{
		rep_start_time = std::chrono::system_clock::now();

		double one = 1;
		GA_Fill(ga_land_grid, &one);
		GA_Update_ghosts(ga_land_grid);
		NGA_Access_ghosts(ga_land_grid, ghost_dims, &ghost_grid_ptr, ghost_grid_ld);
		if (!ghost_grid_ptr)
		{
			char grid_err[] = "NULL pointer for ghost grid.";
			GA_Error(grid_err, 1);
		}

		// zero out mask except local edges
		memset(land_mask_data, 0, sizeof(*land_mask_data));
		for (size_t i = 0; i < local_rows; i++)
		{
			land_mask[i][0] = 1;
			land_mask[i][local_cols - 1] = 1;
		}
		for (size_t i = 0; i < local_cols; i++)
		{
			land_mask[0][i] = 1;
			land_mask[local_rows - 1][i] = 1;
		}

		// invasion rule variables
		cell_type neighborhood[8];
		cell_type inv[8];
		cell_type inv_sum = 0;
		int inv_index = 0;

		// track speciation events
		vector<tuple<size_t, cell_type, cell_type>> speciation_events;

		int row, col;
		for (int step = 0; step < timescale; step++)
		{

			// speciation rule
			int speciation_event_count = speciation_distribution(generator);
			std::set<int> spec_events;
			while (spec_events.size() < speciation_event_count)
			{
				int index = rand() % local_area;
				if (!spec_events.count(index))
				{
					size_t i = index / local_cols; // row
					size_t j = index % local_cols; // col
					spec_events.insert(index);
					int down = rand() % 2;
					float ratio = (1 + RanGen.Random() * mutsize);
					if (down)
					{
						ratio = 1 / ratio;
					}
					float probsuccess = p * ratio / (p * (ratio - 1) + 1);
					if (RanGen.Random() <= probsuccess)
					{
						cell_type old_val = ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS];
						cell_type new_val = old_val * ratio;
						ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS] = new_val;
						speciation_events.push_back(make_tuple(step, old_val, new_val));
					}

					// set mask for [i][j] and surrounding cells unless on edge
					for (int x = -1; x <= 1; x++)
					{
						for (int y = -1; y <= 1; y++)
						{
							row = i + x;
							col = j + y;
							if (row >= local_rows - 1 || row <= 0)
								continue;
							if (col >= local_cols - 1 || col <= 0)
								continue;
							cell_type local_max = ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS];
							for (int xx = -1; xx <= 1; xx++)
							{
								for (int yy = -1; yy <= 1; yy++)
								{
									int neighbor_row = row + xx;
									int neighbor_col = col + yy;
									local_max = std::max(local_max, ghost_grid_ptr[(neighbor_row + GHOSTS) * ghost_grid_ld[0] + neighbor_col + GHOSTS]);
								}
							}
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols)
								land_mask[row][col] = local_max * p / (local_max * p + ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] * (1 - p));
						}
					}
				}
			}

			if (step % 10 == 0)
				GA_Update_ghosts(ga_land_grid);

			// invasion rule
			std::vector<cell_update> updates;
			for (size_t i = 0; i < local_rows; i++)
			{
				for (size_t j = 0; j < local_cols; j++)
				{
					if (land_mask[i][j] != 0)
					{
						double randval = RanGen.Random(); //random_float();
						if (randval < land_mask[i][j])
						{
							inv_sum = 0;
							inv_index = 0;
							for (int x = -1; x <= 1; x++)
							{
								for (int y = -1; y <= 1; y++)
								{
									if ((x != 0 || y != 0))
									{
										neighborhood[inv_index] = ghost_grid_ptr[(i + x + GHOSTS) * ghost_grid_ld[0] + j + y + GHOSTS];
										inv[inv_index] = p * neighborhood[inv_index] / (p * neighborhood[inv_index] + ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS] * (1 - p));
										inv_sum += inv[inv_index];
										inv_index++;
									}
								}
							}
							if (randval <= inv_sum / 8)
							{
								// Get random element with weighted probabilities
								double weighted_rand = RanGen.Random() * inv_sum;
								inv_index = 0;
								while (weighted_rand > inv[inv_index])
								{
									weighted_rand -= inv[inv_index];
									inv_index++;
								}
								if (neighborhood[inv_index] != ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS])
								{
									updates.push_back({i, j, neighborhood[inv_index]});
								}
							}
						}
					}
				}
			}
			
			for (const auto &[i, j, val] : updates)
			{
				ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS] = val;
				for (int x = -1; x <= 1; x++)
				{
					for (int y = -1; y <= 1; y++)
					{
						row = i + x;
						col = j + y;
						if (row >= local_rows - 1 || row <= 0)
							continue;
						if (col >= local_cols - 1 || col <= 0)
							continue;
						cell_type local_max = 0;
						for (int xx = -1; xx <= 1; xx++)
						{
							for (int yy = -1; yy <= 1; yy++)
							{
								row = i + x + xx;
								col = j + y + yy;
								if (row >= local_rows - 1 || row <= 0)
									continue;
								if (col >= local_cols - 1 || col <= 0)
									continue;
								if (ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] != val && ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] > local_max)
								{
									local_max = ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS];
								}
							}
						}
						row = i + x;
						col = j + y;
						if (local_max == 0)
						{
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols)
								land_mask[row][col] = 0;
						}
						else
						{
							if(row >= 0 && col >= 0 && row < local_rows && col < local_cols)
								land_mask[row][col] = local_max * p / (local_max * p + ghost_grid_ptr[(row + GHOSTS) * ghost_grid_ld[0] + col + GHOSTS] * (1 - p));
						}
					}
				}
			}

			// renormalize every nstep steps
			if (step % endtime == endtime - 1)
			{
				cell_type land_grid_mean = 0;
				for (size_t i = 0; i < local_rows; i++)
					for (size_t j = 0; j < local_cols; j++)
						land_grid_mean += ghost_grid_ptr[(i + GHOSTS) * ghost_grid_ld[0] + j + GHOSTS];
				land_grid_mean = land_grid_mean / (size * size);
				MPI_Allreduce(MPI_IN_PLACE, &land_grid_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				println("Global average at step %d = %f", step, land_grid_mean);
				// TODO: pick a useful threshold for when to normalize
				if (land_grid_mean > 100)
				{
					if (me == 0)
						println("Normalizing by %f", land_grid_mean);

					for (tuple<size_t, cell_type, cell_type> &event : speciation_events) {
						get<1>(event) = get<1>(event) / land_grid_mean;
						get<2>(event) = get<2>(event) / land_grid_mean;
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
		outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".csv";
		fp = fopen(outfile.c_str(), "w");
		GA_Print_csv_file(fp, ga_land_grid);
		fclose(fp);

		// gather speciation events
		vector<tuple<size_t, cell_type, cell_type>> global_speciation_events = mxx::gatherv(speciation_events, 0);
		// gather extant species set
		// TODO: make extant_species vector
		debug("Gathering extant species list");
		set<cell_type> extant_species(ghost_grid_ptr+GHOSTS*ghost_grid_ld[0]+GHOSTS, ghost_grid_ptr+(local_rows + GHOSTS) * ghost_grid_ld[0] + local_cols + GHOSTS);
		vector<cell_type> extant_species_vec(extant_species.begin(), extant_species.end());
		debug("Local extant species: %d", extant_species.size());
		vector<cell_type> global_extant_species_vec = mxx::gatherv(extant_species_vec, 0);
		if (me == 0)
		{
			set<cell_type> global_extant_species(global_extant_species_vec.begin(), global_extant_species_vec.end());
			debug("Global extant species: %d", global_extant_species.size());
			// TODO: take this out
			for(cell_type x : global_extant_species) {
				debug("%f", x);
			}
			map<cell_type, speciation_tree_node*> leaves;
			// sort global_speciation_events by timestep
			sort(global_speciation_events.begin(), global_speciation_events.end());
			speciation_tree_node *tree = new speciation_tree_node(get<1>(global_speciation_events[0]), 0, nullptr);
			leaves[get<1>(global_speciation_events[0])] = tree;
			debug("Creating tree");
			for (tuple<size_t, cell_type, cell_type> event : global_speciation_events)
			{
				size_t timestep = get<0>(event);
				cell_type old_val = get<1>(event);
				cell_type new_val = get<2>(event);
				// instead of finding the node by traversing the tree, just look it up in a hash map (for performance)
				// speciation_tree_node *parent_node = find_node(tree, old_val);
				speciation_tree_node *parent_node = leaves[old_val];
				if (!parent_node)
				{
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
			debug("Binary tree: %d", binaryTree(tree));
			debug("Pruning tree");
			// NOTE: the tree is binary before pruning
			tree = prune(tree, global_extant_species);
			debug("Pruned, root time = %d", tree->time);

			// TODO: solve the problem that the tree is not binary

			debug("Binary tree: %d", binaryTree(tree));

			// prune(tree, global_extant_species);
			// TODO: the pruning isn't working, correct it -- there are too many species in the tree
			// TODO: test the pruning

			outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".tree";
			ofstream fout(outfile.c_str());
			fout << toString_final(tree, timescale) << endl;
			delete_tree(tree);
			fout.close();
		}

		out_end_time = std::chrono::system_clock::now();
		std::chrono::duration<double> out_time = out_end_time - out_start_time;
		println("Output run time = %fs", out_time.count());
	}

	end_time = std::chrono::system_clock::now();
	std::chrono::duration<double> total_time = end_time - start_time;
	println("Total run time = %fs", total_time.count());
	println("Wrote results to file %s", outfile.c_str());
	fflush(stdout);

	GA_Destroy(ga_land_grid);

	GA_Terminate();
	MPI_Finalize();
}
