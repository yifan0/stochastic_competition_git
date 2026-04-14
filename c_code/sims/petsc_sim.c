#include "petscmat.h"  
#include "sfmt.h"
#include "sfmt.cpp"
#include "userintf.cpp"
#include <chrono>
#include <random>
#include <set>
#include <map>

typedef PetscScalar cell_type;

static char help[] = "Simulate plant spread on grid\n\t-gridsize <size> : number of cells along an edge of the grid\n\t-nreps <reps> : number of repetitions of the simulation\n\t-specrate <rate> : probability of speciaiton event\n\t-invrate <rate> : probability of invasion event\n\t-outfile <file> : output file name\n\n";

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

int main(int argc,char **argv) {
    PetscErrorCode ierr;
    char           outfile_name[PETSC_MAX_PATH_LEN] = "sim.log";
    PetscInt       size = 100;
    PetscInt       nreps = 1;
    PetscScalar    p = 0.1;
    PetscScalar    mutsize = 0.1;
    PetscScalar    specrate = 0.0001;
    PetscScalar    invrate = 0.2;
    int            timescale;
    int            nsteps = 100;
    int            endtime;
    Mat            matrix = NULL;
    Mat            mean_grid = NULL;
    Mat            rep_grid = NULL;
    PetscScalar*   rep_grid_data = NULL;
    PetscInt       local_rows, local_cols;
    PetscInt       global_rows, global_cols;
    PetscInt       local_start_col, local_end_col;
    PetscInt       mean_start_col, mean_end_col;
    PetscInt       local_start_row, local_end_row;
    PetscInt       mean_start_row, mean_end_row;
    IS             local_row_IS, local_col_IS;
    PetscInt       local_row_count, local_col_count;
    PetscMPIInt    comm_rank, comm_size;
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time, rep_start_time, rep_end_time, out_start_time, out_end_time;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

    // get options
    ierr = PetscOptionsGetInt(NULL,NULL, "-gridsize", &size, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL, "-nreps", &nreps, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(NULL,NULL, "-specrate", &specrate, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(NULL,NULL, "-invrate", &invrate, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL, "-outfile", outfile_name, sizeof(outfile_name), NULL); CHKERRQ(ierr);
    timescale = 100*size/p;
    endtime = timescale/nsteps;
    if(size <= 0 || nreps <= 0 || specrate > 1 || specrate < 0 || invrate > 1 || invrate < 0) {
        fprintf(stderr, help);
        return 1;
    }
    MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&comm_rank);

    // print out options
    if(comm_rank == 0) {
        println("Inputs:");
        println("\trepetitions = %d", nreps);
        println("\tgrid size = %d%s%d", size, "x", size);
        println("\tindividuals per patch = %f", 1/p);
        println("\tmutation size = %f", mutsize);
        println("\tspeciation rate = %f", specrate);
        println("\tinvasion rate = %f", invrate);
        println("\ttimescale = %d", timescale);
        println("\tnsteps = %d", nsteps);
        println("\tend time = %d", endtime);
        println("");
        fflush(stdout);
    }

    start_time = std::chrono::system_clock::now();

    // setup grid for average across reps and for each rep
    ierr = MatCreateDense(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, size, size, NULL, &mean_grid); CHKERRQ(ierr);
    ierr = MatCreateDense(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, size, size, NULL, &rep_grid); CHKERRQ(ierr);
    ierr = MatGetOwnershipRangeColumn(rep_grid, &local_start_col, &local_end_col);CHKERRQ(ierr);
    ierr = MatGetOwnershipRangeColumn(mean_grid, &mean_start_col, &mean_end_col);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(rep_grid, &local_start_row, &local_end_row);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(mean_grid, &mean_start_row, &mean_end_row);CHKERRQ(ierr);

    // RNG
    CRandomSFMT1 RanGen((int)time(NULL)); // Agner Combined generator
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution((local_end_row-local_start_row)*(local_end_col-local_start_col), specrate); // gives local speciation event count

    // invasion rule variables
    cell_type neighborhood[8];
    cell_type inv[8];
    cell_type inv_sum = 0;
    int inv_index = 0;

    for(int rep = 0; rep < nreps; rep++) {
        rep_start_time = std::chrono::system_clock::now();

        //ierr = MatDenseGetArray(rep_grid, &rep_grid_data); CHKERRQ(ierr);
        //ierr = MatGetOwnershipIS(rep_grid, &local_row_IS, &local_col_IS);CHKERRQ(ierr);
        //ierr = ISGetSize(local_row_IS, &local_row_count);CHKERRQ(ierr);
        //ierr = ISGetSize(local_col_IS, &local_col_count);CHKERRQ(ierr);
        //println("Proc %d: Local rows = %d, Local cols = %d, Global rows = %d, Global cols = %d, Local start col = %d, Local end col = %d, Local start row = %d, Local end row = %d, local col count = %d, local row count = %d", comm_rank, local_rows, local_cols, global_rows, global_cols, local_start_col, local_end_col, local_start_row, local_end_row, local_col_count, local_row_count);

        // fill rep_grid with 1s
        ierr = MatAssemblyBegin(rep_grid, MAT_FLUSH_ASSEMBLY);
        ierr = MatAssemblyEnd(rep_grid, MAT_FLUSH_ASSEMBLY);
        for(PetscInt i = local_start_row; i < local_end_row; i++) {
            for(PetscInt j = local_start_col; j < local_end_col; j++) {
                ierr = MatSetValue(rep_grid, i, j, 1, INSERT_VALUES); CHKERRQ(ierr); // TODO: change to MatSetValues
            }
        }

        ierr = MatAssemblyBegin(rep_grid, MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(rep_grid, MAT_FINAL_ASSEMBLY);

        // iteration loop
        for(int step = 0; step < timescale; step++) {
            // speciation event
            int speciation_event_count = speciation_distribution(generator);
            std::set<int> spec_events;
            while(spec_events.size() < speciation_event_count) {
                int i = RanGen.IRandom(local_start_row, local_end_row-1);
                int j = RanGen.IRandom(local_start_col, local_end_col-1);
                int index = i*(local_end_row-local_start_row) + j;
                if(!spec_events.count(index)) {
                    spec_events.insert(index);
                    int down = RanGen.IRandom(0, 1);
                    float ratio = (1+RanGen.Random()*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        PetscScalar val = 1;
                        ierr = MatGetValue(rep_grid, i, j, &val); CHKERRQ(ierr);
                        val = val * ratio;
                        ierr = MatSetValue(rep_grid, i, j, val, INSERT_VALUES); CHKERRQ(ierr);
                        ierr = MatAssemblyBegin(rep_grid, MAT_FINAL_ASSEMBLY);
                        ierr = MatAssemblyEnd(rep_grid, MAT_FINAL_ASSEMBLY);

                    }
                }
            }

            // invasion event
            std::map<int, cell_type> invasion_event_archive;
            for(PetscInt i = local_start_row; i < local_end_row; i++) {
                for(PetscInt j = local_start_col; j < local_end_col; j++) {
                    double randval = RanGen.Random();
                    if(randval < invrate) {
                        PetscScalar value;
                        ierr = MatGetValue(rep_grid, i, j, &value); CHKERRQ(ierr);
                        inv_sum = 0;
                        inv_index = 0;
                        for(int x = -1; x < 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                if((x != 0 || y != 0) && i+x >= 0 && i+x < size && j+y >= 0 && j+y < size) {
                                    if(invasion_event_archive.find( (i+x)*size+j+y ) != invasion_event_archive.end() )
                                        neighborhood[inv_index] = invasion_event_archive[(i+x)*size+j+y];
                                    else {
                                        PetscScalar nvalue;
                                        ierr = MatGetValue(rep_grid, i+x, j+x, &nvalue); CHKERRQ(ierr);
                                        neighborhood[inv_index] = nvalue;
                                    }
                                    inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+value*(1-p));
                                    inv_sum += inv[inv_index];
                                    inv_index++;
                                }
                            }
                        }
                    }
                }
            }
        }

        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        //println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);
    }

    ierr = PetscFinalize();
    return ierr;
}

