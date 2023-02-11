#include "petscmat.h"  
#include "sfmt.h"
#include "sfmt.cpp"
#include "userintf.cpp"
#include <chrono>
#include <random>

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
    PetscInt       local_rows, local_cols;
    PetscInt       global_rows, global_cols;
    PetscInt       local_start_col, local_end_col;
    PetscInt       local_start_row, local_end_row;
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

    // RNG
    CRandomSFMT1 RanGen((int)time(NULL)); // Agner Combined generator
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution(size*size, specrate);

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

    // setup grid for average across reps
    ierr = MatCreateDense(MPI_COMM_WORLD, 100, 25, PETSC_DECIDE, PETSC_DECIDE, NULL, &mean_grid); CHKERRQ(ierr);
    ierr = MatCreateDense(MPI_COMM_WORLD, 100, 25, PETSC_DECIDE, PETSC_DECIDE, NULL, &rep_grid); CHKERRQ(ierr);
    //ierr = MatZeroEntries(mean_grid); CHKERRQ(ierr);

    for(int rep = 0; rep < nreps; rep++) {
        rep_start_time = std::chrono::system_clock::now();

        // TODO: fill the rep_grid with 1s
        ierr = MatGetLocalSize(rep_grid, &local_rows, &local_cols); CHKERRQ(ierr);
        ierr = MatGetSize(rep_grid, &global_rows, &global_cols);CHKERRQ(ierr);
        ierr = MatGetOwnershipRangeColumn(rep_grid, &local_start_col, &local_end_col);CHKERRQ(ierr);
        ierr = MatGetOwnershipRange(rep_grid, &local_start_row, &local_end_row);CHKERRQ(ierr);
        ierr = MatGetOwnershipIS(rep_grid, &local_row_IS, &local_col_IS);CHKERRQ(ierr);
        ierr = ISGetSize(local_row_IS, &local_row_count);CHKERRQ(ierr);
        ierr = ISGetSize(local_col_IS, &local_col_count);CHKERRQ(ierr);
        println("Proc %d: Local rows = %d, Local cols = %d, Global rows = %d, Global cols = %d, Local start col = %d, Local end col = %d, Local start row = %d, Local end row = %d, local col count = %d, local row count = %d", comm_rank, local_rows, local_cols, global_rows, global_cols, local_start_col, local_end_col, local_start_row, local_end_row, local_col_count, local_row_count);

        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        //println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);
    }

    ierr = PetscFinalize();
    return ierr;
}

