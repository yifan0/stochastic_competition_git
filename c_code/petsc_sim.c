#include "petscmat.h"  

static char help[] = "Simulate plant spread on grid\n\
 -gridsize <size> : number of cells along an edge of the grid\n\
 -nreps <reps> : number of repetitions of the simulation\n\
 -outfile <file> : output file name\n\n";

int main(int argc,char **argv) {
    PetscErrorCode ierr;
    char           outfile_name[PETSC_MAX_PATH_LEN] = "sim.log";
    PetscInt       gridsize = 100;
    PetscInt       nreps = 1;
    Mat            matrix = NULL;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

    // get options
    ierr = PetscOptionsGetInt(NULL,NULL, "-gridsize", &gridsize, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL, "-nreps", &nreps, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL, "-outfile", outfile_name, sizeof(outfile_name), NULL); CHKERRQ(ierr);

    // setup matrix
    ierr = MatCreateDense(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, gridsize, gridsize, NULL, &matrix); // TODO: figure out how to fill with 1s
    ierr = MatSetValues(matrix, );

    ierr = PetscFinalize();
    return ierr;
}

