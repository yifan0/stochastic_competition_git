#include <stdlib.h>
#include <argp.h>
#include <stdbool.h>

const char *argp_program_version = "sim 1.0";
static char doc[] = "Simulation of speciation events";
static char args_doc[] = "LIST OF ARGS"; // TODO: add list of args
static struct argp_option options[] = {
    { "tree", 't', 0, 0, "Generate tree representation of speciation events."},
    { "fout", 'f', "FILE", 0, "Output file for results."},
    { "results-file", 'r', 0, 0, "Save comparison results to a file."},
    { "parallel", 'p', 0, 0, "Run in parallel using MPI."},
    { "dims", 'd', "DIMS", 0, "Dimensions of the grid."},
    { "size", 's', "SIZE", 0, "Size of the grid."},
    { "reps", 'n', "REPS", 0, "Number of repetitions."},
    { "specrate", 'c', "SPEC_RATE", 0, "Speciation rate."},
    { "mutsize", 'm', "MUT_SIZE", 0, "Maximum change in mutation event."},
    { "diff",  'u', 0, 0, "Run difference statistic."},
    { "div",   'v', 0, 0, "Run diversity statistic."},
    { "width", 'w', 0, 0, "Run width statistic."},
    { 0 }
};

struct arguments {
    bool tree;
    char* fout;
    bool results_file;
    bool parallel;
    int dims;
    int size;
    int reps;
    double specrate;
    double mutsize;
    bool diff;
    bool div;
    bool width;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = (struct arguments *)state->input;
    switch (key) {
    case 't': arguments->tree = true; break;
    case 'f': arguments->fout = arg; break;
    case 'r': arguments->results_file = true; break;
    case 'p': arguments->parallel = true; break;
    case 'd': arguments->dims = atoi(arg); break;
    case 's': arguments->size = atoi(arg); break;
    case 'n': arguments->reps = atoi(arg); break;
    case 'c': arguments->specrate = atof(arg); break;
    case 'm': arguments->mutsize = atof(arg); break;
    case 'u': arguments->diff = true; break;
    case 'v': arguments->div = true; break;
    case 'w': arguments->width = true; break;
    case ARGP_KEY_ARG: return 0;
    default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

int main(int argc, char *argv[])
{
    struct arguments arguments;

    // default parameter values
    arguments.tree = false;
    arguments.fout = "out.csv";
    arguments.results_file = false;
    arguments.parallel = false;
    arguments.dims = 2;
    arguments.size = 100;
    arguments.reps = 10;
    arguments.specrate = 0.1;
    arguments.mutsize = 0.01;
    arguments.diff = false;
    arguments.div = false;
    arguments.width = false;

    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    fprintf(stdout, "tree = %s\n", arguments.tree ? "true" : "false");
    fprintf(stdout, "fout = %s\n", arguments.fout);
    fprintf(stdout, "results_file = %s\n", arguments.results_file ? "true" : "false");
    fprintf(stdout, "parallel = %s\n", arguments.parallel ? "true" : "false");
    fprintf(stdout, "dims = %d\n", arguments.dims);
    fprintf(stdout, "size = %d\n", arguments.size);
    fprintf(stdout, "reps = %d\n", arguments.reps);
    fprintf(stdout, "specrate = %f\n", arguments.specrate);
    fprintf(stdout, "mutsize = %f\n", arguments.mutsize);
    fprintf(stdout, "diff = %s\n", arguments.diff ? "true" : "false");
    fprintf(stdout, "div = %s\n", arguments.div ? "true" : "false");
    fprintf(stdout, "width = %s\n", arguments.width ? "true" : "false");

    // TODO: add the actual simulation code from other files

}
