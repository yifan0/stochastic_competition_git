#ifndef SIM_GRID
#define SIM_GRID

template <typename T>
class SimRow {
    private:
        T* data;
        SimGrid* grid;
        size_t row;

    public:
        SimRow(T* data_, size_t row_, SimGrid* grid_) : data(data_), row(row_), grid(grid_) {}

        // Subscript operator
        T& operator[](size_t i) const {
            return data[i];
        }
        
        // Subscript operator
        T& operator[](size_t i) {
            if(i < grid->ghosts || i > grid->local_cols-grid->ghosts) {
                grid->updateGhost(i,j);
            }
            return data[i];
        }
};

template <typename T>
class SimGrid {
    friend class SimRow;
    private:
        T* data;             // pointer to data
        SimRow<T>* rows;     // pointer to beginning of each row of data, excluding ghosts
        size_t global_size;  // size of grid
        size_t local_rows;   // local size of grid
        size_t local_cols;   // local size of grid
        size_t ghosts;       // number of ghost cells in each direction
        int* proc_rows;     // number of rows on each process
        int* proc_cols;     // number of cols on each process
        int proc_grid_rows; // number of rows of processes
        int proc_grid_cols; // number of cols of processes
        int rank, nprocs;    // MPI rank and comm size
        int sqrt_nprocs;    // sqrt of the number of processes
        MPI_Comm cart_comm; // cartesian topology communicator

        /*
         * Update remote element at (i,j) where i is the relative row and j is the relative column.
         */
        void updateGhost(size_t i, size_t j) {
            // Find destination process
            int coord[2];
            MPI_Cart_coords(comm, rank, 2, coord);
            int send_to = rank;
            while(i < 0) { // while the row is negative, check if it fits in the process directly above
                coord[0] -= 1;
                MPI_Cart_rank(comm, coord, &send_to);
                i += proc_rows[send_to];
            }
            while(i > proc_rows[send_to]) { // check the process below
                coord[0] += 1;
                i -= proc_rows[send_to];
                MPI_Cart_rank(comm, coord, &send_to);
            }
            while(j < 0) { // check process left
                coord[1] -= 1;
                MPI_Cart_rank(comm, coord, &send_to);
                j += proc_cols[send_to];
            }
            while(j > proc_cols[send_to]) { // check process right
                coord[1] += 1;
                j -= proc_cols[send_to];
                MPI_Cart_rank(comm, coord, &send_to);
            }

            // Find offset at destination including ghost offset
            MPI_Aint offset = i*(proc_cols[send_to]+2)+j+1;

            // TODO: do a put to the remote proc
            // Send to remote proc
        }

    public:
        // Constructor
        SimGrid(size_t size_, size_t ghosts_) : global_size(size_), ghosts(ghosts_) {

            // MPI communication info
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

            // Calculate local dimensions
            sqrt_nprocs = sqrt(nprocs);
            if(sqrt_nprocs*sqrt_nprocs != nprocs) {
                if(nprocs < 4) proc_grid_rows = 1;
                else {
                    for(proc_grid_rows = sqrt_nprocs; proc_grid_rows > 0; proc_grid_rows--) {
                        if(nprocs % proc_grid_rows == 0) {
                            break;
                        }
                    }
                    proc_grid_cols = nprocs / proc_grid_rows;
                }
            }
            local_rows = global_size/proc_grid_rows;
            local_cols = global_size/proc_grid_cols;
            int rows_remainder = global_size % proc_grid_rows;
            int cols_remainder = global_size % proc_grid_cols;
            if(rank % proc_grid_cols < cols_remainder) { // column of process
                local_cols += 1;
            }
            if(rank / proc_grid_rows < rows_remainder) { // row of process
                local_rows += 1;
            }

            // Share information with other processes about dimensions
            proc_rows = new int[nprocs];
            proc_cols = new int[nprocs];
            MPI_Allgather(&local_rows, 1, MPI_INT, &proc_rows, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Allgather(&local_cols, 1, MPI_INT, &proc_cols, 1, MPI_INT, MPI_COMM_WORLD);

            // Allocate memory for the data array, including ghost cells
            data = new T[(local_rows + 2*ghosts) * (local_cols + 2*ghosts)];

            // Allocate memory for the SimRow array, excluding ghost cells
            rows = new SimRow<T>[local_rows];
            for (size_t i = 0; i < local_rows; ++i) {
                rows[i] = SimRow<T>(data + (i+ghosts)*(local_cols+2*ghosts) + ghosts, i, this);
            }
            
            // Create cartesian topology
            int dim[2] = {proc_grid_rows, proc_grid_cols};
            int periods[2] = {1,1};
            MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &cart_comm);

        }

        // Destructor
        ~SimGrid() {
            delete[] data;
            delete[] rows;
            delete proc_rows;
            delete proc_cols;
        }

        // TODO: make sure these operators work accurately to access the ghost elements

        // Subscript operator for rows
        SimRow<T>& operator[](size_t i) const {
            return rows[i];
        }

        // Subscript operator for rows
        SimRow<T>& operator[](size_t i) {
            if (i < ghosts || i > local_rows-ghosts) {
                updateGhost(i,j); // TODO: where does j come from?
            }

            return rows[i];
        }

        // Subscript operator for individual cells
        T& operator()(size_t i, size_t j) const {
            return rows[i][j];
        }

        // Subscript operator for individual cells
        T& operator()(size_t i, size_t j) {
            if(i < ghosts || i > local_rows-ghosts || j < ghosts || j > local_cols-ghosts) {
                updateGhost(i,j);
            }
            return rows[i][j];
        }

};

#endif
