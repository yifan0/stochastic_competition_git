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
            if(i < grid->ghosts || i > grid->size-grid->ghosts) {
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
        size_t size;         // size of grid
        size_t ghosts;       // number of ghost cells in each direction
        int rank, nprocs;    // MPI rank and comm size

        void updateGhost(size_t i, size_t j) {
            // TODO: use MPI communication to update the ghost
        }

    public:
        // Constructor
        SimGrid(size_t size_, size_t ghosts_) : size(size_), ghosts(ghosts_) {

            // MPI communication info
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

            // Allocate memory for the data array, including ghost cells
            data = new T[(size + 2*ghosts) * (size + 2*ghosts)];

            // Allocate memory for the SimRow array, excluding ghost cells
            // TODO: adjust for MPI
            rows = new SimRow<T>[size];
            for (size_t i = 0; i < size; ++i) {
                rows[i] = SimRow<T>(data + (i+ghosts)*(size+2*ghosts) + ghosts, i, this);
            }
        }

        // Destructor
        ~SimGrid() {
            delete[] data;
            delete[] rows;
        }

        // TODO: make sure these operators work accurately to access the ghost elements

        // Subscript operator for rows
        SimRow<T>& operator[](size_t i) const {
            return rows[i];
        }
        
        // Subscript operator for rows
        SimRow<T>& operator[](size_t i) {
            if (i < ghosts || i > size-ghosts) {
                updateGhost(i,j);
            }

            return rows[i];
        }

        // Subscript operator for individual cells
        T& operator()(size_t i, size_t j) const {
            return rows[i][j];
        }

        // Subscript operator for individual cells
        T& operator()(size_t i, size_t j) {
            if(i < ghosts || i > size-ghosts || j < ghosts || j > size-ghosts) {
                updateGhost(i,j);
            }
            return rows[i][j];
        }

};

#endif
