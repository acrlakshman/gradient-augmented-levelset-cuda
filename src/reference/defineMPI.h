/*
 * Required for MPI implementation
 */

#ifndef _DEFINEMPI_H
#define _DEFINEMPI_H

namespace mpi {
    class context {
        int m_rank, m_size;
    public:
        context(int *argc, char **argv[]) : m_rank { -1 } {
            if (MPI_Init(argc, argv) == MPI_SUCCESS) {
                MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
                MPI_Comm_size(MPI_COMM_WORLD, &m_size);
            }
        }
        ~context() {
            if(m_rank >= 0) {
                MPI_Finalize();
            }
        }
        explicit operator bool() const {
            return m_rank >= 0;
        }
        int rank() const noexcept { return m_rank; }
        int size() const noexcept { return m_size; }
    };
}
#endif
