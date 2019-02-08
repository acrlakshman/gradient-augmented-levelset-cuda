/*
 * Data types defined for use in the serial and MPI codes
 * These data types could not be used in the CUDA code.
 */

namespace galsfunctions
{
    typedef std::tuple<double,double> psinode;
    typedef std::vector<std::vector<double>> gridarray;
    typedef std::vector<double> vectorarray;
    typedef std::vector<std::vector<unsigned int>> intgridarray;
    typedef std::vector<unsigned int> intvectorarray;
}
