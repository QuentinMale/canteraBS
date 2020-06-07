#include "cantera/base/ct_defs.h"
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

namespace Cantera {
    typedef Eigen::SparseMatrix<double> SparseMat_fp;
    typedef Eigen::Triplet<double> Triplet_fp;
}
