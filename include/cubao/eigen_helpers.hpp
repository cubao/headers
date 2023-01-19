#ifndef CUBAO_EIGEN_HELPERS_HPP
#define CUBAO_EIGEN_HELPERS_HPP

// should sync
// - https://github.com/cubao/polyline-ruler/blob/master/src/eigen_helpers.hpp
// - https://github.com/cubao/headers/tree/main/include/cubao/eigen_helpers.hpp

// https://github.com/microsoft/vscode-cpptools/issues/9692
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <Eigen/Core>
#include <numeric>

namespace cubao
{
using RowVectors = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RowVectorsNx3 = RowVectors;
using RowVectorsNx2 = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;

// https://stackoverflow.com/a/73799908/5089147
inline Eigen::VectorXd arange(double low, double high, double step,
                              bool with_last = false)
{
    int N = static_cast<int>(std::floor((high - low) / step) + 1);
    Eigen::VectorXi V(N);
    std::iota(V.data(), V.data() + N, 0);
    Eigen::VectorXd ret = V.cast<double>();
    ret.array() *= step;
    ret.array() += low;
    if (with_last && ret[N - 1] != high) {
        ret.conservativeResize(N + 1);
        ret[N] = high;
    }
    return ret;
}

inline Eigen::VectorXi
indexes2mask(const Eigen::Ref<const Eigen::VectorXi> &indexes, int N)
{
    Eigen::VectorXi mask(N);
    mask.setZero();
    for (int c = 0, C = indexes.size(); c < C; ++c) {
        mask[indexes[c]] = 1;
    }
    return mask;
}

inline Eigen::VectorXi
mask2indexes(const Eigen::Ref<const Eigen::VectorXi> &mask)
{
    Eigen::VectorXi indexes(mask.sum());
    for (int i = 0, j = 0, N = mask.size(); i < N; ++i) {
        if (mask[i]) {
            indexes[j++] = i;
        }
    }
    return indexes;
}

inline RowVectors select_by_mask(const Eigen::Ref<const RowVectors> &coords,
                                 const Eigen::Ref<const Eigen::VectorXi> &mask)
{
    RowVectors ret(mask.sum(), coords.cols());
    int N = mask.size();
    for (int i = 0, k = 0; i < N; ++i) {
        if (mask[i]) {
            ret.row(k++) = coords.row(i);
        }
    }
    return ret;
}

inline RowVectors to_Nx3(const Eigen::Ref<const RowVectorsNx2> &coords)
{
    RowVectors _coords(coords.rows(), 3);
    _coords.leftCols(2) = coords;
    _coords.col(2).setZero();
    return _coords;
}

} // namespace cubao

#endif
