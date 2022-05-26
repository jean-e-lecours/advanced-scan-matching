#ifndef LSSM_HPP
#define LSSM_HPP

#include <vector>
#include <iostream>

#include "psr.hpp"

namespace lssm {
    std::vector<double> make_g_vector(std::vector<Correlation>& corrs, std::vector<Point2D>& scan_pts, bool print);

    std::vector<double> make_M_matrix(std::vector<Correlation>& corrs, std::vector<Point2D>& scan_pts, bool print);
     
    std::vector<double> solve_system(std::vector<double>& g, std::vector<double>& G, bool print);

    void update_transform(Transform2D& transform, std::vector<double>& x);
}

#endif