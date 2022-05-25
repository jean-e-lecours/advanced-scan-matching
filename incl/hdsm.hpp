#ifndef HDSM_HPP
#define HDSM_HPP
#include "psr.hpp"

//Function to find Hausdorf Distance of two sets (which have correlations).
//Not specifying a cutoff value will give the true Hausdorf Distance,
//specifying one will give the corresponding partial Hausdorf Distance.
double find_hd(std::vector<Correlation>& corrs, Transform2D& transform, int cutoff = 0);

#endif