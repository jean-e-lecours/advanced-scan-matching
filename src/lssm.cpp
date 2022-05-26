#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "../incl/psr.hpp"
#include "../incl/lssm.hpp"


std::vector<double> lssm::make_g_vector(std::vector<Correlation>& corrs, std::vector<Point2D>& scan_pts, bool print){
    std::vector<double> g = {0,0,0,0};
    double g_const = 0;
    int data_size = corrs.size();
    for (int i = 0; i < data_size; i++){
            g_const = (corrs[i].mcor1.x * corrs[i].norm[0] + corrs[i].mcor1.y * corrs[i].norm[1]);
            //g = [nx⋅(nx⋅qx + ny⋅qy)  ny⋅(nx⋅qx + ny⋅qy)  (nx⋅px + ny⋅py)⋅(nx⋅qx + ny⋅qy)  (-nx⋅py + ny⋅px)⋅(nx⋅qx + ny⋅qy)]
            g[0] += corrs[i].norm[0] * g_const;
            g[1] += corrs[i].norm[1] * g_const;
            g[2] += (corrs[i].norm[0] * scan_pts[i].x) + (corrs[i].norm[1] * scan_pts[i].y) * g_const;
            g[3] += (-corrs[i].norm[0] * scan_pts[i].y) + (corrs[i].norm[1] * scan_pts[i].x) * g_const;
    }
    if (print){
        //Print vector so that it can easily be imported into matlab if need be
        std::cout << "g=[" << g[0] << ',' << g[1] << ',' << g[2] << ',' << g[3] << "]\n\n";
    }
    return g;
}

std::vector<double> lssm::make_M_matrix(std::vector<Correlation>& corrs, std::vector<Point2D>& scan_pts, bool print){
    std::vector<double> M = {0,0,0,0,
                                0,0,0,0,
                                0,0,0,0,
                                0,0,0,0};
    int data_size = corrs.size();
    for (int i = 0; i < data_size; i++){

/*⎡          2                                                                                               ⎤
⎢    0   nx                1 nx⋅ny              2 nx⋅(nx⋅px + ny⋅py)             3  nx⋅(-nx⋅py + ny⋅px)       ⎥
⎢                                                                                                          ⎥
⎢                               2                                                                          ⎥
⎢     4 nx⋅ny               5 ny               6  ny⋅(nx⋅px + ny⋅py)           7   ny⋅(-nx⋅py + ny⋅px)       ⎥
⎢                                                                                                          ⎥
⎢                                                              2                                          ⎥
⎢8nx⋅(nx⋅px + ny⋅py)   9 ny⋅(nx⋅px + ny⋅py)      10(nx⋅px + ny⋅py)         11(nx⋅px + ny⋅py)⋅(-nx⋅py + ny⋅px)  ⎥ 
⎢                                                                                                         ⎥
⎢                                                                                               2         ⎥
⎣12nx⋅(-nx⋅py + ny⋅px) 13 ny⋅(-nx⋅py + ny⋅px)  14-(nx⋅px + ny⋅py)⋅(nx⋅py - ny⋅px)  15(nx⋅py - ny⋅px)           ⎦ */

        M[0] += corrs[i].norm[0] * corrs[i].norm[0];
        M[5] += corrs[i].norm[1] * corrs[i].norm[1];
        M[10] += (corrs[i].norm[0] * scan_pts[i].x + corrs[i].norm[1] * scan_pts[i].y) * (corrs[i].norm[0] * scan_pts[i].x + corrs[i].norm[1] * scan_pts[i].y);
        M[15] += (corrs[i].norm[0] * scan_pts[i].y - corrs[i].norm[1] * scan_pts[i].x) * (corrs[i].norm[0] * scan_pts[i].y - corrs[i].norm[1] * scan_pts[i].x);

        M[1] += corrs[i].norm[0] * corrs[i].norm[1];                                                        M[4] = M[1];
        M[2] += corrs[i].norm[0] * (corrs[i].norm[0] * scan_pts[i].x + corrs[i].norm[1] * scan_pts[i].y);   M[8] = M[2];
        M[3] += corrs[i].norm[0] * (-corrs[i].norm[0] * scan_pts[i].y + corrs[i].norm[1] * scan_pts[i].x);  M[12] = M[3];

        M[6] += corrs[i].norm[1] * (corrs[i].norm[0] * scan_pts[i].x + corrs[i].norm[1] * scan_pts[i].y);   M[9] = M[6];
        M[7] += corrs[i].norm[1] * (-corrs[i].norm[0] * scan_pts[i].y + corrs[i].norm[1] * scan_pts[i].x);  M[13] = M[7];

        M[11] += (corrs[i].norm[0] * scan_pts[i].x + corrs[i].norm[1] * scan_pts[i].y) * (-corrs[i].norm[0] * scan_pts[i].y + corrs[i].norm[1] * scan_pts[i].x); M[14] = M[11];

    }
    //Set values that reappear
    M[4] = M[1]; M[8] = M[2]; M[12] = M[3];
    M[9] = M[6]; M[13] = M[7];
    M[14] = M[11];

    //Print matrix so it can be easily imported into matlab if need be;
    if (print){
        std::cout << std::setprecision(30);
        std::cout << "M=[";
        for (int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                std::cout << M[4*(i) + j];
                if (j+1 < 6){
                    std::cout << ',';
                }
            }
            if (i+1 < 4){
                std::cout << ';';
            }
        }
        std::cout << "]\n\n";
    }

    return M;
}

std::vector<double> lssm::solve_system(std::vector<double>& g, std::vector<double>& G, bool print){
    std::vector<double> x = {0,0,0,0};
    double mult = 0;
    //This will make the G matrix triangular using Gaussian Elimination
    //TODO: swap lines in case of bad pivot point
    for (int j = 0; j < 3; j++){
        for (int i = j+1; i < 4; i++){
            mult = G[4*i+j]/G[4*j+j];
            for (int k = 0; k < 4; k++){
                G[4*i+k] = G[4*i+k] - mult * G[6*j+k];
            }
            g[i] = g[i] - (mult * g[j]);
        }
    }
    //Solve the now triangular system
    double sum = 0;
    x[3] = g[3]/G[15];
    for (int i = 2; i >= 0; i--){
        sum = 0;
        for (int j = 3; j >= i+1; j--){
            sum += G[4*i+j] * x[j];
        }
        x[i] = (g[i] - sum)/G[4*i+i];
    }
    //Print the solution vector to double check if need be
    if (print){
        std::cout << "Solution vector is: ";
        for (int i = 0; i < 4; i++){
            std::cout << x[i] << '\t';
        }
        std::cout << "\n\n";
    }
    return x;
}

void lssm::update_transform(Transform2D& transform, std::vector<double>& x){
    transform.rot_mat[0] = x[2];
    transform.rot_mat[1] = -x[3];
    transform.rot_mat[2] = x[3];
    transform.rot_mat[3] = x[2];
    transform.trans_vec[0] = x[0];
    transform.trans_vec[1] = x[1];

    transform.z_rot = std::atan2(x[3],x[2]);
}

