#include <cmath>
#include <complex>
#include <iostream>
#include <ctime>
#include <unistd.h>
#include <vector>
#include <math.h>
#include <iomanip>

#include "incl/kdt.hpp"
#include "incl/lssm.hpp"
#include "incl/psr.hpp"
#include "incl/test.hpp"

#include "incl/asm.hpp"
#include "incl/hdsm.hpp"

enum alg {LSSM = 0, HDSM = 1, ASM = 2};

int main(int, char**) {

    bool trans = false;
    double map_res = 0.025;

    char algorithm = ASM;
    int max_guesses = 10;
    double threshold = 0.001;

    TextData scan_data;
    if (!scan_data.read_from_file("../dat/test_11.txt")){
        scan_data.read_from_file("dat/test_11.txt");
    }
    //commentddddd
    TextData map_data;
    if (!map_data.read_from_file("../dat/test_map.txt")){
        map_data.read_from_file("dat/test_map.txt");
    }

    Transform2D las_transform(0.25,0.25,0);
    Transform2D map_transform(0,0,0);

    clock_t start, end;

    start = clock();

    std::vector<Point2D> las_vec = map_scan_points(las_transform, scan_data.dat_vec, 2*M_PI/360);
    std::vector<Point2D> pure_las_vec = las_vec;

    std::vector<Point2D> map = make_map(map_data.dat_vec, 100, 100, map_res);

    map_transform.transform(map);

    KdTree map_kdt(map);

    Plot2D plot(true, "temp2d_");
    plot.add_data(map);
    plot.add_data(las_vec);

    Transform2D guess_transf(0,0,0);
    Transform2D total_transf(0,0,0);
    Transform2D basic_transf(0,0,0);
    std::vector<Correlation> correlations;

    int guesses = 0;

    do{
        correlations.clear();
        double corr_mean = 0;
        double corr_stdev = 0;

        for (int i = 0; i < las_vec.size(); i++){
            Correlation corr(map_kdt, las_vec[i], THOROUGH, total_transf);
            corr_mean += corr.corrected_value;
            correlations.push_back(corr);
        }
        
        corr_mean = corr_mean/correlations.size();

        for (int i = 0; i < correlations.size(); i++){
            corr_stdev += std::abs(correlations[i].corrected_value - corr_mean);
        }
        corr_stdev = sqrt(corr_stdev/correlations.size());

        std::vector<double> g = make_g_vector(correlations, corr_stdev, 0.1, false);
        std::vector<double> G = make_G_matrix(correlations, corr_stdev, 0.1, false);
        std::vector<double> x = solve_system(g, G, false);

        guess_transf.update_transform(x);

        /* double angle = std::atan2(guess_transf.rot_mat[2],guess_transf.rot_mat[0]);
        Transform2D rigid_guess(guess_transf.trans_vec[0], guess_transf.trans_vec[1], angle);

        guess_transf = rigid_guess; */

        //Transforms stay affine for now
        total_transf.add_transform(guess_transf);
        
        guesses++; 
        
    } while(guess_transf.is_significant(threshold) && guesses < max_guesses);
    plot.add_data(las_vec,total_transf);

    //Evaluate the distance between data and map points, remember the ones that are "converged"
    std::vector<int> sig_ind;
    for (int i = 0; i < correlations.size(); i++){

        if (correlations[i].get_distance() < map_res){
            sig_ind.push_back(i);
        }
    }
    //Project affine transform into SO2
    double angle = std::atan2(total_transf.rot_mat[2],total_transf.rot_mat[0]);
    Transform2D rigid_transf(total_transf.trans_vec[0], total_transf.trans_vec[1], angle);

    

    //Make new correlations based on the rigid total transform
    correlations.clear();
    for (int i = 0; i < las_vec.size(); i++){
        Correlation corr(map_kdt, las_vec[i], THOROUGH, rigid_transf);
        correlations.push_back(corr);
    }

    plot.add_data(las_vec, rigid_transf);
    plot.add_corrs(correlations, SINGLE);
    //If data points that were converged are now away from the corresponding map point, consider moving
    std::vector<double> sig_vec = {0,0};
    int sig_count = 0;
    double dist_score = 0;
    for (int i = 0; i < sig_ind.size(); i++){
        if (correlations[i].get_distance() > map_res){
            //Using an average right now, may need a max value in the future
            sig_vec[0] += correlations[i].get_trans()[0];
            sig_vec[1] += correlations[i].get_trans()[1];
            sig_count++;
        }
    }
    sig_vec[0] = sig_vec[0]/sig_count;
    sig_vec[1] = sig_vec[1]/sig_count;
    std::cout << "Sig vec is: [" << sig_vec[0] << ", " << sig_vec[1] << "]\n";

    //Use obtained sig vec as a maximum value for a bisection search
    double bisec = sqrt(sig_vec[0]*sig_vec[0] + sig_vec[1]*sig_vec[1]);
    std::cout << bisec << '\n';
    Transform2D base_t = rigid_transf;
    Transform2D far_t = rigid_transf;
    std::vector<double> dist_vec = sig_vec;
    far_t.trans_vec[0] -= dist_vec[0];
    far_t.trans_vec[1] -= dist_vec[1];

    //Evaluate both sides
    double dist_score_base = 0;
    double dist_score_far = 0;
    for (int i = 0; i < sig_ind.size(); i++){
        Correlation corr_base(map_kdt, las_vec[sig_ind[i]], THOROUGH, base_t);
        dist_score_base += corr_base.get_distance();
        Correlation corr_far(map_kdt, las_vec[sig_ind[i]], THOROUGH, far_t);
        dist_score_far += corr_far.get_distance();
    }

    Transform2D mid_t = rigid_transf;
    std::cout << "Entering bisection...\n";
    while(bisec > map_res/2){
        mid_t.trans_vec[0] = (base_t.trans_vec[0] + far_t.trans_vec[0])/2;
        mid_t.trans_vec[1] = (base_t.trans_vec[1] + far_t.trans_vec[1])/2;
        //Evaluate the middle
        double dist_score_mid = 0;
        for (int i = 0; i < sig_ind.size(); i++){
            Correlation corr_mid(map_kdt, las_vec[sig_ind[i]], THOROUGH, mid_t);
            dist_score_mid += corr_mid.get_distance();
        }
        if (dist_score_mid < dist_score_base && dist_score_mid > dist_score_far){
            //Here the answer lies between mid and far, re-assign mid to base and try again
            base_t = mid_t;
            bisec = bisec/2;
        }
        else if (dist_score_mid > dist_score_base && dist_score_mid < dist_score_far){
            //Here the answer lies between base and mid, re-assign far to mid and try again
            far_t = mid_t;
            bisec = bisec/2;
        }
        else if ((dist_score_mid < dist_score_base && dist_score_mid < dist_score_far) ||\
         (dist_score_mid > dist_score_base && dist_score_mid > dist_score_far)){
            //There are multiple minima in the range, need to pick the best one
            if (dist_score_base < dist_score_far){
                far_t = mid_t;
                bisec = bisec/2;
            }
            else{
                base_t = mid_t;
                bisec = bisec/2;
            }
            break;
        }
    }
    plot.add_data(las_vec, mid_t);

    end = clock();
    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    std::cout << "Time for execution: " << time_taken << "s\n";

    plot.plot_data();

    return 0;
}
