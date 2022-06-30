#include <cmath>
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
    if (!scan_data.read_from_file("../dat/test_6.txt")){
        scan_data.read_from_file("dat/test_6.txt");
    }
    //commentddddd
    TextData map_data;
    if (!map_data.read_from_file("../dat/test_map.txt")){
        map_data.read_from_file("dat/test_map.txt");
    }

    Transform2D las_transform(2,0.5,1.5);
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
        plot.add_corrs(correlations, SINGLE);
        corr_mean = corr_mean/correlations.size();

        for (int i = 0; i < correlations.size(); i++){
            corr_stdev += std::abs(correlations[i].corrected_value - corr_mean);
        }
        corr_stdev = sqrt(corr_stdev/correlations.size());

        std::vector<double> g = make_g_vector(correlations, corr_stdev, 0.1, false);
        std::vector<double> G = make_G_matrix(correlations, corr_stdev, 0.1, false);
        std::vector<double> x = solve_system(g, G, false);

        guess_transf.update_transform(x);

        double angle = std::atan2(guess_transf.rot_mat[2],guess_transf.rot_mat[0]);
        Transform2D rigid_guess(guess_transf.trans_vec[0], guess_transf.trans_vec[1], angle);

        guess_transf = rigid_guess;
        guess_transf.print_transform();
        total_transf.add_transform(guess_transf);
        
        guesses++; 
        
    } while(guess_transf.is_significant(threshold) && guesses < max_guesses);
    total_transf.transform(las_vec);

    plot.add_data(las_vec);

    std::vector<std::vector<double>> sig_vec;
    for (int i = 0; i < correlations.size(); i++){
        sig_vec.push_back(correlations[i].get_trans());
    }
    Plot2D sig_plot(true, "densig_");
    sig_plot.add_vec_data_2d(sig_vec);
    sig_plot.plot_data();

    std::vector<double> sig_trans = {0,0};
    double sign_x;
    double sign_y;
    for (int i = 0; i < sig_vec.size(); i++){
        if (std::abs(sig_vec[i][0]) > map_res){
            sig_trans[0] += sig_vec[i][0];
            sign_x ++;
        }
        if (std::abs(sig_vec[i][1]) > map_res){
            sig_trans[1] += sig_vec[i][1];
            sign_y ++;
        }
    }
    sig_trans[0] = sig_trans[0]/sign_x;
    sig_trans[1] = sig_trans[1]/sign_y;

    Transform2D add_transf(-sig_trans[0], -sig_trans[1], 0);
    total_transf.add_transform(add_transf);
    add_transf.transform(las_vec);
    plot.add_data(las_vec);
    end = clock();
    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    std::cout << "Time for execution: " << time_taken << "s\n";

    plot.plot_data();

    return 0;
}
