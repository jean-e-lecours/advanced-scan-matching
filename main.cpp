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
    if (!scan_data.read_from_file("../dat/test_1.txt")){
        scan_data.read_from_file("dat/test_1.txt");
    }
    //commentddddd
    TextData map_data;
    if (!map_data.read_from_file("../dat/test_map.txt")){
        map_data.read_from_file("dat/test_map.txt");
    }

    Transform2D las_transform(0.3,0.3,0);
    Transform2D map_transform(0,0,0);

    clock_t start, end;

    start = clock();

    std::vector<Point2D> las_vec = map_scan_points(las_transform, scan_data.dat_vec, 2*M_PI/360);
    std::vector<Point2D> pure_las_vec = las_vec;

    std::vector<Point2D> map = make_map(map_data.dat_vec, 100, 100, map_res);

    map_transform.transform(map);

    KdTree map_kdt(map);

    Plot2D plot(true);
    plot.add_data(map);
    plot.add_data(las_vec);

    Transform2D guess_transf(0,0,0);
    Transform2D total_transf(0,0,0);
    std::vector<Correlation> correlations;

    bool hdt_done = false;
    int guesses = 0;

    guess_transf.transform(las_vec);
    total_transf.add_transform(guess_transf);
    do{
        

        correlations.clear();
        double corr_mean = 0;
        double corr_stdev = 0;

        for (int i = 0; i < las_vec.size(); i++){
            Correlation corr(map_kdt, las_vec[i], THOROUGH, guess_transf);
            corr_mean += corr.corrected_value;
            correlations.push_back(corr);
        }
        corr_mean = corr_mean/correlations.size();

        for (int i = 0; i < correlations.size(); i++){
            corr_stdev += std::abs(correlations[i].corrected_value - corr_mean);
        }
        corr_stdev = sqrt(corr_stdev/correlations.size());
        corr_stdev = 1;

        std::vector<double> g = make_g_vector(correlations, las_vec, corr_stdev, 0.1, false);
        std::vector<double> G = make_G_matrix(correlations, las_vec, corr_stdev, 0.1, false);
        std::vector<double> x = solve_system(g, G, false);

        guess_transf.update_transform(x);
        guess_transf.print_transform();
        

        guess_transf.transform(las_vec);
        total_transf.add_transform(guess_transf);
        //plot.add_data(las_vec);
        
        guesses++;
        continue;
        
    } while(guess_transf.is_significant(threshold) && guesses < max_guesses);

    plot.add_corrs(correlations, SINGLE);

    correlations.clear();
    for (int i = 0; i < las_vec.size(); i++){
        Correlation corr(map_kdt, las_vec[i], THOROUGH, guess_transf);
        correlations.push_back(corr);
    }
    std::vector<double> sig_vec;
    for (int i = 0; i < correlations.size(); i++){
        sig_vec.push_back(correlations[i].get_distance());
    }

    double angle = std::atan2(total_transf.rot_mat[2],total_transf.rot_mat[0]);
    Transform2D rigid_transf(total_transf.trans_vec[0], total_transf.trans_vec[1], angle);
    total_transf = rigid_transf;
    total_transf.print_transform();
    total_transf.transform(pure_las_vec);

    correlations.clear();
    for (int i = 0; i < las_vec.size(); i++){
        Correlation corr(map_kdt, pure_las_vec[i], THOROUGH, guess_transf);
        correlations.push_back(corr);
    }
    plot.add_data(pure_las_vec);
    plot.add_corrs(correlations, SINGLE);

    std::vector<double> corr_trans = {0,0};
    std::vector<Point2D> sig_points;

    int sig_dist = 0;
    for (int i = 0; i < correlations.size(); i++){
        double new_dist = correlations[i].get_distance();
        if (sig_vec[i] < map_res && new_dist > map_res ){
            sig_points.push_back(correlations[i].scan);
            sig_dist++;
            corr_trans[0] += abs(correlations[i].scan.x - correlations[i].mcor1.x);
            corr_trans[1] += abs(correlations[i].scan.y - correlations[i].mcor1.y);
        }
    }
    corr_trans[0] = corr_trans[0]/sig_dist;
    corr_trans[1] = corr_trans[1]/sig_dist;
    Transform2D add_transf(corr_trans[0],corr_trans[1],0);
    add_transf.transform(pure_las_vec);
    plot.add_data(sig_points);
    plot.add_data(pure_las_vec);
    total_transf.add_transform(add_transf);
    total_transf.print_transform();
    //std::cout << "Affine distance score is: " << score << '\n';
    std::cout << sig_dist << " significant distances, consensus around " << corr_trans[0] << ", "<< corr_trans[1] << "\n";

    end = clock();
    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    std::cout << "Time for execution: " << time_taken << "s\n";

    plot.plot_data();

    return 0;
}
