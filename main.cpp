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

    Transform2D las_transform(2,0.7,1.5);
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

    Transform2D guess_transf(0,0,0);
    Transform2D total_transf(0,0,0);
    std::vector<Correlation> correlations;

    bool hdt_done = false;
    int guesses = 0;
    do{
        guess_transf.transform(las_vec);
        total_transf.add_transform(guess_transf);

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
        double angle = std::atan2(guess_transf.rot_mat[2],guess_transf.rot_mat[0]);
        Transform2D rigid_transf(guess_transf.trans_vec[0], guess_transf.trans_vec[1], angle);
        guess_transf = rigid_transf;
        guess_transf.print_transform();

        plot.add_data(las_vec);
        if (guesses == 4){
            plot.add_corrs(correlations, SINGLE);
        }
        
        guesses++;
        continue;
        
    } while(guess_transf.is_significant(threshold) && guesses < max_guesses);
    

    double score = 0;
    int sig_dist = 0;
    std::vector<double> sig_vec = {0,0};
    for (int i = 0; i < correlations.size(); i++){
        double dist = correlations[i].get_distance();
        score += dist;
        if (dist > 1.5*map_res){
            sig_vec[0] += correlations[i].scan.x - correlations[i].mcor1.x;
            sig_vec[1] += correlations[i].scan.y - correlations[i].mcor1.y;
            sig_dist++;
        }
    }
    sig_vec[0] = sig_vec[0]/correlations.size();
    sig_vec[1] = sig_vec[1]/correlations.size();
    total_transf.print_transform();
    std::cout << "Affine distance score is: " << score << '\n';
    std::cout << sig_dist << " significant distances, consensus around " << sig_vec[0] << ", "<< sig_vec[1] << "\n";

    /* if (abs(sig_vec[0]) > map_res || sig_vec[1] > map_res){
        std::cout << "Correcting translation...\n";
        Transform2D sig_correct(-sig_vec[0], -sig_vec[1],0);
        sig_correct.transform(las_vec);
        plot.add_data(las_vec);
    }

    correlations.clear();
        double corr_mean = 0;
        double corr_stdev = 0;

        for (int i = 0; i < las_vec.size(); i++){
            Correlation corr(map_kdt, las_vec[i], THOROUGH, guess_transf);
            corr_mean += corr.corrected_value;
            correlations.push_back(corr);
        }

    score = 0;
    sig_dist = 0;
    sig_vec = {0,0};
    for (int i = 0; i < correlations.size(); i++){
        double dist = correlations[i].get_distance();
        score += dist;
        if (dist > 1.5*map_res){
            sig_vec[0] += correlations[i].scan.x - correlations[i].mcor1.x;
            sig_vec[1] += correlations[i].scan.y - correlations[i].mcor1.y;
            sig_dist++;
        }
    }
    sig_vec[0] = sig_vec[0]/correlations.size();
    sig_vec[1] = sig_vec[1]/correlations.size();
    total_transf.print_transform();
    std::cout << "Affine distance score is: " << score << '\n';
    std::cout << sig_dist << " significant distances, consensus around " << sig_vec[0] << ", "<< sig_vec[1] << "\n"; */



    if (trans){
        correlations.clear();
        double angle = std::atan2(total_transf.rot_mat[2],total_transf.rot_mat[0]);
        Transform2D rigid_transf(total_transf.trans_vec[0], total_transf.trans_vec[1], angle);
        Transform2D pure_trans(0,0,0);
        
        Transform2D trans_correction(0,0,0);

        rigid_transf.transform(pure_las_vec);
        do{
            std::vector<double> mean_trans = {0,0};
            correlations.clear();
            for (int i = 0; i < pure_las_vec.size(); i++){
                Correlation corr(map_kdt, pure_las_vec[i], THOROUGH, pure_trans);
                std::vector<double> temp_trans = corr.get_trans();
                mean_trans[0] -= temp_trans[0];
                mean_trans[1] -= temp_trans[1];
                correlations.push_back(corr);
            }

            trans_correction.trans_vec[0] = mean_trans[0]/correlations.size();
            trans_correction.trans_vec[1] = mean_trans[1]/correlations.size();

            trans_correction.transform(pure_las_vec);
            rigid_transf.add_transform(trans_correction);
            
        }while(trans_correction.is_significant(0.001));

        plot.add_data(pure_las_vec);
         
        score = 0;
        for (int i = 0; i < correlations.size(); i++){
            score += correlations[i].get_distance();
        }
        rigid_transf.print_transform();
        std::cout << "Rigid distance score is: " << score << '\n';
    }

    end = clock();
    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    std::cout << "Time for execution: " << time_taken << "s\n";

    plot.plot_data();

    return 0;
}
