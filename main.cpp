#include <cmath>
#include <iostream>
#include <ctime>
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

    char algorithm = ASM;
    int max_guesses = 5;
    double threshold = 0.01;

    TextData scan_data;
    if (!scan_data.read_from_file("../dat/test_5.txt")){
        scan_data.read_from_file("dat/test_5.txt");
    }

    TextData map_data;
    if (!map_data.read_from_file("../dat/test_map.txt")){
        map_data.read_from_file("dat/test_map.txt");
    }

    Transform2D las_transform(2,0,0.3);
    Transform2D map_transform(0,0,0);

    clock_t start, end;

    start = clock();

    std::vector<Point2D> las_vec = map_scan_points(las_transform, scan_data.dat_vec, 2*M_PI/360);

    std::vector<Point2D> map = make_map(map_data.dat_vec, 100, 100, 0.025);

    map_transform.transform(map);

    KdTree map_kdt(map);

    Plot2D plot(true);
    plot.add_data(map);

    Transform2D guess_transf(0,0,0);
    Transform2D total_transf(0,0,0);

    bool hdt_done = false;
    int guesses = 0;
    do{
        guess_transf.transform(las_vec);
        total_transf.add_transform(guess_transf);

        std::vector<Correlation> correlations;
        double corr_mean = 0;
        double corr_stdev = 0;

        for (int i = 0; i < las_vec.size(); i++){
            Correlation corr(map_kdt, las_vec[i], THOROUGH, las_transform);
            corr_mean += corr.corrected_value;
            correlations.push_back(corr);
        }
        switch (algorithm) {
            case LSSM:{
                std::vector<double> g = lssm::make_g_vector(correlations, las_vec, false);
                std::vector<double> M = lssm::make_M_matrix(correlations, las_vec, false);
                std::vector<double> x = lssm::solve_system(g, M, false);
                
                lssm::update_transform(guess_transf, x);
                
                if (!guess_transf.is_significant(0.001)){
                    guess_transf.print_transform();

                    plot.add_data(las_vec);
                    plot.add_corrs(correlations, SINGLE);
                    break;
                }

                guess_transf.print_transform();

                plot.add_data(las_vec);
                //plot.add_corrs(correlations, SINGLE);
                continue;
            }
            case HDSM:{

                if (!hdt_done){
                    double tgrad_step = 0.1;
                    std::vector<Transform2D> tgrad_transfs {Transform2D(-tgrad_step, tgrad_step, 0),  Transform2D(0, tgrad_step, 0),  Transform2D(tgrad_step, tgrad_step, 0),\
                                                        Transform2D(-tgrad_step, 0, 0),          Transform2D(0, 0, 0),          Transform2D(tgrad_step, 0, 0),\
                                                        Transform2D(-tgrad_step, -tgrad_step, 0), Transform2D(0, -tgrad_step, 0), Transform2D(tgrad_step, -tgrad_step, 0)};
                    //Translation only step conducted first
                    int transf_ind;
                    double best_hd = INFINITY;
                    double new_hd;
                    for (int t = 0; t < tgrad_transfs.size(); t++){
                        new_hd = find_hd(correlations, tgrad_transfs[t], 5);
                        if (new_hd < best_hd){
                            best_hd = new_hd;
                            transf_ind = t;
                        }
                    }
                    if (transf_ind == 4){
                        std::cout << "As good as it gets for translation!\n";
                        hdt_done = true;
                        break;
                    }
                    guess_transf = tgrad_transfs[transf_ind];
                    
                    continue;
                }
                else{
                    double agrad_step = 0.01;

                    std::vector<Transform2D> agrad_transfs {Transform2D(0, 0, -agrad_step), Transform2D(0, 0, 0), Transform2D(0, 0, agrad_step)};

                    //Rotation only step conducted after
                    int transf_ind;
                    double best_hd = INFINITY;
                    double new_hd;
                    for (int t = 0; t < agrad_transfs.size(); t++){
                        new_hd = find_hd(correlations, agrad_transfs[t], 50);
                        if (new_hd < best_hd){
                            best_hd = new_hd;
                            transf_ind = t;
                        }
                    }
                    if (transf_ind == 1){
                        std::cout << "As good as it gets for rotation!\n";
                        hdt_done = true;
                        break;
                    }
                    guess_transf = agrad_transfs[transf_ind];
                    
                    continue;
                }
            }
            case ASM:{
                corr_mean = corr_mean/correlations.size();

                for (int i = 0; i < correlations.size(); i++){
                    corr_stdev += std::abs(correlations[i].corrected_value - corr_mean);
                }
                corr_stdev = sqrt(corr_stdev/correlations.size());
                std::cout << corr_stdev << '\n';

                std::vector<double> g = make_g_vector(correlations, las_vec, corr_stdev, 0.1, false);
                std::vector<double> G = make_G_matrix(correlations, las_vec, corr_stdev, 0.1, false);
                std::vector<double> x = solve_system(g, G, false);

                guess_transf.update_transform(x);
                guess_transf.print_transform();

                plot.add_data(las_vec);
                plot.add_corrs(correlations, SINGLE);
                guesses++;
                continue;
            }
        }
        
    } while(guess_transf.is_significant(threshold) && guesses < max_guesses);
    
    plot.plot_data();

    end = clock();

    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    std::cout << "Time for execution: " << time_taken << "s\n";

    return 0;
}
