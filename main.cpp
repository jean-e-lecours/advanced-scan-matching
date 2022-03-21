#include <bits/types/clock_t.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>

#include "hdr/mat.hpp"
#include "hdr/asm_extra.hpp"
#include "hdr/cpplot.hpp"

double* make_g_vector(Correlation* corrs, Set* data_set, double std_dev, double correntropy_factor, int data_size, bool print){
    double* g = new double[6];
    double delta = 0;
    double g_delta = 0;
    for (int i = 0; i < data_size; i++){
            delta = exp(- corrs[i].corrected_value * corrs[i].corrected_value / (2 * correntropy_factor * std_dev));
            g_delta = (corrs[i].model_corr_1->val[0] * corrs[i].norm_x + corrs[i].model_corr_1->val[1] * corrs[i].norm_y) * delta;
            //g = [nx⋅px⋅(mx⋅nx + my⋅ny)  nx⋅py⋅(mx⋅nx + my⋅ny)  ny⋅px⋅(mx⋅nx + my⋅ny)  ny⋅py⋅(mx⋅nx + my⋅ny)  nx⋅(mx⋅nx + my⋅ny)  ny⋅(mx⋅nx+ my⋅ny)]
            g[0] += corrs[i].norm_x * data_set->points[i].val[0] * g_delta;
            g[1] += corrs[i].norm_x * data_set->points[i].val[1] * g_delta;
            g[2] += corrs[i].norm_y * data_set->points[i].val[0] * g_delta;
            g[3] += corrs[i].norm_y * data_set->points[i].val[1] * g_delta;
            g[4] += corrs[i].norm_x * g_delta;
            g[5] += corrs[i].norm_y * g_delta;
    }
    if (print){
        //Print vector so that it can easily be imported into matlab if need be
        std::cout << "g=[" << g[0] << ',' << g[1] << ',' << g[2] << ',' << g[3] << ',' << g[4] << ',' << g[5] << "]\n\n";
    }
    return g;
}

double* make_G_matrix(Correlation* corrs, Set* data_set, double std_dev, double correntropy_factor, int data_size, bool print){
    double* G = new double[36];
    double delta = 0;
    for (int i = 0; i < data_size; i++){
        delta = exp(- corrs[i].corrected_value * corrs[i].corrected_value / (2 * correntropy_factor * std_dev));

    /*    ⎡     2   2        2                   2                       2                    ⎤
        ⎢ 0 nx ⋅px     1 nx ⋅px⋅py   2 nx⋅ny⋅px   3 nx⋅ny⋅px⋅py    4 nx ⋅px       5 nx⋅ny⋅px⎥
        ⎢                                                                                   ⎥
        ⎢     2              2   2                            2        2                    ⎥
        ⎢ 6 nx ⋅px⋅py    7 nx ⋅py    8 nx⋅ny⋅px⋅py  9 nx⋅ny⋅py     10 nx ⋅py     11 nx⋅ny⋅py⎥
        ⎢                                                                                   ⎥
        ⎢           2                        2   2         2                           2    ⎥
        ⎢12 nx⋅ny⋅px   13 nx⋅ny⋅px⋅py   14 ny ⋅px     15 ny ⋅px⋅py  16 nx⋅ny⋅px   17 ny ⋅px ⎥
        ⎢                                                                                   ⎥
        ⎢                           2        2               2   2                     2    ⎥
        ⎢18 nx⋅ny⋅px⋅py  19 nx⋅ny⋅py    20 ny ⋅px⋅py    21 ny ⋅py    22 nx⋅ny⋅py  23 ny ⋅py ⎥
        ⎢                                                                                   ⎥
        ⎢     2               2                                            2                ⎥
        ⎢24 nx ⋅px       25 nx ⋅py     26 nx⋅ny⋅px    27 nx⋅ny⋅py     28 nx       29 nx⋅ny  ⎥
        ⎢                                                                                   ⎥
        ⎢                                    2              2                           2   ⎥
        ⎣30 nx⋅ny⋅px    31 nx⋅ny⋅py     32 ny ⋅px      33 ny ⋅py     34 nx⋅ny     35 ny     ⎦ */

        G[0] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[0] * data_set->points[i].val[0] * delta;
        G[7] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[1] * data_set->points[i].val[1] * delta;
        G[14] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[0] * delta;
        G[21] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[1] * data_set->points[i].val[1] * delta;
        G[28] += corrs[i].norm_x * corrs[i].norm_x * delta;
        G[35] += corrs[i].norm_y * corrs[i].norm_y * delta;

        G[1] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[0] * data_set->points[i].val[1] * delta; G[6] = G[1];//xx xy
        G[2] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[0] * delta; G[12] = G[2];//xy xx
        G[3] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[1] * delta; G[18] = G[3];//xy xy
        G[4] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[0] * delta;                              G[24] = G[4];//xx x
        G[5] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[0] * delta;                              G[30] = G[5];//xy x

        //G[8]           -- This is equal to G[3] so we don't compute it here --                                      //xy xy
        G[9] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[1] * data_set->points[i].val[1] * delta; G[19] = G[9];//xy yy
        G[10] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[1] * delta;                             G[25] = G[10];//xx y
        G[11] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[1] * delta;                             G[31] = G[11];//xy y

        G[15] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[1] * delta; G[20] = G[15];//yy xy
        //G[16]          -- This is equal to G[5] so we don't compute it here --                         //xy x
        G[17] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[0] * delta;                             G[32] = G[17];//yy x
        
        //G[22]          -- This is equal to G[11] so we don't compute it here --                         //xy y
        G[23] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[1] * delta;                             G[33] = G[23];//yy y
        G[29] += corrs[i].norm_x * corrs[i].norm_y * delta;                                                          G[34] = G[29];//xy

    }
    //Set values that reappear
    G[8] = G[3]; G[13] = G[3]; G[18] = G[3];
    G[16] = G[5]; G[26] = G[5]; G[30] = G[5];
    G[22] = G[11]; G[27] = G[11]; G[31] = G[11];

    //Print matrix so it can be easily imported into matlab if need be;
    if (print){
        std::cout << std::setprecision(30);
        std::cout << "G=[";
        for (int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                std::cout << G[6*(i) + j];
                if (j+1 < 6){
                    std::cout << ',';
                }
            }
            if (i+1 < 6){
                std::cout << ';';
            }
        }
        std::cout << "]\n\n";
    }

    return G;
}

double* solve_system(double* g, double* G, bool print){
    double* x = new double[6];
    double mult = 0;
    //This will make the G matrix triangular using Gaussian Elimination
    //TODO: swap lines in case of bad pivot point
    for (int j = 0; j < 5; j++){
        for (int i = j+1; i < 6; i++){
            mult = G[6*i+j]/G[6*j+j];
            for (int k = 0; k < 6; k++){
                G[6*i+k] = G[6*i+k] - mult * G[6*j+k];
            }
            g[i] = g[i] - (mult * g[j]);
        }
    }
    //Solve the now triangular system
    double sum = 0;
    x[5] = g[5]/G[35];
    for (int i = 4; i >= 0; i--){
        sum = 0;
        for (int j = 5; j >= i+1; j--){
            sum += G[6*i+j] * x[j];
        }
        x[i] = (g[i] - sum)/G[6*i+i];
    }
    //Print the solution vector to double check if need be
    if (print){
        std::cout << "Solution vector is: ";
        for (int i = 0; i < 6; i++){
            std::cout << x[i] << '\t';
        }
        std::cout << "\n\n";
    }
    return x;
}

int Point::dims = 2;
int Set::dims = 2;
int KdTree::dims = 2;

int main(){
    //Parameters that would be set on ROS
    int max_guesses = 100;
    double correntropy_factor = 0.1;
    double transform_tresh = 0.00001;

    clock_t start, end;
    
    //Parameters we just need initialized
    int guesses = 0;
    double transform_diff = INFINITY;
    //This will serve as an initializer and a way to keep track of the total transforms
    Transform2D total_transform(0, 0, 0);
    Transform2D base_transform(0, 0, 0);
    Transform2D true_transform(0, 0, 0.6);

    start = clock();

    //Info we would get from the laser scan message
    double scan_period = (2*M_PI)/360;
    TextLaserScanData laser_scan_1(false);
    laser_scan_1.read_from_file("dat/test_dat/test_3_1.txt");
    //TextLaserScanData laser_scan_2(false);
    //laser_scan_2.read_from_file("incl/test_dat/test_3_2.txt");
    
    
    //Will need something to interpret the laser scan message
    int data_size = laser_scan_1.usable_las_size;
    //laser_scan_1.filter_scan_points();

    Set model_set = laser_scan_1.map_scan_points(&base_transform, scan_period);
    Set data_set = laser_scan_1.map_scan_points(&base_transform, scan_period);

    true_transform.transform_set(&data_set);

    Plot2D plot(true, max_guesses);

    KdTree model_tree(&model_set);
    model_set.print_set();

    Transform2D guess_transform(0, 0, 0);

    //Get first and second correlations in data set to model set

    plot.add_data(&model_set);
    plot.add_data(&data_set);


    while (transform_diff > transform_tresh && guesses < max_guesses){

        //Apply guess transform on the set and add it to the total transform
        //The guess transform is with respect to the previous guess
        guess_transform.transform_set(&data_set);
        total_transform.add_transform(&guess_transform);

        Correlation* corrs = new Correlation[data_size];
        //Get correlation between points and average for the standard deviation
        double corr_mean = 0;
        for (int i = 0; i < data_size; i++){
            corrs[i] = Correlation(&model_tree, &data_set, i, &guess_transform);
            corr_mean += corrs[i].corrected_value;
            if (guesses > 13){
                plot.add_line(corrs[i].model_corr_1, &corrs[i].current_point);
            }
            //plot.add_line(corrs[i].model_corr_1, &corrs[i].current_point);
        }
        std::cout << guesses << '\n';
        corr_mean = corr_mean/data_size;
        plot.add_data(&data_set);

        //Calculate the standard deviation of the corrected value estimate
        double std_dev = 0;
        for (int i = 0; i < data_size; i++){
            std_dev += abs(corrs[i].corrected_value - corr_mean);
        }
        std_dev = sqrt(std_dev/data_size);

        //For two dimensions, g is a 1x6 vector and G is a 6x6 matrix
        double* g = make_g_vector(corrs, &data_set, std_dev, correntropy_factor, data_size, false);
        double* G = make_G_matrix(corrs, &data_set, std_dev, correntropy_factor, data_size, false);
        //Solve the system using gaussian elimination
        //TODO: skip or/and error message on receiving unsolvable system
        double* x = solve_system(g, G, false);

        //Evaluate difference and update guess Transform2D and update
        transform_diff = guess_transform.compare_transform(x);
        guess_transform.update_transform(x);

        std::cout << "Guessed translation of: " << guess_transform.trans_vec[0] << ", " << guess_transform.trans_vec[1] << '\n';
        guesses++;
    }
    guess_transform.transform_set(&data_set);
    total_transform.add_transform(&guess_transform);
    plot.add_data(&data_set);

    double angle = std::atan2(total_transform.rot_mat[2],total_transform.rot_mat[0]);
    

    Transform2D simple_transform(total_transform.trans_vec[0],total_transform.trans_vec[1],angle);
    simple_transform.print_transform();
    std::cout << angle << '\n';

    /* double* AtA = affine_to_rot(total_transform.rot_mat, dimensions, true);

    Transform2D test_transform(total_transform.trans_vec[0],total_transform.trans_vec[1],0);
    test_transform.rot_mat[0] = AtA[0];
    test_transform.rot_mat[1] = AtA[1];
    test_transform.rot_mat[2] = AtA[2];
    test_transform.rot_mat[3] = AtA[3];
    double angle2 = std::atan2(test_transform.rot_mat[2],test_transform.rot_mat[0]);
    //This is a transform with a proper rotation matrix
    test_transform.print_transform();
    std::cout << angle2 << '\n'; */
    //TODO: adjust translation vector for the new proper rotation matrix


    std::cout << "Final transform is:\n";
    total_transform.print_transform();
    end = clock();
    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    std::cout << "Time for execution: " << time_taken << "s\n";
    plot.plot_data();
    return 0;
}