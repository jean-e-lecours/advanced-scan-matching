#include <cmath>
#include <ctime>
#include <string>

#include "hdr/mat.hpp"
#include "hdr/asm_extra.hpp"
#include "hdr/cpplot.hpp"

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
            std_dev += std::abs(corrs[i].corrected_value - corr_mean);
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