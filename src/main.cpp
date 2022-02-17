#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gnuplot.h>
#include <math.h>
#include <string>
#include <sys/types.h>
#include <vector>

struct NODE{
    double pos_val[2];
    bool dimension;
    NODE* small_link;
    NODE* big_link;
    NODE* parent;
    int layer;
};

struct POINT{
    double x;
    double y;
};

class Set{

    public:
        int data_size;
        int max_layer = 0;
        double* x_vals;
        double* y_vals;
        double** axes_ptr[2] = {&x_vals, &y_vals};

        Set* copy_set(){
            Set* new_set = new Set(data_size, x_vals, y_vals);
            return new_set;
        }

        void print_data(){
            for (int i = 0; i < data_size; i++){
                std::cout << *(*axes_ptr[0]+i) << ", " << *(*axes_ptr[1]+i) << '\n';
            }
        }

    Set(int data_size, double* x_vals, double* y_vals){
        //Constructor
        this->data_size = data_size;
        this->x_vals = new double[data_size];
        this->y_vals = new double[data_size];

        for (int i = 0; i < data_size; i++){
            this->x_vals[i] = x_vals[i];
            this->y_vals[i] = y_vals[i];
        }
    }
    ~Set(){
        //Deconstructor
    }
};

class KdTree{
    NODE* kd_node_array;
    int node_index = 0;

    void populate_node(bool dimension, int start, int points_left, NODE* this_node, NODE* starting_node){
        int index = start + std::ceil(points_left/2.0)- 1;
        this_node->dimension = dimension;
        this_node->pos_val[0] = x_vals[index];
        this_node->pos_val[1] = y_vals[index];
        this_node->parent = starting_node;
        bool size = false;
        if (starting_node != NULL){
            this_node->layer = starting_node->layer + 1;
            //Node linking, not needed for the first node (has a NULL starting_node)
            if (this_node->pos_val[!dimension] < starting_node->pos_val[!dimension]){
                starting_node->small_link = this_node;
            }
            else{
                starting_node->big_link = this_node;
            }
        }
        else{
            this_node->layer = 0;
        }
    }

    public:
        int data_size;
        int max_layer = 0;
        double* x_vals;
        double* y_vals;
        double** axes_ptr[2] = {&x_vals, &y_vals};
        bool first_dimension = false;

        void print_kd_tree(){
            for (int i = 0; i < data_size; i++) {
                std::cout << kd_node_array[i].pos_val[0] << ", " << kd_node_array[i].pos_val[1];
                std::cout << "[Layer " << kd_node_array[i].layer << "] ";
                if (kd_node_array[i].big_link != NULL){
                    std::cout << " =+=> " << kd_node_array[i].big_link->pos_val[0] << ", " << kd_node_array[i].big_link->pos_val[1];
                }
                if (kd_node_array[i].small_link != NULL){
                    std::cout << " =-=> " << kd_node_array[i].small_link->pos_val[0] << ", " << kd_node_array[i].small_link->pos_val[1];
                }
                std::cout << '\n';
            }
            std::cout << max_layer << '\n';
        }

        NODE* find_approximate_closest_point(double* target_point, bool dimension){
            NODE* node_ptr = NULL;
            NODE* next_node_ptr = &kd_node_array[0];

            while (next_node_ptr != NULL){
                node_ptr = next_node_ptr;
                if (target_point[0] < node_ptr->pos_val[dimension]){
                    next_node_ptr = node_ptr->small_link;
                }
                else{
                    next_node_ptr = node_ptr->big_link;
                }
                dimension = !dimension;
            }
            return node_ptr;
        }

        NODE** find_closest_point(double* target_point, bool dimension, bool find_second_closest){
            NODE** possible_node_ptr = new NODE*[(max_layer * (max_layer + 1))/2];
            NODE** best_node_ptr = new NODE*[2];
            best_node_ptr[0] = NULL;
            best_node_ptr[1] = NULL;
            int possible_node_index = 0;
            NODE* node_ptr = &kd_node_array[0];
            double distance = 0;
            double border_distance = 0;
            double closest_distance = INFINITY;
            double second_closest_distance = INFINITY;

            while (true){
                while (node_ptr != NULL){
                    //We only go forwards and end when the pointer dies, unless there is
                    //another branch to inspect.
                    distance = sqrt((target_point[0] - node_ptr->pos_val[0]) * (target_point[0] - node_ptr->pos_val[0]) + (target_point[1] - node_ptr->pos_val[1]) * (target_point[1] - node_ptr->pos_val[1]));
                    if (distance < closest_distance){
                        if (find_second_closest){
                            best_node_ptr[1] = best_node_ptr[0];
                            second_closest_distance = closest_distance;
                        }
                        best_node_ptr[0] = node_ptr;
                        closest_distance = distance;
                    }
                    else if (find_second_closest && distance < second_closest_distance){
                        best_node_ptr[1] = node_ptr;
                        second_closest_distance = distance;
                    }

                    dimension = node_ptr->dimension;
                    border_distance = target_point[dimension] - node_ptr->pos_val[dimension];
                    //If the distance to the node is larger than the distance to the border,
                    //there might be a closer point on the other side.
                    if (border_distance < 0){
                        //The point is to the left of the node, need to refer to small link
                        if (distance > abs(border_distance)){
                            //Point is close to the border, will keep the big link in mind for later as well
                            possible_node_ptr[possible_node_index] = node_ptr->big_link;
                            possible_node_index++;
                        }
                        node_ptr = node_ptr->small_link;
                    }
                    else{
                        //The point is equal to or to the right of the node, need to refer to big link
                        if (distance > abs(border_distance)){
                            //Point is close to the border, will keep the small link in mind for later as well
                            possible_node_ptr[possible_node_index] = node_ptr->small_link;
                            possible_node_index++;
                        }
                        node_ptr = node_ptr->big_link;
                    }
                }

                if (possible_node_index != 0){
                    //The possible node pointer will always be 1 ahead of where it should be.
                    //If it is 0, then we would want to look at -1 and that is no good.
                    node_ptr = possible_node_ptr[possible_node_index - 1];
                    possible_node_index--;
                }
                else{
                    //This is the exit condition for the loop, better than a do while in my
                    //opinion (saves a branch).
                    break;
                }
            }
            if (find_second_closest){
                return best_node_ptr;
            }
            else{
                return &best_node_ptr[0];
            }
            
        }

        void make_kd_tree(bool dimension, int start, int points_left, NODE* starting_node){
            for (int i = start; i < (points_left + start); i++){
                for (int j = i; j < (points_left + start); j++){
                    if (*(*axes_ptr[dimension] + i) > *(*axes_ptr[dimension] + j)){
                        double temp_x = x_vals[i];
                        double temp_y = y_vals[i];
                        x_vals[i] = x_vals[j];
                        y_vals[i] = y_vals[j];
                        x_vals[j] = temp_x;
                        y_vals[j] = temp_y;
                    }
                }
            }
            populate_node(dimension, start, points_left, &kd_node_array[node_index], starting_node);
            NODE* node = &kd_node_array[node_index];
            if (node->layer > max_layer){
                max_layer = node->layer;
            }
            node_index++;

            if (points_left <= 1){
                //If we hit a dead end we return and leave the links NULL
                //std::cout << kd_node_array[node_index].x_val << '\n';
                return;
            }
            else if (points_left <= 2) {
                make_kd_tree(!dimension, start + 1, 1, node);
            }
            else{
                make_kd_tree(!dimension, start, std::ceil(points_left/2.0 - 1), node);
                make_kd_tree(!dimension, start + std::ceil(points_left/2.0), std::floor(points_left/2.0), node);
                return;
            }
        }

    KdTree(Set parent_set, bool first_dimension){
        //Constructor
        this->x_vals = parent_set.x_vals;
        this->y_vals = parent_set.y_vals;
        this->data_size = parent_set.data_size;
        this->first_dimension = first_dimension;

        this->kd_node_array = new NODE[data_size];

        //Calls a recursive function to make create the object and all its properties
        make_kd_tree(false, 0, data_size, NULL);
    }
    ~KdTree(){
        //Destructor
    }
};

class Transform{
    public:
        double rot_mat[4];
        double trans_vec[2];
        double z_rot;

        void transform_set(Set* set){
            for (int i = 0; i < set->data_size; i++){
                //Applies a single matrix multiplication with the provided x, y vector from the set
                double temp_x = rot_mat[0] * set->x_vals[i] + rot_mat[1] * set->y_vals[i] + trans_vec[0];
                double temp_y = rot_mat[2] * set->x_vals[i] + rot_mat[3] * set->y_vals[i] + trans_vec[1];
                set->x_vals[i] = temp_x;
                set->y_vals[i] = temp_y;
            }
        }
        NODE transform_referenced_point(double x, double y){
            double temp_x = rot_mat[0] * x + rot_mat[1] * y + trans_vec[0];
            double temp_y = rot_mat[2] * x + rot_mat[3] * y + trans_vec[1];
            NODE point;
            point.pos_val[0] = temp_x;
            point.pos_val[1] = temp_y;
            return point;
        }

    Transform(double x_trans, double y_trans, double z_rot){
        //Constructor
        //Sets the two elements of the transforms as per the desired movement
        rot_mat[0] = cos(z_rot);
        rot_mat[1] = -sin(z_rot);
        rot_mat[2] = -rot_mat[1];
        rot_mat[3] = rot_mat[0];
        trans_vec[0] = x_trans;
        trans_vec[1] = y_trans;

        this->z_rot = z_rot;
    }
    ~Transform(){
        //Deconstructor
    }
};

class Plot2D{
    bool data_permanence = false;
    int number_of_plots = 0;
    double padding = 1;
    double extremas[4] = {INFINITY, 0, INFINITY, 0};

    std::string gnup_line = "plot ";

    public:
        std::string* filenames;
        int total_plots;

        void add_data(Set* set){
            //Add file name to the array and create the file itself
            std::string filename = "dat/temp2d_" + std::to_string(number_of_plots)+ ".dat";
            std::ofstream osetf{filename};
            //Go through all the data to find what is the maximum and minimum in x
            //and y so that we can size the plot window properly. This is not
            //necessary but we are going through all the data anyways so why not.
            for (int i = 0; i < set->data_size; i++){ 
                if (set->x_vals[i] < extremas[0]){
                    extremas[0] = set->x_vals[i]; 
                }
                else if (set->x_vals[i] > extremas[1]){
                    extremas[1] = set->x_vals[i];
                }
                if (set->y_vals[i] < extremas[2]){
                    extremas[2] = set->y_vals[i]; 
                }
                else if (set->y_vals[i] > extremas[3]){
                    extremas[3] = set->y_vals[i];
                }
                //This is where we actually place the data into the file
                osetf << set->x_vals[i] << ' ' << set->y_vals[i] << '\n';
            }
            number_of_plots++;
            gnup_line += "\"" + filename + "\",";
        }

        void add_line(NODE* node1, NODE* node2){
            double x1 = node1->pos_val[0];
            double x2 = node2->pos_val[0];
            double y1 = node1->pos_val[1];
            double y2 = node2->pos_val[1];

            double a = (y2 - y1)/(x2 - x1);
            double b = y1 - a*x1;

            if (a != INFINITY && b != INFINITY){
                gnup_line += "[" + std::to_string(x1) + ":" + std::to_string(x2) + "] " + std::to_string(a) + "*x+" + std::to_string(b) + " notitle,";
            }
        }

        void plot_data(){
            //The plot is drawn when the pipe is destroyed. It is created here so that it
            //goes out of scope right after the function is done and before the files are
            //deleted when the Plot2D object is destroyed.
            GnuplotPipe gp;
            
            //Set graph bounds from the data gathered previously
            std::string gnup_range = "set xrange[" + std::to_string(extremas[0] - padding) + ":" + std::to_string(extremas[1] + padding) + "]\n" + \
                                     "set yrange[" + std::to_string(extremas[2] - padding) + ":" + std::to_string(extremas[3] + padding) + "]\n";
            gp.sendLine(gnup_range);
            gnup_line.erase(gnup_line.end() - 1);

            std::ofstream gnup_file{"gnup_line"};
            gnup_file << gnup_line << '\n';
            gnup_file.close();

            gp.sendLine(gnup_line);
        }

    Plot2D(bool keep_data){
        //Constructor
        this->data_permanence = keep_data;
    }
    ~Plot2D(){
        if (!data_permanence){
            for (int i = 0; i < total_plots; i++){
                remove(filenames[i].c_str());
            }
        }
    }
};

class LaserScanData{
    std::vector<double> las_vec;

    public:
        //Total number of points to compute xy position
        int las_size;
        //Amount of valid data
        int usable_las_size = 0;

        int read_from_file(std::string filename){
            std::ifstream ilasf{filename};
            if (!ilasf){
                std::cerr << "Could not read laser scan data!\n";
                return 0;
            }
            else{
                std::cout << "Reading data\n";
            }
            while (ilasf){
                std::string las_read;
                ilasf >> las_read;
                int las_length = las_read.length()-1;
                if (las_read[las_length] == ','){
                    las_read.erase(las_read.end()-1);
                }
                
                double las_val = 0;
                if (las_length > 2){
                    las_val = std::stod(las_read);
                    las_vec.push_back(las_val);
                }
                if (las_val < INFINITY){
                    usable_las_size++;
                }
                
            }
            ilasf.close();
            las_size = las_vec.size();
            return 1;
        }
        Set map_scan_points(Transform* transform, double scan_period){
            double* las_x = new double[usable_las_size];
            double* las_y = new double[usable_las_size];
            int j = 0;
            for (int i = 0; i < las_size; i++){
                if (las_vec[i] < 5){
                    double temp_scan_x = las_vec[i] * cos(scan_period * i);
                    double temp_scan_y = las_vec[i] * sin(scan_period * i);
                    las_x[j] = transform->trans_vec[0] + temp_scan_x * cos(transform->z_rot) - temp_scan_y * sin(transform->z_rot);
                    las_y[j] = transform->trans_vec[1] + temp_scan_y * cos(transform->z_rot) - temp_scan_x * sin(transform->z_rot);
                    j++;
                }
            }
            Set las_set(usable_las_size, las_x, las_y);
            return las_set;
        }

    LaserScanData(){
    //Constructor
    }
    ~LaserScanData(){
    //Destructor
    }
};

int main(){
    double guesses = 0;
    double correntropy_factor = 0.1;
    double transform_tresh = 0.001 ;
    double transform_diff = INFINITY;

    Transform base_transform(0, 0, 0);
    Transform true_transform(-0.1, -0.1, 0.1);

    /* double x_vals[9] = {-1, -0.5, 0, 0.5, 1, 1, 1, 0.5, 1};
    double y_vals[9] = {1, 0.5, 1, 1, 1, 0.5, 0, 0, -1};
    int data_size = 9;

    Set model_set(data_size, x_vals, y_vals);
    Set data_set(data_size, x_vals, y_vals);
    true_transform.transform_set(&data_set); */

    double scan_period = (2*M_PI)/360;
    LaserScanData laser_scan;
    laser_scan.read_from_file("incl/test_1_range_data.txt");
    int data_size = laser_scan.usable_las_size;

    Set model_set = laser_scan.map_scan_points(&base_transform, scan_period);
    Set data_set = laser_scan.map_scan_points(&true_transform, scan_period);


    Plot2D plot(false);

    KdTree model_tree(model_set, false);

    Transform guess_transform(0, 0, 0);

    //Get first and second correlations in data set to model set
    NODE** model_corr_1 = new NODE*[data_size];
    NODE** model_corr_2 = new NODE*[data_size];
    NODE** model_temp;

    plot.add_data(&model_set);
    

    while (transform_diff > transform_tresh && guesses < 25){

        //  --  This section works fine  --
        //Set transformed_set = *data_set.copy_set();
        Set transformed_set = data_set;
        guess_transform.transform_set(&transformed_set);
        for (int i = 0; i < data_size; i++){
            //Making copies of individual point
            double current_point_arr[2] = {transformed_set.x_vals[i], transformed_set.y_vals[i]};
            NODE current_point;
            current_point.pos_val[0] = transformed_set.x_vals[i];
            current_point.pos_val[1] = transformed_set.y_vals[i];
            //Doing some weird pointer thing. Looks like this is working.
            model_temp = model_tree.find_closest_point(current_point_arr, false, true);
            model_corr_1[i] = model_temp[0]; //mi1
            model_corr_2[i] = model_temp[1]; //mi2

            plot.add_line(model_corr_1[i], &current_point);
            //plot.add_line(data_corr_2[i], &current_point);
        }
        plot.add_data(&transformed_set);

        //Calculate norms from the correlations for p2l and directly calculate all the corrected point values n*[Ap+t-m]
        double* norm_x = new double[data_size];
        double* norm_y = new double[data_size];
        double norm_size = 0;
        double* corrected_value = new double[data_size];
        for (int i = 0; i < data_size; i++){
            //Sign of the norm does not matter, using cross product with [0;0;1]
            norm_x[i] = model_corr_1[i]->pos_val[1] - model_corr_2[i]->pos_val[1];
            norm_y[i] = - (model_corr_1[i]->pos_val[0] - model_corr_2[i]->pos_val[0]);
            norm_size = sqrt(norm_x[i]*norm_x[i] + norm_y[i]*norm_y[i]);
            norm_x[i] = norm_x[i]/norm_size;
            norm_y[i] = norm_y[i]/norm_size;
            
            //CHECK: I assume this is correct for now, to be investigated further if stuff doesn't work as expected
            corrected_value[i] = norm_x[i]*(guess_transform.rot_mat[0]*transformed_set.x_vals[i] + guess_transform.rot_mat[1]*transformed_set.y_vals[i] + guess_transform.trans_vec[0] - model_corr_1[i]->pos_val[0]) + \
                                norm_y[i]*(guess_transform.rot_mat[2]*transformed_set.x_vals[i] + guess_transform.rot_mat[3]*transformed_set.y_vals[i] + guess_transform.trans_vec[1] - model_corr_1[i]->pos_val[1]);
        }

        //Calculate the standard deviation of the corrected value estimate
        double mean = 0;
        for (int i = 0; i < data_size; i++){
            mean += corrected_value[i];
        }
        mean = mean/data_size;
        double std_dev = 0;
        for (int i = 0; i < data_size; i++){
            std_dev += abs(corrected_value[i] - mean);
        }
        std_dev = sqrt(std_dev/data_size);

        //Compute matrices g and G for the final equation
        double g[6] = {0, 0, 0, 0, 0, 0};
        double G[36] = {0, 0, 0, 0, 0, 0, \
                        0, 0, 0, 0, 0, 0, \
                        0, 0, 0, 0, 0, 0, \
                        0, 0, 0, 0, 0, 0, \
                        0, 0, 0, 0, 0, 0, \
                        0, 0, 0, 0, 0, 0};
        double delta = 0;
        double g_delta = 0;
        for (int i = 0; i < data_size; i++){
            delta = exp(- corrected_value[i] * corrected_value[i] / (2 * correntropy_factor * std_dev));
            g_delta = (model_corr_1[i]->pos_val[0] * norm_x[i] + model_corr_1[i]->pos_val[1] * norm_y[i]) * delta;
            //g = [nx⋅px⋅(mx⋅nx + my⋅ny)  nx⋅py⋅(mx⋅nx + my⋅ny)  ny⋅px⋅(mx⋅nx + my⋅ny)  ny⋅py⋅(mx⋅nx + my⋅ny)  nx⋅(mx⋅nx + my⋅ny)  ny⋅(mx⋅nx+ my⋅ny)]
            g[0] += norm_x[i] * transformed_set.x_vals[i] * g_delta;
            g[1] += norm_x[i] * transformed_set.y_vals[i] * g_delta;
            g[2] += norm_y[i] * transformed_set.x_vals[i] * g_delta;
            g[3] += norm_y[i] * transformed_set.y_vals[i] * g_delta;
            g[4] += norm_x[i] * g_delta;
            g[5] += norm_y[i] * g_delta;

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

            G[0] += norm_x[i] * norm_x[i] * transformed_set.x_vals[i] * transformed_set.x_vals[i] * delta;
            G[7] += norm_x[i] * norm_x[i] * transformed_set.y_vals[i] * transformed_set.y_vals[i] * delta;
            G[14] += norm_y[i] * norm_y[i] * transformed_set.x_vals[i] * transformed_set.x_vals[i] * delta;
            G[21] += norm_y[i] * norm_y[i] * transformed_set.y_vals[i] * transformed_set.y_vals[i] * delta;
            G[28] += norm_x[i] * norm_x[i] * delta;
            G[35] += norm_y[i] * norm_y[i] * delta;

            G[1] += norm_x[i] * norm_x[i] * transformed_set.x_vals[i] * transformed_set.y_vals[i] * delta; G[6] = G[1];//xx xy
            G[2] += norm_x[i] * norm_y[i] * transformed_set.x_vals[i] * transformed_set.x_vals[i] * delta; G[12] = G[2];//xy xx
            G[3] += norm_x[i] * norm_y[i] * transformed_set.x_vals[i] * transformed_set.y_vals[i] * delta; G[18] = G[3];//xy xy
            G[4] += norm_x[i] * norm_x[i] * transformed_set.x_vals[i] * delta;                              G[24] = G[4];//xx x
            G[5] += norm_x[i] * norm_y[i] * transformed_set.x_vals[i] * delta;                              G[30] = G[5];//xy x

            //G[8]           -- This is equal to G[3] so we don't compute it here --                                      //xy xy
            G[9] += norm_x[i] * norm_y[i] * transformed_set.y_vals[i] * transformed_set.y_vals[i] * delta; G[19] = G[9];//xy yy
            G[10] += norm_x[i] * norm_x[i] * transformed_set.y_vals[i] * delta;                             G[25] = G[10];//xx y
            G[11] += norm_x[i] * norm_y[i] * transformed_set.y_vals[i] * delta;                             G[31] = G[11];//xy y

            G[15] += norm_y[i] * norm_y[i] * transformed_set.x_vals[i] * transformed_set.y_vals[i] * delta; G[20] = G[15];//yy xy
            //G[16]          -- This is equal to G[5] so we don't compute it here --                         //xy x
            G[17] += norm_y[i] * norm_y[i] * transformed_set.x_vals[i] * delta;                             G[32] = G[17];//yy x
        
            //G[22]          -- This is equal to G[11] so we don't compute it here --                         //xy y
            G[23] += norm_y[i] * norm_y[i] * transformed_set.y_vals[i] * delta;                             G[33] = G[23];//yy y
        
            G[29] += norm_x[i] * norm_y[i] * delta;                                                          G[34] = G[29];//xy

        }
        //Set values that reappear
        G[8] = G[3]; G[13] = G[3]; G[18] = G[3];
        G[16] = G[5]; G[26] = G[5]; G[30] = G[5];
        G[22] = G[11]; G[27] = G[11]; G[31] = G[11];


        //  --  This section works fine  --

        //Print matrix so it can be easily imported into matlab if need be;
        /* std::cout << std::setprecision(30);
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

        std::cout << "g=[" << g[0] << ',' << g[1] << ',' << g[2] << ',' << g[3] << ',' << g[4] << ',' << g[5] << "]\n\n"; */

        double x[6] = {0, 0, 0, 0, 0, 0};

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
        /* std::cout << "Solution vector is: ";
        for (int i = 0; i < 6; i++){
            std::cout << x[i] << '\t';
        }
        std::cout << "\n\n"; */

        //Evaluate difference and update guess transform
        transform_diff = 0;
        transform_diff += (guess_transform.rot_mat[0] - x[0]) * (guess_transform.rot_mat[0] - x[0]);
        transform_diff += (guess_transform.rot_mat[1] - x[1]) * (guess_transform.rot_mat[1] - x[1]);
        transform_diff += (guess_transform.rot_mat[2] - x[2]) * (guess_transform.rot_mat[2] - x[2]);
        transform_diff += (guess_transform.rot_mat[3] - x[3]) * (guess_transform.rot_mat[3] - x[3]);
        transform_diff += (guess_transform.trans_vec[0] - x[4]) * (guess_transform.trans_vec[0] - x[4]);
        transform_diff += (guess_transform.trans_vec[1] - x[5]) * (guess_transform.trans_vec[1] - x[5]);
        transform_diff = transform_diff/6;

        guess_transform.rot_mat[0] = x[0];
        guess_transform.rot_mat[1] = x[1];
        guess_transform.rot_mat[2] = x[2];
        guess_transform.rot_mat[3] = x[3];
        guess_transform.trans_vec[0] = x[4];
        guess_transform.trans_vec[1] = x[5];
        std::cout << "Guessed translation of: " << guess_transform.trans_vec[0] << ", " << guess_transform.trans_vec[1] << '\n';
        guesses++;
    }
    plot.plot_data();
}