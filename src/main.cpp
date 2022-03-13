#include <bits/types/clock_t.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "../incl/gnuplot.h"
#include <math.h>
#include <string>
#include <vector>

#include "mat.hpp"

int dimensions;

class Point{
    public:
        double* val;

        explicit Point(int dims){
            this->val = new double[dims];
        }
        explicit Point(){
            this->val = new double[dimensions];
        }
};

class Node: public Point{
    public:
        int dim;
        int layer;
        Node* small_link;
        Node* big_link;
        Node* parent;
};

class Set{

    public:
        int data_size;
        int max_layer = 0;
        Point* points;

        Set* copy_set(){
            Set* new_set = new Set(data_size, points);
            return new_set;
        }

        void print_data(){
            for (int i = 0; i < data_size; i++){
                for (int d = 0; d < dimensions; d++){
                    std::cout << points[i].val[d];
                    if (d+1 == dimensions){
                        std::cout << '\n';
                    }
                    else{
                        std::cout << ", ";
                    }
                }
            }
        }

    Set(int data_size, double** point_vals){
        //Constructor
        this->data_size = data_size;
        this->points = new Point[data_size];

        for (int i = 0; i < data_size; i++){
            for (int d = 0; d < dimensions; d++){
                this->points[i].val[d] = point_vals[d][i];
            }
        }
    }
    Set(int data_size, Point* points){
        //Constuctor for copy only
        this->data_size = data_size;
        this->points = new Point[data_size];

        for (int i = 0; i < data_size; i++){
            for (int d = 0; d < dimensions; d++){
                this->points[i].val[d] = points[i].val[d];
            }
        }
    }
    ~Set(){
        //Deconstructor
    }
};

class KdTree{
    Node* kd_node_array;
    int node_index = 0;
    int dims;

    void populate_node(int dim, int start, int points_left, Node* this_node, Node* starting_node){
        int index = start + std::ceil(points_left/2.0)- 1;
        this_node->dim = dim;
        for (int d = 0; d < dimensions; d++){
            this_node->val[d] = points[index].val[d];
        }
        this_node->parent = starting_node;
        bool size = false;
        if (starting_node != NULL){
            this_node->layer = starting_node->layer + 1;
            //Node linking, not needed for the first node (has a NULL starting_node)
            dim = cycle_dimensions(dim, false);
            if (this_node->val[dim] < starting_node->val[dim]){
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

    int cycle_dimensions(int dim, bool dir){
        if (dir){
            dim++;
            if (dim >= dimensions){
                dim = 0;
            }
        }
        else{
            dim--;
            if (dim < 0){
                dim = dimensions - 1;
            }
        }
        return dim;
    }

    public:
        int data_size;
        int max_layer = 0;
        Point* points;
        int dim = 0;

        void print_kd_tree(){
            for (int i = 0; i < data_size; i++) {
                std::cout << kd_node_array[i].val[0] << ", " << kd_node_array[i].val[1];
                std::cout << "[Layer " << kd_node_array[i].layer << "] ";
                if (kd_node_array[i].big_link != NULL){
                    std::cout << " =+=> " << kd_node_array[i].big_link->val[0] << ", " << kd_node_array[i].big_link->val[1];
                }
                if (kd_node_array[i].small_link != NULL){
                    std::cout << " =-=> " << kd_node_array[i].small_link->val[0] << ", " << kd_node_array[i].small_link->val[1];
                }
                std::cout << '\n';
            }
        }

        Node* find_approximate_closest_point(double* target_point, int dim){
            Node* node_ptr = NULL;
            Node* next_node_ptr = &kd_node_array[0];

            while (next_node_ptr != NULL){
                node_ptr = next_node_ptr;
                if (target_point[0] < node_ptr->val[dim]){
                    next_node_ptr = node_ptr->small_link;
                }
                else{
                    next_node_ptr = node_ptr->big_link;
                }
                dim = cycle_dimensions(dim, true);
            }
            return node_ptr;
        }

        Node** find_closest_point(Point* target_point, int dim, bool find_second_closest){
            Node** possible_node_ptr = new Node*[(max_layer * (max_layer + 1))/2];
            Node** best_node_ptr = new Node*[2];
            best_node_ptr[0] = NULL;
            best_node_ptr[1] = NULL;
            int possible_node_index = 0;
            Node* node_ptr = &kd_node_array[0];
            double distance = 0;
            double border_distance = 0;
            double closest_distance = INFINITY;
            double second_closest_distance = INFINITY;

            while (true){
                while (node_ptr != NULL){
                    //We only go forwards and end when the pointer dies, unless there is
                    //another branch to inspect.
                    distance = sqrt((target_point->val[0] - node_ptr->val[0]) * (target_point->val[0] - node_ptr->val[0]) + (target_point->val[1] - node_ptr->val[1]) * (target_point->val[1] - node_ptr->val[1]));
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

                    dim = node_ptr->dim;
                    border_distance = target_point->val[dim] - node_ptr->val[dim];
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

        void make_kd_tree(int dim, int start, int points_left, Node* starting_node){
            for (int i = start; i < (points_left + start); i++){
                
                for (int j = i; j < (points_left + start); j++){
                    if (points[i].val[dim] > points[j].val[dim]){
                        double* temp = new double[dimensions];
                        for (int d = 0; d < dimensions; d++){
                            temp[d] = points[i].val[d];
                        }
                        for (int d = 0; d < dimensions; d++){
                            points[i].val[d] = points[j].val[d];
                        }
                        for (int d = 0; d < dimensions; d++){
                            points[j].val[d] = temp[d];
                        }
                    }
                }
            }
            populate_node(dim, start, points_left, &kd_node_array[node_index], starting_node);
            Node* node = &kd_node_array[node_index];
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
                dim = cycle_dimensions(dim, true);
                make_kd_tree(dim, start + 1, 1, node);
            }
            else{
                dim = cycle_dimensions(dim, true);
                make_kd_tree(dim, start, std::ceil(points_left/2.0 - 1), node);
                make_kd_tree(dim, start + std::ceil(points_left/2.0), std::floor(points_left/2.0), node);
                return;
            }
        }

    KdTree(Set* parent_set){
        //Constructor
        this->data_size = parent_set->data_size;
        this->dims = dimensions-1;
        this->dim = 0;
        this->points = new Point[data_size];

        for (int i = 0; i < data_size; i++){
            for (int d = 0; d < dimensions; d++){
                this->points[i].val[d] = parent_set->points[i].val[d];
            }
        }

        this->kd_node_array = new Node[data_size];

        //Calls a recursive function to make create the object and all its properties
        std::cout << "making tree\n";
        make_kd_tree(false, 0, data_size, NULL);
    }
    ~KdTree(){
        //Destructor
    }
};

class Transform2D{
    public:
        double rot_mat[4];
        double trans_vec[2];
        double z_rot;

        double compare_transform(double* vector){
            double transform_diff = 0;
            transform_diff += (this->rot_mat[0] - vector[0]) * (this->rot_mat[0] - vector[0]);
            transform_diff += (this->rot_mat[1] - vector[1]) * (this->rot_mat[1] - vector[1]);
            transform_diff += (this->rot_mat[2] - vector[2]) * (this->rot_mat[2] - vector[2]);
            transform_diff += (this->rot_mat[3] - vector[3]) * (this->rot_mat[3] - vector[3]);
            transform_diff += (this->trans_vec[0] - vector[4]) * (this->trans_vec[0] - vector[4]);
            transform_diff += (this->trans_vec[1] - vector[5]) * (this->trans_vec[1] - vector[5]);
            return transform_diff/6;
        }

        void update_transform(double* vector){
            this->rot_mat[0] = vector[0];
            this->rot_mat[1] = vector[1];
            this->rot_mat[2] = vector[2];
            this->rot_mat[3] = vector[3];
            this->trans_vec[0] = vector[4];
            this->trans_vec[1] = vector[5];
        }

        void print_transform(){
            std::cout << "⎡" << rot_mat[0] << "\t " << rot_mat[1] << "\t " << trans_vec[0] << "⎤\n" \
                      << "|" << rot_mat[2] << "\t " << rot_mat[3] << "\t " << trans_vec[1] << "|\n" \
                      << "⎣0\t 0\t 1⎦\n";
        }

        void transform_set(Set* set){
            for (int i = 0; i < set->data_size; i++){
                //Applies a single matrix multiplication with the provided x, y vector from the set
                double temp_x = rot_mat[0] * set->points[i].val[0] + rot_mat[1] * set->points[i].val[1] + trans_vec[0];
                double temp_y = rot_mat[2] * set->points[i].val[0] + rot_mat[3] * set->points[i].val[1] + trans_vec[1];
                set->points[i].val[0] = temp_x;
                set->points[i].val[1] = temp_y;
            }
        }
        Point transform_referenced_point(double x, double y){
            double temp_x = rot_mat[0] * x + rot_mat[1] * y + trans_vec[0];
            double temp_y = rot_mat[2] * x + rot_mat[3] * y + trans_vec[1];
            Point point(2);
            point.val[0] = temp_x;
            point.val[1] = temp_y;
            return point;
        }
        void add_transform(Transform2D* added_transform){
            Transform2D temp_transform(0, 0, 0);
            temp_transform.rot_mat[0] = rot_mat[0] * added_transform->rot_mat[0] + rot_mat[1] * added_transform->rot_mat[2];
            temp_transform.rot_mat[1] = rot_mat[0] * added_transform->rot_mat[1] + rot_mat[1] * added_transform->rot_mat[3];
            temp_transform.rot_mat[2] = rot_mat[2] * added_transform->rot_mat[0] + rot_mat[3] * added_transform->rot_mat[2];
            temp_transform.rot_mat[3] = rot_mat[2] * added_transform->rot_mat[1] + rot_mat[3] * added_transform->rot_mat[3];
            temp_transform.trans_vec[0] = rot_mat[0] * added_transform->trans_vec[0] + rot_mat[1] * added_transform->trans_vec[1] + trans_vec[0];
            temp_transform.trans_vec[1] = rot_mat[2] * added_transform->trans_vec[0] + rot_mat[3] * added_transform->trans_vec[1] + trans_vec[1];
            this->rot_mat[0] = temp_transform.rot_mat[0];
            this->rot_mat[1] = temp_transform.rot_mat[1];
            this->rot_mat[2] = temp_transform.rot_mat[2];
            this->rot_mat[3] = temp_transform.rot_mat[3];
            this->trans_vec[0] = temp_transform.trans_vec[0];
            this->trans_vec[1] = temp_transform.trans_vec[1];
        }

    Transform2D(double x_trans, double y_trans, double z_rot){
        //Constructor
        //Sets the two elements of the Transform2Ds as per the desired movement
        rot_mat[0] = cos(z_rot);
        rot_mat[1] = -sin(z_rot);
        rot_mat[2] = -rot_mat[1];
        rot_mat[3] = rot_mat[0];
        trans_vec[0] = x_trans;
        trans_vec[1] = y_trans;

        this->z_rot = z_rot;
    }
    ~Transform2D(){
        //Deconstructor
    }
};

class Plot2D{
    bool data_permanence = false;
    int max_plots = 0;
    int number_of_plots = 0;
    double padding = 1;
    double extremas[4] = {INFINITY, 0, INFINITY, 0};

    std::string gnup_line = "plot ";
    std::string filenames = "";

    public:
        int total_plots;

        void add_data(Set* set){
            //Add file name to the array and create the file itself
            std::string filename = "dat/temp2d_" + std::to_string(number_of_plots)+ ".dat";
            filenames += filename + "|";

            std::ofstream osetf{filename};
            //Go through all the data to find what is the maximum and minimum in x
            //and y so that we can size the plot window properly. This is not
            //necessary but we are going through all the data anyways so why not.
            for (int i = 0; i < set->data_size; i++){ 
                if (set->points[i].val[0] < extremas[0]){
                    extremas[0] = set->points[i].val[0]; 
                }
                else if (set->points[i].val[0] > extremas[1]){
                    extremas[1] = set->points[i].val[0];
                }
                if (set->points[i].val[1] < extremas[2]){
                    extremas[2] = set->points[i].val[1]; 
                }
                else if (set->points[i].val[1] > extremas[3]){
                    extremas[3] = set->points[i].val[1];
                }
                //This is where we actually place the data into the file
                osetf << set->points[i].val[0] << ' ' << set->points[i].val[1] << '\n';
            }
            number_of_plots++;
            gnup_line += "\"" + filename + "\",";
        }

        void add_line(Point* point1, Point* point2){
            double x1 = point1->val[0];
            double x2 = point2->val[0];
            double y1 = point1->val[1];
            double y2 = point2->val[1];

            double a = (y2 - y1)/(x2 - x1);
            double b = y1 - a*x1;

            if (a != INFINITY && b != INFINITY && a == a){
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

    Plot2D(bool keep_data, int max_guesses){
        //Constructor
        this->data_permanence = keep_data;
        this->max_plots = max_guesses;
    }
    ~Plot2D(){
        //Deconstructor, removes data files if desired
        if (!data_permanence){
            int s = 0;
            for (int i = 0; i < max_plots+1; i++){
                std::string filename = "";
                while(filenames[s] != '|'){
                    filename += filenames[s];
                    s++;
                }
                remove(filename.c_str());
                s++;
            }
        }
    }
};

class TextLaserScanData{
    std::vector<double> las_vec;
    bool done = false;
    bool file_permanence;

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
            while (!done){
                std::string las_read;
                ilasf >> las_read;
                int las_length = las_read.length()-1;
                if (las_read[las_length] == ','){
                    las_read.erase(las_read.end()-1);
                }
                else{
                    done = true;
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
            std::cout << "done reading\n";
            ilasf.close();
            las_size = las_vec.size();
            return 1;
        }
        Set map_scan_points(Transform2D* transform, double scan_period){
            double* las_x = new double[usable_las_size];
            double* las_y = new double[usable_las_size];
            int j = 0;
            for (int i = 0; i < las_size; i++){
                if (las_vec[i] < 5){
                    double temp_scan_x = las_vec[i] * cos(scan_period * i);
                    double temp_scan_y = las_vec[i] * sin(scan_period * i);
                    las_x[j] = transform->trans_vec[0] + temp_scan_x * cos(transform->z_rot) - temp_scan_y * sin(transform->z_rot);
                    las_y[j] = transform->trans_vec[1] + temp_scan_y * cos(transform->z_rot) + temp_scan_x * sin(transform->z_rot);
                    j++;
                }
            }
            double* las_ar[2] = {las_x, las_y};
            Set las_set(usable_las_size, las_ar);
            return las_set;
        }

    TextLaserScanData(bool file_permanence){
    //Constructor
        this->file_permanence = file_permanence;
    }
    ~TextLaserScanData(){
    //Destructor
    }
};

class Correlation{
    
    Point* model_corr_2;
    Node** temp_node;
    Point** temp_point;

    double norm_size;

    

    public:
        Point current_point;
        Point* model_corr_1;
        double norm_x;
        double norm_y;
        double corrected_value;

        static double std_dev;
        static double mean;

        Correlation(KdTree* model_tree, Set* data_set, int index, Transform2D* transform){
            //Making copies of individual point
            this->current_point.val[0] = data_set->points[index].val[0];
            this->current_point.val[1] = data_set->points[index].val[1];
            //Doing some weird pointer thing. Looks like this is working.
            this->temp_node = model_tree->find_closest_point(&current_point, false, true);
            model_corr_1 = new Point();
            model_corr_2 = new Point();
            this->model_corr_1->val[0] = temp_node[0]->val[0]; //mi1
            this->model_corr_1->val[1] = temp_node[0]->val[1]; //mi1
            this->model_corr_2->val[0] = temp_node[1]->val[0]; //mi2
            this->model_corr_2->val[1] = temp_node[1]->val[1]; //mi2
            
            //Sign of the norm does not matter, using cross product with [0;0;1]
            this->norm_x = model_corr_1->val[1] - model_corr_2->val[1];
            this->norm_y = - (model_corr_1->val[0] - model_corr_2->val[0]);
            this->norm_size = sqrt(norm_x*norm_x + norm_y*norm_y);
            this->norm_x = norm_x/norm_size;
            this->norm_y = norm_y/norm_size;

            //CHECK: I assume this is correct for now, to be investigated further if stuff doesn't work as expected
            this->corrected_value = norm_x*(transform->rot_mat[0]*current_point.val[0] + transform->rot_mat[1]*current_point.val[1] + transform->trans_vec[0] - model_corr_1->val[0]) + \
                                norm_y*(transform->rot_mat[2]*current_point.val[0] + transform->rot_mat[3]*current_point.val[1] + transform->trans_vec[1] - model_corr_1->val[1]);
        }
        Correlation(){
        }

};

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

int main(){
    //Parameters that would be set on ROS
    dimensions = 2;
    int max_guesses = 10;
    double correntropy_factor = 0.1;
    double transform_tresh = 0.00001;

    clock_t start, end;
    
    //Parameters we just need initialized
    int guesses = 0;
    double transform_diff = INFINITY;
    //This will serve as an initializer and a way to keep track of the total transforms
    Transform2D total_transform(0, 0, 0);
    Transform2D base_transform(0, 0, 0);

    start = clock();

    //Info we would get from the laser scan message
    double scan_period = (2*M_PI)/360;
    TextLaserScanData laser_scan_1(false);
    laser_scan_1.read_from_file("incl/test_dat/test_3_1.txt");
    TextLaserScanData laser_scan_2(false);
    laser_scan_2.read_from_file("incl/test_dat/test_3_2.txt");
    
    
    //Will need something to interpret the laser scan message
    int data_size = laser_scan_1.usable_las_size;

    Set model_set = laser_scan_1.map_scan_points(&base_transform, scan_period);
    Set data_set = laser_scan_2.map_scan_points(&base_transform, scan_period);

    Plot2D plot(true, max_guesses);

    KdTree model_tree(&model_set);

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
            //plot.add_line(corrs[i].model_corr_1, &corrs[i].current_point);
        }
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

    double* AtA = affine_to_rot(total_transform.rot_mat, dimensions, true);

    Transform2D test_transform(total_transform.trans_vec[0],total_transform.trans_vec[1],0);
    test_transform.rot_mat[0] = AtA[0];
    test_transform.rot_mat[1] = AtA[1];
    test_transform.rot_mat[2] = AtA[2];
    test_transform.rot_mat[3] = AtA[3];
    //This is a transform with a proper rotation matrix
    test_transform.print_transform();
    //TODO: adjust translation vector for the new proper rotation matrix


    std::cout << "Final transform is:\n";
    total_transform.print_transform();
    end = clock();
    double time_taken = double(end-start) / double(CLOCKS_PER_SEC);
    std::cout << "Time for execution: " << time_taken << "s\n";
    plot.plot_data();
    return 0;
}