#ifndef ASM_H
#define ASM_H

#include <iostream>
#include <vector>
#include <fstream>

#include "kdtree.hpp"

class Transform2D{
    public:
        double rot_mat[4];
        double trans_vec[2];
        double z_rot;

        double compare_transform(double* vector);

        void update_transform(double* vector);

        void print_transform();

        void transform_set(Set* set);

        Point transform_referenced_point(double x, double y){
            double temp_x = rot_mat[0] * x + rot_mat[1] * y + trans_vec[0];
            double temp_y = rot_mat[2] * x + rot_mat[3] * y + trans_vec[1];
            Point point(2);
            point.val[0] = temp_x;
            point.val[1] = temp_y;
            return point;
        }
        void add_transform(Transform2D* added_transform);

        Transform2D(double x_trans, double y_trans, double z_rot);
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

        Correlation(KdTree* model_tree, Set* data_set, int index, Transform2D* transform);
        Correlation();
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

        int read_from_file(std::string filename);
        
        Set map_scan_points(Transform2D* transform, double scan_period);

        void filter_scan_points();

        TextLaserScanData(bool file_permanence);
};

#endif