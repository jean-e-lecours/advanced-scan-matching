#include "../hdr/data.hpp"
#include "../hdr/asm_extra.hpp"

double Transform2D::compare_transform(double* vector){
    double transform_diff = 0;
    transform_diff += (this->rot_mat[0] - vector[0]) * (this->rot_mat[0] - vector[0]);
    transform_diff += (this->rot_mat[1] - vector[1]) * (this->rot_mat[1] - vector[1]);
    transform_diff += (this->rot_mat[2] - vector[2]) * (this->rot_mat[2] - vector[2]);
    transform_diff += (this->rot_mat[3] - vector[3]) * (this->rot_mat[3] - vector[3]);
    transform_diff += (this->trans_vec[0] - vector[4]) * (this->trans_vec[0] - vector[4]);
    transform_diff += (this->trans_vec[1] - vector[5]) * (this->trans_vec[1] - vector[5]);
    return transform_diff/6;
}
void Transform2D::update_transform(double *vector){
    this->rot_mat[0] = vector[0];
    this->rot_mat[1] = vector[1];
    this->rot_mat[2] = vector[2];
    this->rot_mat[3] = vector[3];
    this->trans_vec[0] = vector[4];
    this->trans_vec[1] = vector[5];
}
void Transform2D::print_transform(){
    std::cout << "⎡" << rot_mat[0] << "\t " << rot_mat[1] << "\t " << trans_vec[0] << "⎤\n" \
              << "|" << rot_mat[2] << "\t " << rot_mat[3] << "\t " << trans_vec[1] << "|\n" \
              << "⎣0\t 0\t 1⎦\n";
}
void Transform2D::transform_set(Set *set){
    for (int i = 0; i < set->data_size; i++){
        //Applies a single matrix multiplication with the provided x, y vector from the set
        double temp_x = rot_mat[0] * set->points[i].val[0] + rot_mat[1] * set->points[i].val[1] + trans_vec[0];
        double temp_y = rot_mat[2] * set->points[i].val[0] + rot_mat[3] * set->points[i].val[1] + trans_vec[1];
        set->points[i].val[0] = temp_x;
        set->points[i].val[1] = temp_y;
    }
}
void Transform2D::add_transform(Transform2D *added_transform){
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
Transform2D::Transform2D(double x_trans, double y_trans, double z_rot){
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

Correlation::Correlation(KdTree* model_tree, Set* data_set, int index, Transform2D* transform){
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
Correlation::Correlation(){
}

int TextLaserScanData::read_from_file(std::string filename){
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
        if (las_val < 5){
            usable_las_size++;
        }
                
    }
    std::cout << "done reading\n";
    ilasf.close();
    las_size = las_vec.size();
    return 1;
}
Set TextLaserScanData::map_scan_points(Transform2D *transform, double scan_period){
    double* las_x = new double[usable_las_size];
    double* las_y = new double[usable_las_size];
    int j = 0;
    bool skip = true;
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
void TextLaserScanData::filter_scan_points(){
    for (int i = 1; i < las_vec.size(); i++){
        if (std::abs(las_vec[i] - std::abs(las_vec[i-1])) < 0.05){
            las_vec[i] *= -1;
            usable_las_size--;
        }
    }
}
TextLaserScanData::TextLaserScanData(bool file_permanence){
    this->file_permanence = file_permanence;
}