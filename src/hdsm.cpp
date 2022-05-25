#include "../incl/hdsm.hpp"
double find_hd(std::vector<Correlation>& corrs, Transform2D& transform, int cutoff){
    std::vector<double> distances;
    
    //Creating vector of distances to nearest neighbour
    for (int i = 0; i < corrs.size(); i++){
        //Adjusting the value of the scan point based on the transform
        Point2D t_scan(transform.rot_mat[0]*corrs[i].scan.x+transform.rot_mat[1]*corrs[i].scan.y+transform.trans_vec[0], transform.rot_mat[2]*corrs[i].scan.x+transform.rot_mat[3]*corrs[i].scan.y+transform.trans_vec[1]); 
        distances.push_back(sqrt((t_scan.x - corrs[i].mcor1.x) * (t_scan.x - corrs[i].mcor1.x) + (t_scan.y - corrs[i].mcor1.y) * (t_scan.y - corrs[i].mcor1.y)));
    }
    //Use bubble sort to sort the distances vector
    for (int i = 0; i < distances.size(); i++){   
        for (int j = i; j < distances.size(); j++){
            if (distances[i] < distances[j]){
                double temp = distances[i];
                distances[i] = distances[j];
                distances[j] = temp;
            }
        }
    }
    return distances[cutoff];
}