#include<vector>
#include "main.h"

std::vector<coordinate> set_mesh(double x_start, double y_start, int num_nodes_x, int num_nodes_y, int nj, double h){
    std::vector<coordinate> mesh(nj);
    double x_coordinate, y_coordinate;
    int counter = 0;

    y_coordinate = y_start;
    for (int i = 0; i < num_nodes_y; ++i) {
        x_coordinate = x_start;

        for (int j = 0; j < num_nodes_x; ++j) {
            mesh[counter].y = y_coordinate;
            mesh[counter].x = x_coordinate;
            x_coordinate += h;
            counter++;
        }
        y_coordinate += h;
    }
    return mesh;
}
