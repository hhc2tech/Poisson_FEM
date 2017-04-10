//
// Created by Pushkar Kumar Jain on 4/10/17.
//
#include<vector>

std::vector<int> set_bc(int num_nodes_x, int num_nodes_y){
    int counter = 0;
    std::vector<int> id(num_nodes_x*num_nodes_y);

    for (int i = 0; i < num_nodes_y; ++i) {
        for (int j = 0; j < num_nodes_x; ++j) {
            (i == 0 || i == num_nodes_y - 1 || j == 0 || j == num_nodes_x - 1) ? id[counter] = 1 : id[counter] = 0;
            counter++;
        }
    }
    return id;
}

