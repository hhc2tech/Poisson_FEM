//
// Created by Pushkar Kumar Jain on 4/10/17.
//
#include<vector>
#include "main.h"

std::vector<connectivity> set_connectivity(int nelem, int num_nodes_x){
    std::vector<connectivity> inode(nelem);
    int iter1 = 0;
    int iter2 = 0;
    for (int element_num = 0; element_num < nelem; ++element_num) {
        inode[element_num].node[0] = iter1 + iter2 * (num_nodes_x);
        inode[element_num].node[1] = inode[element_num].node[0] + 1;
        inode[element_num].node[2] = inode[element_num].node[1] + num_nodes_x;
        inode[element_num].node[3] = inode[element_num].node[2] - 1;
        iter1++;
        if ((element_num + 1) % (num_nodes_x - 1) == 0) {
            iter1 = 0;
            iter2++;
        }
    }
    return inode;
}


