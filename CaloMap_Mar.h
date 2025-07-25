#ifndef CALO_MAP_H
#define CALO_MAP_H

#include <map>
#include <vector>
using namespace std;

map<pair<int,int>, vector<int>> GetCaloChMap(void)
{
    map<pair<int,int>, vector<int>> chMap;

    // 각 vector<int> 구성:
    // [0] L/R: 0 = Left, 1 = Right
    // [1] Module ID
    // [2] Column (0~7)
    // [3] Row (0~3)

    // MID 41
    chMap[{41, 1}] = {0, 1, 7, 0};
    chMap[{41, 2}] = {1, 1, 7, 0};
    chMap[{41, 3}] = {0, 6, 7, 1};
    chMap[{41, 4}] = {1, 6, 7, 1};
    chMap[{41, 5}] = {1, 5, 4, 0};
    chMap[{41, 6}] = {1, 26, 4, 1};
    chMap[{41, 7}] = {1, 30, 4, 2};
    chMap[{41, 8}] = {1, 29, 4, 3};
    chMap[{41, 9}] = {1, 3, 5, 0};
    chMap[{41,10}] = {1, 18, 5, 1};
    chMap[{41,11}] = {1, 17, 5, 2};
    chMap[{41,12}] = {1, 31, 5, 3};
    chMap[{41,13}] = {1, 2, 6, 0};
    chMap[{41,14}] = {1, 27, 6, 1};
    chMap[{41,15}] = {1, 22, 6, 2};
    chMap[{41,16}] = {1, 28, 6, 3};
    chMap[{41,17}] = {0, 7, 7, 2};
    chMap[{41,18}] = {1, 7, 7, 2};
    chMap[{41,19}] = {0, 15, 7, 3};
    chMap[{41,20}] = {1, 15, 7, 3};
    chMap[{41,21}] = {0, 5, 4, 0};
    chMap[{41,22}] = {0, 26, 4, 1};
    chMap[{41,23}] = {0, 30, 4, 2};
    chMap[{41,24}] = {0, 29, 4, 3};
    chMap[{41,25}] = {0, 3, 5, 0};
    chMap[{41,26}] = {0, 18, 5, 1};
    chMap[{41,27}] = {0, 17, 5, 2};
    chMap[{41,28}] = {0, 31, 5, 3};
    chMap[{41,29}] = {0, 2, 6, 0};
    chMap[{41,30}] = {0, 27, 6, 1};
    chMap[{41,31}] = {0, 22, 6, 2};
    chMap[{41,32}] = {0, 28, 6, 3};

    // MID 42
    chMap[{42, 1}] = {1, 8, 0, 0};
    chMap[{42, 2}] = {1, 9, 0, 1};
    chMap[{42, 3}] = {1, 24, 0, 2};
    chMap[{42, 4}] = {1, 12, 0, 3};
    chMap[{42, 5}] = {1, 14, 1, 0};
    chMap[{42, 6}] = {1, 13, 1, 1};
    chMap[{42, 7}] = {1, 21, 1, 2};
    chMap[{42, 8}] = {1, 10, 1, 3};
    chMap[{42, 9}] = {1, 4, 2, 0};
    chMap[{42,10}] = {1, 25, 2, 1};
    chMap[{42,11}] = {1, 23, 2, 2};
    chMap[{42,12}] = {1, 32, 2, 3};
    chMap[{42,13}] = {1, 11, 3, 0};
    chMap[{42,14}] = {1, 33, 3, 1};
    chMap[{42,15}] = {1, 20, 3, 2};
    chMap[{42,16}] = {1, 19, 3, 3};
    chMap[{42,17}] = {0, 8, 0, 0};
    chMap[{42,18}] = {0, 9, 0, 1};
    chMap[{42,19}] = {0, 24, 0, 2};
    chMap[{42,20}] = {0, 12, 0, 3};
    chMap[{42,21}] = {0, 14, 1, 0};
    chMap[{42,22}] = {0, 13, 1, 1};
    chMap[{42,23}] = {0, 21, 1, 2};
    chMap[{42,24}] = {0, 10, 1, 3};
    chMap[{42,25}] = {0, 4, 2, 0};
    chMap[{42,26}] = {0, 25, 2, 1};
    chMap[{42,27}] = {0, 23, 2, 2};
    chMap[{42,28}] = {0, 32, 2, 3};
    chMap[{42,29}] = {0, 11, 3, 0};
    chMap[{42,30}] = {0, 33, 3, 1};
    chMap[{42,31}] = {0, 20, 3, 2};
    chMap[{42,32}] = {0, 19, 3, 3};

    return chMap;
}


#endif // CALO_MAP_H
