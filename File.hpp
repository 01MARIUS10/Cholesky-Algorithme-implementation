#ifndef FILE_HPP
#define FILE_HPP

#include <string>
#include <iostream>
#include <vector>

using namespace std;

class File {
    private:
        string name;
    public :
        File();
        File(string name);
        ~File();
        vector<vector<double>>getMatricefromFile(int dim);
        vector<double> getVectorFromFile(int dim);
};
#endif