#ifndef MATRICE_HPP
#define MATRICE_HPP

#include <iostream>
#include <vector>
#include "File.hpp"

using namespace std;

class Matrice {
    private :
        vector<vector<double>> matrice;
        int dim;
    public :
        Matrice(int n);
        Matrice(vector<vector<double>>M);
        Matrice(File f,int n);
        ~Matrice();
        vector<vector<double>> getMatrice();
        void setMatrice(int i,int j,double val);
        void setDimension(int n);
        Matrice transposer();
        void display();
};
#endif
