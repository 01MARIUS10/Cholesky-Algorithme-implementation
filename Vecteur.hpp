#ifndef VECTEUR_HPP
#define VECTEUR_HPP

#include <iostream>
#include <vector>
#include "File.hpp"

using namespace std;

class Vecteur {
    private :
        vector<double> b;
        int lenth;
    public :
        Vecteur(int n);
        Vecteur(vector<double> v);
        Vecteur(File f,int n);
        ~Vecteur();
        vector<double>getVecteur();
        void setVecteur(int i,double val);
        void setDimension(int n);
        void display();
};

#endif