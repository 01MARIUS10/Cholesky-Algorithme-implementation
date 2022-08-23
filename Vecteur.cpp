#include <vector>
#include "Vecteur.hpp"
#include "File.hpp"

using namespace std;
Vecteur::Vecteur(int n){
    for(int i=0;i<n;i++){
        b.push_back(0);
    }
}
Vecteur::Vecteur(vector<double> v){
    b = v;
};
Vecteur::Vecteur(File f,int n){
   setDimension(n);
    b = f.getVectorFromFile(lenth);
}
Vecteur::~Vecteur(){};
vector<double> Vecteur::getVecteur(){
    return b;
};
void Vecteur::setVecteur(int i,double val){
    b[i] = val; 
};
void Vecteur::setDimension(int n){
    lenth = n;
};
void Vecteur::display(){
    cout<<"----Display vector --------"<<endl;
    cout<<"( ";
    for(int i=0;i<lenth;i++){
        cout<<b[i]<<" ";
    }
    cout<<")"<<endl;
}