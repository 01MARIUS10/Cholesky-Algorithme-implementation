#include <vector>
#include <iomanip>
#include "Matrice.hpp"
#include "File.hpp"

using namespace std;

Matrice::Matrice(int n){
    dim = n;
    vector<double> ligne;
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            ligne.push_back(0);
        }
        matrice.push_back(ligne);
        ligne.clear();
    }
};
Matrice::Matrice(vector<vector<double>> M){
    matrice = M;
    dim =(int)M.size();
};
Matrice::Matrice(File f,int n){
    setDimension(n);
    matrice = f.getMatricefromFile(dim);
}
Matrice::~Matrice(){}
vector<vector<double>> Matrice::getMatrice(){
    return matrice;
}
void Matrice::setMatrice(int i,int j,double val){
    matrice[i][j]=val;
};
void Matrice::setDimension(int n){
    dim = n;
}
void Matrice::display(){
    cout<<"-------AFFICHAGE DE LA MATRICE----------"<<endl;
    for(int i=0;i<dim;i++){
        cout<<"\t";
        for(int j=0;j<dim;j++){
        cout<<setw(3)<<matrice[i][j];
        }
        cout<<endl;
    }
    cout<<"----------------------------------------"<<endl;
}
Matrice Matrice::transposer(){
    Matrice T(dim);
    double val = 0;
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            val = j<i? 0:matrice[j][i];
            T.setMatrice(i,j,val);
            val = 0;
        }
    }
    return T;
};
