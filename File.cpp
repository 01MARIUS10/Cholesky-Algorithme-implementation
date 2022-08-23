#include "File.hpp"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

File::File(){
    cout<<"entrer le nom du fichier : ";
    cin>>name;
};
File::File(string nom){
    name=nom;
}
File::~File(){}
vector<double> File::getVectorFromFile(int dim){
    cout<<"Getting vector from file......."<<endl;
    double doubleTmp;
    vector<double> vector;
    ifstream f (name);

    if (!f.is_open())
        cout << "Impossible d'ouvrir le fichier en lecture ! \n sachez que le vecteur colonne du second membre doit etre dans le fichier vecteur.txt sur le meme racine que le programme" << endl;
    else
    {
        for (int i = 0; i < dim; i++){
            f >> doubleTmp;
            vector.push_back(doubleTmp);
            doubleTmp=0;
        }
    }
    f.close();
    return vector;
}
vector<vector<double>> File::getMatricefromFile(int dim){
    cout<<"Getting matrice from file......."<<endl;
    double doubleTmp;
    vector<double> vectorTmp;
    vector<vector<double>> matrice ;
    ifstream f (name);

    if (!f.is_open())
        cout << "Impossible d'ouvrir le fichier en lecture ! \n sachez que la matrice doit etre dans le fichier matrice.txt sur le meme racine que le programme " << endl;
    else
    {
        for (int i = 0; i < dim; i++){
            for(int j=0;j<dim;j++){
                f >> doubleTmp;
                vectorTmp.push_back(doubleTmp);
                doubleTmp=0;
            }
            matrice.push_back(vectorTmp);
            vectorTmp.clear();
        }
    }
    f.close();
    return matrice;
};

