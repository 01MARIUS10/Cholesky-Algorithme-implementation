#include "Matrice.hpp"
#include "Vecteur.hpp"
#include "File.hpp"
#include <iomanip>
#include "math.h"

using namespace std;

/*
@param Matrice A
permettant de trouver la Matrice triangle inferieur B tel que B.Bt = A
avec Bt ,la transposer de B (triangulaire superieur)
*/
Matrice cholesky(Matrice A){
    int n =(int)A.getMatrice().size();
    Matrice B(n);
    double s = 0;
    double res = 0;

    for(int i=0 ;i<n;i++){
        for(int j= 0;j<=i;j++){
            if(i == j){
               for(int k=0;k<i;k++){
                 s += B.getMatrice()[i][k]*B.getMatrice()[i][k];
               } 
               res = sqrt(A.getMatrice()[i][i] - s);
               B.setMatrice(i,j,res);
               res = 0;
               s=0;
            }
            else{
                for(int k=0;k<j;k++){
                    s += B.getMatrice()[i][k]*B.getMatrice()[j][k];
                } 
               res = (A.getMatrice()[i][j] - s) / B.getMatrice()[j][j];
               B.setMatrice(i,j,res);
               res = 0;
               s=0;
            }
        }
    }
    return B;
}

/*
@param Matrice B et vecteur b
resout le systeme d'equation B.Y=b avec B matrice triangulaire inferieur
retour le vecteur solution Y
*/
Vecteur solve_TriangleInf(Matrice B,Vecteur b){
    int n= (int)B.getMatrice().size();
    Vecteur y(n);
    double s = 0;
    double val = 0;
    for(int i=0;i<n;i++){
		s = 0;
		for(int j=0;j<i;j++){
			s += B.getMatrice()[i][j] * y.getVecteur()[j];
		}		
 		val = ((b.getVecteur()[i]-s)/B.getMatrice()[i][i]);
        y.setVecteur(i,val);
        val=0;
	}
    return y;
}
/*
@param Matrice Bt et vecteur y
resout le systeme d'equation Bt.X=Y avec B matrice triangulaire superieur
retour le vecteur solution X
*/
Vecteur solve_TriangleSup(Matrice Bt , Vecteur y){
    int n= (int)Bt.getMatrice().size();
    Vecteur x(n);
    double s = 0;
    double val = 0;

    for(int i= n-1;i>=0;i--){
		s=0;
		for(int j=i;j<n;j++){
			s += Bt.getMatrice()[i][j] * x.getVecteur()[j];
		}		
 		val = ((y.getVecteur()[i]-s)/Bt.getMatrice()[i][i]);
        x.setVecteur(i,val);
        val =0;
 	}
    return x;
}
/*
@param Matrice matrice et vecteur vect
rendent une meilleur affichage de l'equation matrice.X=vect
*/
void EquationDisplay(Matrice matrice , Vecteur vect){
    cout<<endl<<"-------AFFICHAGE DE L' EQUATION'----------"<<endl;
    int dim = (int)matrice.getMatrice().size();
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            cout<< matrice.getMatrice()[i][j] <<"  X"<<j<<setw(10);
        }
        cout<<setw(8)<<"= "<<setw(3)<<vect.getVecteur()[i]<<endl;
    }
    cout<<"----------------------------------------"<<endl;
};

int main(){
    cout<<"Veuillez bien verifier que les donnees de la matrice soit dans le matrice.txt "<<endl;
    cout<<"Et aussi que les donnees de la matrice colonne du second membre soit dans le vecteur.txt "<<endl;
    
    File f1("matrice.txt");
    File f2("vecteur.txt");
    int n;
    cout<<"Entrer la taille du systeme d'equation = ";
    cin>>n;
    Matrice A(f1,n);
    Vecteur b(f2,n);

    EquationDisplay(A,b);

    Matrice B = cholesky(A);
    Matrice Bt = B.transposer();
    cout<<endl<<"La matrice B : "<<endl;B.display();
    cout<<"La matrice Bt reansposer de B est  : "<<endl;Bt.display();
    cout<<endl<<"Solution : "<<endl;
    solve_TriangleSup(Bt, solve_TriangleInf(B,b)).display();
    return 0;
}