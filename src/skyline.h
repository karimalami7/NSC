#ifndef SKYLINE_H_INCLUDED
#define SKYLINE_H_INCLUDED

#include "declarations.h"

//calcule le score de chaque tuple d'une table de données
DataType* subspaceScore(vector<Point> &donnees, vector<Space> &subspace){
    DataType n=donnees.size();
    DataType* scores=new DataType[n];
    DataType i;
    Space j;
    for (i=0;i<n;i++){
        scores[i]=0;
        for (j=0;j<(Space)subspace.size();j++){
            scores[i]+=donnees[i][subspace[j]];
        }
    }
    return scores;
}

//retourne un vecteur de booleens valant 'true' pour les tuples skyline
bool* subspaceSkyline_NAIF(vector<Point> &donnees, vector<Space> &subspace){
    DataType n=donnees.size();
    DataType* scores= subspaceScore(donnees, subspace);
    bool* tableSkyline=new bool[n];
    DataType i, j;
    Space l;
    for (i=0;i<n;i++) tableSkyline[i]=true;

    for (i=0;i<n-1;i++){
        for (j=i+1;j<n && tableSkyline[i];j++){
            if (tableSkyline[j]){
                if (scores[i]<scores[j]){
                    l=0;
                    while(l<(Space)subspace.size() && donnees[i][subspace[l]]<=donnees[j][subspace[l]]) l++;
                    if (l==(Space)subspace.size()) tableSkyline[j]=false;
                }else if (scores[i]>scores[j]){
                    l=0;
                    while(l<(Space)subspace.size() && donnees[i][subspace[l]]>=donnees[j][subspace[l]]) l++;
                    if (l==(Space)subspace.size()) tableSkyline[i]=false;
                }

            }
        }
    }
    delete[] scores;
    return tableSkyline;
}

//retourne la taille exacte du skyline d'un jeu de données
DataType subspaceSkylineSize_NAIF(vector<Point> &donnees, vector<Space> &subspace){
    DataType skySize;
    bool* skyline=subspaceSkyline_NAIF(donnees, subspace);

    skySize=compteLesMarques(skyline, donnees.size());
    delete[] skyline;
    return skySize;
}






#endif // SKYLINE_H_INCLUDED
