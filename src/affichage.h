#ifndef AFFICHAGE_H_INCLUDED
#define AFFICHAGE_H_INCLUDED

#include "declarations.h"

void displaySubspace(Space subspace, Space d){
//écrit un sous espace sous forme alphabétique: 1->A, 12->CD
    bool table[d];
    Space j;
    Space pow=1;
    Space aux=subspace;
    for (j=0;j<d;j++) {
        table[j]=false;
        pow*=2;
    }
    for (j=d;j>=1;j--){
        pow/=2;
        if (aux>=pow){
            table[j-1]=true;
            aux-=pow;
        }
    }
    for (j=0;j<d;j++) if (table[j]) cout<<(char)(j+65);
}

void displayDualSpace(DualSpace sortie, Space d){
	cout<<"(";
	displaySubspace(sortie.dom,d);
	cout<<";";
	displaySubspace(sortie.equ,d);
	cout<<")";
}

void afficheUSetDualSpace(USetDualSpace &uSetDualSpace, Space d){
    USetDualSpace::iterator it;
    for (it=uSetDualSpace.begin();it!=uSetDualSpace.end();it++){
        displayDualSpace(*it, d);cout<<" ";
    }
    cout<<endl;
}

void afficheVectorVectorTupleSpace(VectorVectorTupleSpace &structure, Space d){
//affiche la structure de données OPTI
    DataType nb, nbvid=0;
    long nbSubspace=1<<d;
    VectorTupleSpace::iterator it;
    Space i;
    for (i=1;i<nbSubspace;i++){
        displaySubspace(i, d);
        nb=0;
        cout<<": ";
        for (it=structure[i].begin();it!=structure[i].end();it++){
            cout<<"t"<<(it->idTuple)+1<<"(";
            displaySubspace(it->eqSubspace, d);
            cout<<")  ";
            nb++;
        }
        if (nb==0) nbvid++;
        cout<<endl;
    }
   //cout<<"Nb cellules vides: "<<nbvid<<endl;
}

#endif // AFFICHAGE_H_INCLUDED
