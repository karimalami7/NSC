#ifndef DECLARATIONS_H_INCLUDED
#define DECLARATIONS_H_INCLUDED

using namespace std;
#include <iostream>
#include <list>
#include <vector>
#include <bitset>
#include <ctime>
#include <fstream>
#include <sstream>
#include <stack>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <omp.h>
#include <queue>

#define SPACESIZE	30

#define NB_METHOD 12
#define NAIF 0
#define TREE 1
#define CSC 2
#define BUA 3
#define UBA 4
#define DPI 5
#define NSC 6
#define NSC2 7
#define NSCk 8
#define NSC3 9
#define NSCv2 10
#define NSCwM 11

string methodNames[]         ={"NAIF", "TREE", "CSC" , "BUA" , "UBA" , "DPI" , "NSC" , "NSC2" , "NSCk", "NSC3", "NSCv2", "NSCwM"};
string methodNamesDisplayed[]={"NAIF", "TREE", "CSC_", "BUA_", "UBA_", "DPI_", "NSC_", "NSC2" , "NSCk", "NSC3", "NSCv2", "NSCwM"};


struct DualSpace;
struct TupleSpace;


typedef int DataType;
typedef DataType* Point;
typedef int  Space;
typedef vector<map<int,vector<int>>> values_map;
typedef unordered_set<DualSpace> USetDualSpace;
typedef unordered_map<DualSpace, int> mapDualSpace;
typedef vector<DualSpace> VectorDualSpace;
typedef vector<TupleSpace> VectorTupleSpace;
typedef vector<VectorTupleSpace> VectorVectorTupleSpace;
typedef vector<Point> TableTuple;
typedef map<int,Point> TableTuplev2;
typedef unordered_set<DataType> USetId;
typedef vector<USetId> VectorUSetId;
typedef map<Space, unordered_map<Space,vector<DataType>>> NegSkyStrAux;
typedef vector<pair<Space,vector<pair<Space,vector<DataType>>>>> NegSkyStr;
// NSCv2
typedef map<int,map<string,pair<bool,pair<DualSpace,vector<int>>>>> index_structure;


inline Space spaceSize(Space subspace);

struct DualSpace{
// couple de deux sous espaces (dom, equ), résultat de comparaison entre deux tuples t_i et t_j, où:
//  - dom est le sous-espace dans lequel le premier tuple domine strictement le second (t_i[dom]<t_j[dom])
//  - equ est le sous-espace dans lequel les deux tuples sont égaux t_i[equ]=t_j[equ]
    Space dom;
    Space equ;
    Space poids;
    bool operator==(const DualSpace &ds)const{
        return (this->dom==ds.dom && this->equ==ds.equ);
    }
    bool operator!=(const DualSpace &ds)const{
        return (this->dom!=ds.dom || this->equ!=ds.equ);
    }
    Space dsSize(Space d){
        return spaceSize(dom+equ);
    }
};

bool operator<(const DualSpace &a, const DualSpace &b){
    if (a.poids < b.poids) return true;
    if (a.dom < b.dom) return true;
    else return a.equ<b.equ;
}

namespace std {
    template <>
    struct hash<DualSpace>  {
        public:
            size_t operator()(const DualSpace &s) const{
                //polynome de Cantor f(x,y)= ((x+y)^2+3x+y)/2
                return (pow(s.dom+s.equ, 2)+3*s.dom+s.equ)/2;
            }
    };
}

struct TupleSpace{
// couple formé de l'identifiant d'un tuple idTuple et d'un sous-espace eqSubspace où il n'est pas strictement dominé
//  relativement au couple (dom, eqSubspace) qui l'a généré.
    DataType idTuple;
    Space eqSubspace;
    const bool operator<(const TupleSpace &ts2) const{
        return (this->idTuple<ts2.idTuple);
    }
};

inline double debut(){
    return omp_get_wtime();
}

inline double duree(double start){
// calcule la durée en millisecondes entre l'instant 'start' et maintenant
;
    //return (clock()-start)*(1000.0/CLOCKS_PER_SEC);
    return (omp_get_wtime()-start)*(1000.0);
}

inline bool disjoints(const Space &spc1, const Space &spc2){
/*  teste si spc1 et spc2 sont disjoints tous codés en décimal
    exemples (codage) quand d=4, ABCD->15, AD->9, A->1, B->2, C->4, D->8, BD->10
    retourne vrai s'ils sont disjoints et faux sinon                */
    return ((spc1 & spc2)==0);
}

inline bool estInclusDans(const Space &spc1, const Space &spc2){
/*  test si spc1 est inclus dans spc2 tous codés en décimal
    exemples (codage) quand d=4, ABCD->15, AD->9, A->1, B->2, C->4, D->8, BD->10
    retourne vrai si cette inclusion est vérifiée et faux sinon                  */
    return ((spc1 & spc2)==spc1);
}

inline bool estEgale(const Space &spc1, const Space &spc2){
/*  test si spc1 est inclus dans spc2 tous codés en décimal
    exemples (codage) quand d=4, ABCD->15, AD->9, A->1, B->2, C->4, D->8, BD->10
    retourne vrai si cette inclusion est vérifiée et faux sinon                  */
    return (spc1==spc2);
}

inline Space soustraction(const Space &spc1, const Space &spc2){
/*  retourne spc1\spc2                  */
    return (spc1 & ~spc2);
}

inline int quiReduitQui_DS(const DualSpace &ds1, const DualSpace &ds2){
    //teste entre deux DualSpace lequel est susceptible de simplifier l'autre (suivant la règle décrite par Sofian)
    if ((ds1.dom | ds2.dom)==0) return 0;
    if (estInclusDans(ds2.dom+ds2.equ,ds1.dom)) return 1;
    if (estInclusDans(ds1.dom+ds1.equ,ds2.dom)) return 2;
    return 0;
}

inline bool estCouvertPar(const Space &sp, const DualSpace &ds){
//teste si un sous-espace est couvert par un couple (DualSPace)
    return ((sp & (ds.dom+ds.equ))==sp && !disjoints(sp, ds.dom));
}

inline bool estCouvertPar(const DualSpace &ds1, const DualSpace &ds2){
//teste si un couple d'espaces est couvert par un autre couple (DualSPace)
    return estInclusDans(ds1.dom+ds1.equ, ds2.dom+ds2.equ) && estInclusDans(ds1.dom, ds2.dom);
}

inline void sortIndexes(const vector<DataType> &tab, vector<DataType> &index) {
//Cette procédure trie le vecteur index suivant l'ordre de leur valeur associées contenues dans le vecteur tab
  sort(index.begin(), index.end(), [&tab](DataType i1, DataType i2) {return tab[i1] < tab[i2];});
}

inline void sortIndexes(const vector<double> &tab, vector<DataType> &index) {
//Cette procédure trie le vecteur index suivant l'ordre de leur valeur associées contenues dans le vecteur tab
  sort(index.begin(), index.end(), [&tab](DataType i1, DataType i2) {return tab[i1] < tab[i2];});
}

inline vector<DataType>::iterator maxIndexes(const vector<DataType> &tab, vector<DataType> &index) {
//Cette procédure trie le vecteur index suivant l'ordre de leur valeur associées contenues dans le vecteur tab
  return max_element(index.begin(), index.end(), [&tab](DataType i1, DataType i2) {return tab[i1] < tab[i2];});
}

DataType compteLesMarques(bool* table, DataType length){
/*  function qui compte le nombre de cellules valant 'vrai' dans une table de booleens
    utile pour calculer la taille du skyline*/
    DataType i, total=0;
    for (i=0; i<length; i++) total+=table[i]?1:0;
    return total;
}

void listeAttributsPresents(Space subspace, Space d, vector<Space> &result){
/*  retourne la liste des attributs d'un sous espace
    le résultat se présente sous la forme d'un vector
    où chaque cellule correspond au numéro de l'attribut 1, 2, 3, ...,d
    la première valeur [0] est le nombre d'attributs */
    vector<Space> sortie;
	bitset<SPACESIZE> subspace_aux=subspace;
	Space j;
    for (j=0;j<d;j++) if (subspace_aux[j]) sortie.push_back(j+1);
    result.swap(sortie);
}

void listeAttributsManquants(Space subspace, Space d, vector<Space> &result){
/*  retourne la liste des attributs absents dans un sous espace
    le résultat se présente sous la forme d'un vector
    où chaque cellule correspond au numéro de l'attribut 1, 2, 3, ..., d */
    vector<Space> sortie;
	bitset<SPACESIZE> subspace_aux=subspace;
	Space j;
    for (j=0;j<d;j++) if (!subspace_aux[j]) sortie.push_back(j+1);
    result.swap(sortie);
}

inline Space spaceSize(const Space subspace){
/*  retourne la taille d'un sous-espace (nombre d'attributs */
    return __builtin_popcount(subspace);
}

void listSousEspacesR(vector<Space> &listSubspace, Space actualSpace, vector<Space> &missingAttrib, Space nbDiscarded){
//  aide pour la recherche des sous-espaces d'un espace donné.
//  doit être appelée par la procédure ci-dessous
    Space i, newActual;
    for (i=nbDiscarded;i<(Space)missingAttrib.size();i++){
        newActual=actualSpace+(1<<(missingAttrib[i]-1));
        listSubspace.push_back(newActual);
        listSousEspacesR(listSubspace, newActual, missingAttrib, i+1);
    }
}

void inline listSousEspaces(vector<Space> &listSubspace, Space &originSpace, const Space &d){
//  Cette procédure liste de façon récursive les sous-espaces d'un espace donné
//  appelle la procédure ci-dessus
    vector<Space> listSubspaceR;
    vector<Space> missingAttrib;
    listeAttributsPresents(originSpace, d, missingAttrib);
    listSousEspacesR(listSubspaceR, 0, missingAttrib, 0);
    listSubspace.swap(listSubspaceR);
}


void inline listCouvertsR(Space actualSpace, vector<Space> &missingAttrib, Space nbDiscarded, vector<Space> &listEquSubspace, vector<Space> &listSpace){
//  aide à retrouver les espacescouverts par un DualSPace donné
//  doit être appelé par la procédure ci-dessous
    Space newActual;
    Space i, j;
    for (i=nbDiscarded;i<(Space)missingAttrib.size();i++){
        newActual=actualSpace+(1<<(missingAttrib[i]-1));
        listSpace.push_back(newActual);
        for (j=0;j<(Space)listEquSubspace.size();j++){
            listSpace.push_back(newActual+listEquSubspace[j]);
        }
        listCouvertsR(newActual, missingAttrib, i+1, listEquSubspace, listSpace);
    }
}

void listCouverts(Space domSpace, Space equSpace, Space &d, vector<Space> &listSpace){
//  Liste l'ensemble des espaces couverts par un dualspace donné
//  appelle la procédure ci-dessus
    vector<Space> listSubspace;
    listSousEspaces(listSubspace, equSpace, d);
    vector<Space> missingAttrib;
    listeAttributsPresents(domSpace, d, missingAttrib);
	Space s1=0, s2=0;
    listCouvertsR(s1, missingAttrib, s2, listSubspace, listSpace);
}


void listCouverts2(Space domSpace, Space equSpace, Space d, vector<Space> &listSpace){
    vector<Space> listDomSubspace, listEquSubspace;
    listSousEspaces(listDomSubspace, domSpace, d);
    listSousEspaces(listEquSubspace, equSpace, d);
    listEquSubspace.push_back(0);
    for (auto itD=listDomSubspace.begin(); itD!=listDomSubspace.end();itD++){
        for (auto itE=listEquSubspace.begin(); itE!=listEquSubspace.end(); itE++){
            listSpace.push_back((*itD)+(*itE));
        }
    }
}


#endif // DECLARATIONS_H_INCLUDED
