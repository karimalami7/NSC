#ifndef KDOMSKYCUBE_H_INCLUDED
#define KDOMSKYCUBE_H_INCLUDED

#include "../common/declarations.h"

/*
Space kCompare(const Point &t1, const Point &t2, const Space &d){
//  retourne la taille du sous-espace maximal dans lequel t1 domine t2
    Space dom=0;
    Space equ=0;
	for(auto j = 1; j <= d ; j++){
		if(t1[j] < t2[j]) {dom++;}
		else if(t1[j] == t2[j]) {equ++;}
	}
	if (dom>0) return dom+equ;
	else return 0;
}

void bua_kp1(TableTuple &donnees, TableTuple &Rk, TableTuple &Tk, DataType n, Space d, Space k){
    TableTuple Tkp1;
    if (Rk.size()==n) return;
    for (auto i=0;i<Tk.size();i++){
        bool kDominated=false;
        for (auto j=0;!kDominated && j<donnees.size();j++){
            if (kCompare(donnees[j], Tk[i], d)>=k) kDominated=true;
        }
        if (kDominated) Tkp1.push_back(Tk[i]);
        else Rk.push_back(Tk[i]);
    }
    Tk.swap(Tkp1);
}

DataType bua_cube_size(TableTuple &donnees, Space d){
    DataType n=donnees.size(), cubeSize=0;
    TableTuple Rk, Tk=donnees;
    Space k;
    for (k=1;k<=d;k++){
        bua_kp1(donnees, Rk, Tk, n, d, k);
        cubeSize+=Rk.size();
    }
    return cubeSize;
}

void uba_km1(TableTuple &donnees, TableTuple &Rk, DataType n, Space d, Space k){
    TableTuple Rkp1;
    if (Rk.size()==0) return;
    for (auto i=0;i<Rk.size();i++){
        bool kDominated=false;
        for (auto j=0;!kDominated && j<donnees.size();j++){
            if (kCompare(donnees[j], Rk[i], d)>=k) kDominated=true;
        }
        if (!kDominated) Rkp1.push_back(Rk[i]);
    }
    Rk.swap(Rkp1);
}

DataType uba_cube_size(TableTuple &donnees, Space d){
    DataType n=donnees.size(), cubeSize=0;
    TableTuple Rk=donnees;
    Space k;
    for (k=d;k>=1;k--){
        uba_km1(donnees, Rk, n, d, k);
        cubeSize+=Rk.size();
    }
    return cubeSize;
}*/



Space kCompare(const Point &t1, const Point &t2, const vector<Space> &att){
//  retourne la taille du sous-espace (de l'espace formé par les attributs contenus dans att) maximal dans lequel t1 domine t2
    Space dom=0;
    Space equ=0;
	for(auto j = 0; j <att.size()  ; j++){
		if(t1[att[j]] < t2[att[j]]) {dom++;}
		else if(t1[att[j]] == t2[att[j]]) {equ++;}
	}
	if (dom>0) return dom+equ;
	else return 0;
}

void BUA_kp1(TableTuple &donnees, TableTuple &Rk, TableTuple &Tk, DataType n, vector<Space> &att, Space k){
    TableTuple Tkp1;
    if (Rk.size()==n) return;
    for (auto i=0;i<Tk.size();i++){
        bool kDominated=false;
        for (auto j=0;!kDominated && j<donnees.size();j++){
            if (kCompare(donnees[j], Tk[i], att)>=k) kDominated=true;
        }
        if (kDominated) Tkp1.push_back(Tk[i]);
        else Rk.push_back(Tk[i]);
    }
    Tk.swap(Tkp1);
}

DataType BUA_SkylineCubeSize(TableTuple &donnees, vector<Space> &att){
    DataType n=donnees.size(), cubeSize=0;
    TableTuple Rk, Tk=donnees;
    Space k;
    for (k=1;k<=att.size();k++){
        BUA_kp1(donnees, Rk, Tk, n, att, k);
        cubeSize+=Rk.size();
    }
    return cubeSize;
}

void UBA_km1(TableTuple &donnees, TableTuple &Rk, DataType n, vector<Space> &att, Space k){
    TableTuple Rkp1;
    if (Rk.size()==0) return;
    for (auto i=0;i<Rk.size();i++){
        bool kDominated=false;
        for (auto j=0;!kDominated && j<donnees.size();j++){
            if (kCompare(donnees[j], Rk[i], att)>=k) kDominated=true;
        }
        if (!kDominated) Rkp1.push_back(Rk[i]);
    }
    Rk.swap(Rkp1);
}

DataType UBA_SkylineCubeSize(TableTuple &donnees, vector<Space> &att){
    DataType n=donnees.size(), cubeSize=0;
    TableTuple Rk=donnees;
    Space k;
    for (k=att.size();k>=1;k--){
        UBA_km1(donnees, Rk, n, att, k);
        cubeSize+=Rk.size();
    }
    return cubeSize;
}



void triDpScore(TableTuple &donnees, vector<Space> &att){
//  Cette procédure trie les données suivant le score (somme des valeurs des attributs)
    DataType i, n=donnees.size(), minValue;
    Space j;
    vector<DataType> sco(n);
    vector<double> dp(n);
    TableTuple auxT=donnees;

    for (i=0;i<n;i++) dp[i]=0;
    for (j=0;j<att.size();j++){
        minValue=donnees[0][att[j]];
        for (i=1;i<n;i++){
            if (donnees[i][att[j]]<minValue) minValue=donnees[i][att[j]];
        }
        for (i=0;i<n;i++){
            if (minValue==donnees[i][att[j]]) dp[i]++;
        }
    }

    for (i=0;i<n;i++){
        sco[i]=donnees[i][att[0]];
        for (j=1;j<att.size();j++) sco[i]+=donnees[i][att[j]];
    }

    for (i=0;i<n;i++) dp[i]=1.0/(1.0/(dp[i]+sco[i]+1));

    vector<DataType> index(donnees.size());
    for (i = 0; i != (DataType)index.size(); ++i) index[i] = i;
    sortIndexes(dp, index);
    for (i=0;i<n;i++){
        donnees[i]=auxT[index[i]];
    }
}

void DPI_kDom(TableTuple &donnees, vector<Space> &att, Space k, TableTuple &result){
    TableTuple kDomSky(0);
    triDpScore(donnees, att);
    bool isDominant;
    Space k1;
    DataType j;
    vector<Point>::iterator it;
    for (auto i=0;i<donnees.size();i++){
        isDominant = true;
        it=kDomSky.begin();
        while (it!=kDomSky.end() && isDominant){
            k1 = kCompare(*it, donnees[i], att);
            if (k1>=k){
                isDominant = false;
                it++;
            }else{
                k1 = kCompare(donnees[i], *it, att);
                if (k1>=k){
                    it=kDomSky.erase(it);
                }else it++;
            }
        }
        if (isDominant){
            kDomSky.push_back(donnees[i]);
        }
    }
    for (auto i=0;i<donnees.size();i++){
        it=kDomSky.begin();
        while (it!=kDomSky.end()){
            k1=kCompare(donnees[i], *it, att);
            if (k1>=k && donnees[i][0]!=(*it)[0]) {it=kDomSky.erase(it);}
            else it++;
        }
    }
    result.swap(kDomSky);
}



DataType DPI_SkylineCubeSize(TableTuple &donnees, vector<Space> &att){
    DataType cubeSize=0;
    TableTuple Rk;
    Space k;
    for (k=1;k<=att.size();k++){
        DPI_kDom(donnees, att, k, Rk);
        cubeSize+=Rk.size();
        //cout<<Rk.size()<<" ";
    }
    //cout<<endl;
    return cubeSize;
}

#endif // KDOMSKYCUBE_H_INCLUDED
