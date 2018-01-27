#ifndef NEGATIVE_H_INCLUDED
#define NEGATIVE_H_INCLUDED

namespace NEG{

void triScore(TableTuple &donnees, Space d){
//  Cette procédure trie les données suivant le score (somme des valeurs des attributs)
    DataType i, n=donnees.size();
    Space j;
    vector<DataType> aux(n);
    TableTuple auxT=donnees;
    for (i=0;i<n;i++){
        aux[i]=donnees[i][1];
        for (j=2;j<=d;j++) aux[i]+=donnees[i][j];
    }
    vector<DataType> index(donnees.size());
    for (i = 0; i != (DataType)index.size(); ++i) index[i] = i;
    sortIndexes(aux, index);
    for (i=0;i<n;i++){
        donnees[i]=auxT[index[i]];
    }
}

void visualisation_pairs(vector<USetDualSpace> listUSetDualSpace){

    cout <<"*****visualisation_pairs*****"<<endl;
    for(int i=0; i<listUSetDualSpace.size(); i++){
        cout <<"t"<<i<<": ";
        for (auto it_uset = listUSetDualSpace[i].begin(); it_uset!=listUSetDualSpace[i].end(); it_uset++){
            cout <<it_uset->dom <<" "<<it_uset->equ <<" ; ";
        }
        cout <<endl;
    }

}

bool pet_pair(const DualSpace &sp1, const DualSpace &sp2 ){
    auto n1=sp1.dom + sp1.equ;
    auto n2=sp2.dom + sp2.equ;
    return __builtin_popcount(n1) >__builtin_popcount(n2);
}

inline DualSpace domDualSubspace_1(const Point &t1, const Point &t2, const Space &d){
//  retourne le sous espace dans lequel t1 domine t2//je pense que c'est t2 qui domine t1
//  les sous espaces sont codés en un nombre décimal
//  exemples (codage) quand d=4, ABCD->15, AD->9, A->1, B->2, C->4, BC -> 5
    Space j;
    Space poids1=0, poids2=0;
    DualSpace sortie;
    sortie.dom=0;
    sortie.equ=0;
    long pow=1;
    unsigned int dec=1;
    if (t1[0]==t2[0]) return sortie;
    for(j = 1; j <= d ; ++j){
        if(t1[j] < t2[j]) {sortie.dom+=pow;++poids1;}
        else if(t1[j] == t2[j]) {sortie.equ+=pow;++poids2;}
        pow*=2;
	   //pow = (pow << 1);
    }
    sortie.poids=(1<<(poids1+poids2))-(1<<poids2);
    return sortie;
}

void visualisation_total_pairs(TableTuple& donnees, Space d){

    cout <<"*****visualisation_total_pairs*****"<<endl;
    for(int i=0; i<donnees.size(); i++){
        cout <<"t"<<i<<": ";
        for(int j=0; j<donnees.size(); j++){
            DualSpace ds = domDualSubspace_1(donnees[j], donnees[i], d);
            cout << ds.dom<<" "<<ds.equ<<"; ";
        }
        cout <<endl;
    }

}

void insertDualSpaceToUSet(DualSpace sp, const Space &d, USetDualSpace &uSetDualSpace, Space &all, Space &maxDom, Space &sizeMaxDom, bool &sortie){
//retourne true si le tuple est completement dominé (ne plus le comparer à d'autres tuples)
    if (sp.dom==all){
        USetDualSpace sp1;
        sp1.insert(sp);
        uSetDualSpace.swap(sp1);
        sortie=true;
        return;
    }else if (sp.dom!=0){
        if (!estInclusDans(sp.dom+sp.equ, maxDom)){
            if (spaceSize(sp.dom)>sizeMaxDom){
                maxDom=sp.dom;
                sizeMaxDom=spaceSize(sp.dom);
            }
            uSetDualSpace.insert(sp);
        }
    }
}




// remplit une structure ou l element de base est des espaces 
long creationStructureNSC(NegSkyStrAux &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, vector<USetDualSpace> &listUSetDualSpace, Space d){
    Space all=(1<<d)-1;
    long structSize=0;
    DataType i;
    DataType nbTuples=0;

    for (i=0;i<(DataType)listUSetDualSpace.size();++i){ // on boucle sur tous les tuples taille n
        if (listUSetDualSpace[i].size()!=1 || (listUSetDualSpace[i].begin())->dom<all){  
            DataType idTuple=nbTuples;
            for (auto it=listUSetDualSpace[i].begin();it!=listUSetDualSpace[i].end();++it){ // on boucle sur tous les paris (X|Y) de ce tuple
                Space spaceXY=it->dom+it->equ;
                Space spaceY=it->equ;
                auto it2=structure.find(spaceXY);
                if (it2==structure.end()){
                    unordered_map<Space,vector<DataType>> mapAux;
                    vector<DataType> vectAux;
                    vectAux.push_back(idTuple);
                    mapAux.insert(pair<Space, vector<DataType>>(spaceY,vectAux));
                    structure.insert(pair<Space,unordered_map<Space,vector<DataType>>>(spaceXY,mapAux));
                }else{
                    auto it3=(it2->second).find(spaceY);
                    if (it3==(it2->second).end()){
                        vector<DataType> vectAux;
                        vectAux.push_back(idTuple);
                        (it2->second).insert(pair<Space, vector<DataType>>(spaceY,vectAux));
                    }else{
                        (it3->second).push_back(idTuple);
                    }
                }
                structSize++;
            }
            newIndexes[i]=nbTuples;
            prvIndexes[nbTuples]=i;
            nbTuples++;
        }
    }
    return structSize;
}



bool pet(const pair<Space, TableTuple> &p1, const pair<Space, TableTuple> &p2 ){
    return __builtin_popcount(p1.first) >__builtin_popcount(p2.first);
}

void choixPivot(TableTuple &topmost, Space d){
    DataType n= topmost.size();
    DataType iMinMax=-1, minMax=n, maxVal;
    for (auto i=0;i<n;i++){
        maxVal=topmost[i][0];
        for (auto j=1;j<d;j++) if (maxVal<topmost[i][j]) maxVal=topmost[i][j];
        if (maxVal<minMax){
            iMinMax=i;
            minMax=maxVal;
        }
    }
    Point pointAux=topmost[0];
    topmost[0]=topmost[iMinMax];
    topmost[iMinMax]=pointAux;
}


void negativequery(vector<USetDualSpace> matrix ){

    // to respond to negative query, i need the subspace and matrix
    Space querySubSpace = 14121;
    int responseSize=0;
    bool covered;
    for (int i = 0; i < matrix.size(); i++){
        covered=false;
        // check if querySubSpace is covered for each tuple
        for (auto it=matrix[i].begin();it!=matrix[i].end() && covered==false;++it){ // on boucle sur tous les paris (X|Y) de ce tuple
            Space spaceXY=it->dom+it->equ;

            if (estInclusDans(querySubSpace,spaceXY)){
                covered=true;
            }     
        }
        if (covered==false){
            responseSize++;
        }
    }
    cout<<"Response size: "<<responseSize<< endl;
}

void negativeSkycubeAux(vector<USetDualSpace> &listUSetDualSpace, TableTuple& donnees, TableTuple& topmost, Space d){
//  Cette procédure construit à partir des données la structure de compression de skycube suivant notre idéologie
    long structSize;
    DataType i;
    DataType n=donnees.size();
    Space j;
    vector<Space> attList;
    for (j=1;j<=d;j++) attList.push_back(j);
    ExecuteBSkyTree(attList, donnees, topmost);
    triScore(topmost, d);
    Space all=(1<<d)-1;


   // Space X,Y, delta;

    //Début construction de la Map
    unordered_map<Space, TableTuple> maStructMap;

    //choixPivot(topmost, d);
    Point pivot=topmost[0];
    DualSpace ds1;
    double debut=omp_get_wtime();
    for (i=1;i<topmost.size();i++){
        ds1=domDualSubspace_1(topmost[i], pivot, d);
        auto it1=maStructMap.find(ds1.dom);
        if (it1==maStructMap.end()){
            TableTuple tab(1);
            tab[0]=topmost[i];
            maStructMap[ds1.dom]=tab;
        }else{
            (it1->second).push_back(topmost[i]);
        }
    }
    //cerr<<"la construction de la map a pris "<<omp_get_wtime()-debut<< " sa taille est de "<<maStructMap.size()<<endl;
    //cerr<<"le topmost a une taille de "<<topmost.size()<<endl;
    //cerr<<"je suis ici \n";
    //Fin construction de la Map


    //Début construction du vecteur auxiliaire
    vector<pair<Space, TableTuple>> maStructVec;
    vector<pair<Space, vector<DataType>>> maStructVec2;

    for (auto it1=maStructMap.begin();it1!=maStructMap.end();++it1){
        maStructVec.push_back(pair<Space, TableTuple>(it1->first, it1->second));
    }

    //Fin construction du vecteur auxiliaire
    sort( maStructVec.begin(),  maStructVec.end(),pet);
    for(auto i=0; i< maStructVec.size();i++) triScore(maStructVec[i].second, d);

    //Space maxDom, sizeMaxDom;
    //bool sortie;
    #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
    for (i=0;i<n;++i){
        bool sortie;
        Space maxDom, sizeMaxDom;
         Space X,Y, delta;
        DualSpace ds;
        ds = domDualSubspace_1(pivot, donnees[i], d);
        listUSetDualSpace[i].insert(ds);
        sortie=false;
        maxDom = ds.dom;
        sizeMaxDom = spaceSize(ds.dom);

        if(ds.dom!=all){
            for (int l=0;l<maStructVec.size() && !sortie;l++){
                X=maStructVec[l].first;
                delta=X&ds.equ;
                //if (true){
                if (!estInclusDans(X, ds.dom+ds.equ)){
                    for (auto j=0; j<maStructVec[l].second.size() && !sortie;++j){
                        auto dds=domDualSubspace_1(maStructVec[l].second[j], donnees[i], d);
                        insertDualSpaceToUSet(dds, d, listUSetDualSpace[i], all, maxDom, sizeMaxDom, sortie);
                    }
                }else if (delta){
                    for (auto j=0; j<maStructVec[l].second.size() && !sortie;++j){
                        auto dds=domDualSubspace_1(maStructVec[l].second[j], donnees[i], d);
                            //DualSpace dds;dds.dom=delta;dds.equ=0;
                        insertDualSpaceToUSet(dds, d, listUSetDualSpace[i], all, maxDom, sizeMaxDom, sortie);
                    }
                }
            }
        }

        if(!sortie){
            //CompressionByInclusion
            bool trouve;
            auto it=listUSetDualSpace[i].begin();
            while (it!=listUSetDualSpace[i].end()){//pour chaque dual space, s'il est couvert alors on le supprime
                trouve=false;
                Space spdom=(*it).dom, spequ=(*it).equ;
                for(auto it1=listUSetDualSpace[i].begin();it1!=listUSetDualSpace[i].end()&& !trouve;it1++){
                    if(it!=it1){
                        if (estInclusDans(spdom,it1->dom) && (estInclusDans(spequ, it1->dom + it1->equ))){
                                trouve=true;//(*it) est entièrement couvert par (*it1)
                                it1=listUSetDualSpace[i].erase(it);
                                it=it1;
                        }
                    }
                }
                if(!trouve) it++;
            }
            //CompressionByGreedy
            // size_t t1 = listUSetDualSpace[i].size();
            // int nb=0;
            // for(auto it=listUSetDualSpace[i].begin(); it!=listUSetDualSpace[i].end() ;++it) if(it->equ!=0) {nb++; break;}
            // if((t1>1) && (nb!=0)){
            //     fusionGloutonne(listUSetDualSpace[i], d);
            // }
        }
    }

    // visualisation pairs

    //visualisation_pairs(listUSetDualSpace);


}




void display_NSC(NegSkyStr &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, Space d){
    for (auto itXY=structure.rbegin();itXY!=structure.rend();++itXY){
        displaySubspace(itXY->first, d);cout<<endl;
        for (auto itY=(itXY->second).begin();itY!=(itXY->second).end();++itY) {
            cout<<"   ";displaySubspace(itY->first, d);cout<<": ";
            for (auto itId=(itY->second).begin();itId!=(itY->second).end();++itId){
                cout<<prvIndexes[*itId]+1<<" ";
            }
            cout<<endl;
        }
    }
}

// down     ,  transform NegSkyStrAux to NegSkyStr,   map -> vector
long negativeSkycube(NegSkyStr &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, vector<USetDualSpace> listUSetDualSpace, Space d){

    long structSize;
    Space spXY, spY;

    NegSkyStrAux structure0;
    structSize = creationStructureNSC(structure0, newIndexes, prvIndexes, listUSetDualSpace, d);
    
    for (auto itXY=structure0.begin();itXY!=structure0.end();itXY++){
        spXY=itXY->first;
        vector<pair<Space,vector<DataType>>> vY;
        for (auto itY=(itXY->second).begin();itY!=(itXY->second).end();itY++){
            spY=itY->first;
            vector<Space> vId;
            for (auto itId=(itY->second).begin();itId!=(itY->second).end();itId++){
                vId.push_back(*itId);
            }
            vY.push_back(pair<Space, vector<DataType>>(spY, vId));
        }
        structure.push_back(pair<Space, vector<pair<Space,vector<DataType>>>>(spXY, vY));
    }
    
    //display_NSC(structure, newIndexes, prvIndexes, d);
    return structSize;
}



bool* subspaceSkyline_NSC(NegSkyStr &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, const Space subspace){
    DataType i;
    DataType m=newIndexes.size();
    bool* tableSkyline=new bool[m];
    for (i=0;i<m;i++) tableSkyline[i]=true;
    for (auto itXY=structure.rbegin();itXY!=structure.rend() && subspace<=(itXY->first);++itXY){
        if (estInclusDans(subspace, itXY->first)){
            for (auto itY=(itXY->second).begin();itY!=(itXY->second).end();++itY) {
                if (!estInclusDans(subspace, itY->first)) {
                    for (auto itId=(itY->second).begin();itId!=(itY->second).end();++itId){
                        tableSkyline[*itId]=false;
                    }
                }
            }
        }
    }
    // les vrais indices résultats de la requête sont les prvIndexes[i] tels que tableSkyline[i] vaut "true" et ce pour i allant de 0 à m
    return tableSkyline;
}

DataType subspaceSkylineSize_NSC(NegSkyStr &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, Space subspace){
    DataType skySize;

    bool* skyline=subspaceSkyline_NSC(structure, newIndexes, prvIndexes, subspace);
    skySize=compteLesMarques(skyline, newIndexes.size());
    delete[] skyline;
    return skySize;
}


bool* subspaceSkyline_NSC_kDom(NegSkyStr &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, const Space subspace, const Space k){
    DataType i;
    Space Z;
    DataType m=newIndexes.size();
    bool* tableSkyline=new bool[m];
    for (i=0;i<m;i++) tableSkyline[i]=true;
    for (auto itXY=structure.rbegin();itXY!=structure.rend();++itXY){
        Z = itXY->first & subspace;
        if (spaceSize(Z)>=k){
            for (auto itY=(itXY->second).begin();itY!=(itXY->second).end();++itY) {
                if (!estInclusDans(Z, itY->first)) {
                    for (auto itId=(itY->second).begin();itId!=(itY->second).end();++itId){
                        tableSkyline[*itId]=false;
                    }
                }
            }
        }
    }
    // les vrais indices résultats de la requête sont les prvIndexes[i] tels que tableSkyline[i] vaut "true" et ce pour i allant de 0 à m
    return tableSkyline;
}

DataType subspaceSkylineSize_NSC_kDom(NegSkyStr &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, Space subspace, Space k){
    DataType skySize;
    bool* skyline=subspaceSkyline_NSC_kDom(structure, newIndexes, prvIndexes, subspace, k);
    skySize=compteLesMarques(skyline, newIndexes.size());
    delete[] skyline;
    return skySize;
}



void displayResultv2(string dataName, DataType n, Space d, DataType k, string step, long structSize, double timeToPerform, int method){
    cout<<dataName<<" "<<n<<" "<<d<<" "<<k<<" "<<methodNamesDisplayed[method]<<" "<<NB_THREADS<<" "<<step<<" "<<structSize<<" "<<timeToPerform<<endl;
}

void skylinequery_bySH_option(string dataName, TableTuple &donnees, NegSkyStr structure0, map<DataType, DataType> newIndexes0, map<DataType, DataType> prvIndexes0, Space d, DataType k,vector<Space> &subspaceN, vector<Space> &subspaceAll){

    Space subspace=14121;

    //cerr << "Choose a subspace for the skyline query: ";

    //cin >> subspace;

    //query on 1

    int structSize=0;
    double timeToPerform=debut();
    structSize+=subspaceSkylineSize_NSC(structure0, newIndexes0, prvIndexes0, subspace);
    timeToPerform=duree(timeToPerform);
    displayResultv2(dataName, donnees.size(), d, k, "N=1", structSize, timeToPerform, NSC); 

    //query on all

    Space All=(1<<d)-1;
    structSize=0;
    timeToPerform=debut();
    for (int i=0;i<All;i++){
        structSize+=subspaceSkylineSize_NSC(structure0, newIndexes0, prvIndexes0, subspaceAll[i]);
    }
    timeToPerform=duree(timeToPerform);
    displayResultv2(dataName, donnees.size(), d, k, "SKYCUBE", structSize, timeToPerform, NSC);
}

}
#endif // NEGATIVE_H_INCLUDED
