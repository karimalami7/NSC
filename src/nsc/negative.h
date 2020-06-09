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

bool pairs_file_exists(string dataName, Space d, int n){

    string const nomFichier("../structPairs/"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    ifstream f(nomFichier.c_str());
    return f.good();

}

void load_pairs(string dataName, Space d, int n, vector<USetDualSpace> &listUSetDualSpace){

    cout <<"*****loading_pairs*****"<<endl<<endl;
    string const nomFichier("../structPairs/"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    ifstream file(nomFichier.c_str());
    string line;
    int line_counter=0;
    while (std::getline(file, line)){
        //cout << line<<endl;
        
        int pointer=0;
        while (pointer <line.size()){
            
            // First split by ";"
            string first_delimiter = ";";
            int first_delimiter_position= line.find(first_delimiter,pointer);
            //cout <<pointer<< " "<<first_delimiter_position<<endl;
            string token = line.substr(pointer, first_delimiter_position-pointer);
            
            // Second split by "|"
            {    
                string second_delimiter = "|";
                int second_delimiter_position= token.find(second_delimiter,0);
                string X=token.substr(0, second_delimiter_position);
                string Y=token.substr(second_delimiter_position+1, token.size());
                DualSpace ds;
                ds.dom=stoi(X);
                ds.equ=stoi(Y);
                listUSetDualSpace[line_counter].insert(ds);
                //cout << X << " " << Y<<endl;
            }

            pointer=first_delimiter_position+1;
            //cout <<token<<endl;
        }
        line_counter++;
    }
}

void print_pairs(string dataName, vector<USetDualSpace> &listUSetDualSpace, NegSkyStr &structureNSC, Space d){

    cout <<"*****printing_pairs*****"<<endl;

    int n=listUSetDualSpace.size();
    string const nomFichier1("../structPairs/"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    ofstream monFlux1(nomFichier1.c_str());
    if (monFlux1){
        for(int i=0; i<listUSetDualSpace.size(); i++){
            //monFlux1 <<"t"<<i<<": ";
            for (auto it_uset = listUSetDualSpace[i].begin(); it_uset!=listUSetDualSpace[i].end(); it_uset++){
                monFlux1 <<it_uset->dom <<"|"<<it_uset->equ <<";";
            }
            monFlux1 <<endl;
        }
        monFlux1.close();
    }
    else{
        cout << "ERROR: Couldn't open the file." << endl;
    }

    // string const nomFichier2("../structPairs/"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n)+"-T");
    // ofstream monFlux2(nomFichier2.c_str());
    // if (monFlux2){
    //     for(int i=0; i<structureNSC.size(); i++){
    //         monFlux2 <<structureNSC[i].first <<": ";
    //         for (auto Y = structureNSC[i].second.begin(); Y!=structureNSC[i].second.end(); Y++){
    //             monFlux2 << "( "<< Y->first << " | " ;
    //             for (auto id: Y->second){
    //                 monFlux2 << id << " ";
    //             }
    //             monFlux2 << ")"; 
    //         }
    //         monFlux2 <<endl;
    //     }
    //     monFlux2.close();
    // }
    // else{
    //     cout << "ERROR: Couldn't open the file." << endl;
    // }


    cout <<endl<<"*****done"<<endl;
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





long creationStructureNSC(NegSkyStrAux &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, vector<USetDualSpace> &listUSetDualSpace, Space d){
    Space all=(1<<d)-1;
    long structSize=0;
    DataType i;
    DataType nbTuples=0;

    for (i=0;i<(DataType)listUSetDualSpace.size();++i){ // on boucle sur tous les tuples taille n
        //if (listUSetDualSpace[i].size()!=1 || (listUSetDualSpace[i].begin())->dom<all){  // to discard completely dominated tuples 
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
        //}
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

    // Build a map that reduce the comparisons

    unordered_map<Space, TableTuple> maStructMap;

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




    vector<pair<Space, TableTuple>> maStructVec;
    vector<pair<Space, vector<DataType>>> maStructVec2;

    for (auto it1=maStructMap.begin();it1!=maStructMap.end();++it1){
        maStructVec.push_back(pair<Space, TableTuple>(it1->first, it1->second));
    }


    sort( maStructVec.begin(),  maStructVec.end(),pet);
    for(auto i=0; i< maStructVec.size();i++) triScore(maStructVec[i].second, d);

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
            size_t t1 = listUSetDualSpace[i].size();
            int nb=0;
            for(auto it=listUSetDualSpace[i].begin(); it!=listUSetDualSpace[i].end() ;++it) if(it->equ!=0) {nb++; break;}
            if((t1>1) && (nb!=0)){
                fusionGloutonne(listUSetDualSpace[i], d);
            }
        }
    }

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



bool* subspaceSkyline_NSC(NegSkyStr &structure, int data_size, const Space subspace){
    DataType i;
    bool* tableSkyline=new bool[data_size];
    for (i=0;i<data_size;i++) tableSkyline[i]=true;
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
    //for (int j=0; j<m;j++) if (tableSkyline[j]) cout << j << endl; // to print ids in skyline
    return tableSkyline;
}

DataType subspaceSkylineSize_NSC(NegSkyStr &structure, int data_size, Space subspace){
    
    DataType skySize;
    bool* skyline=subspaceSkyline_NSC(structure, data_size, subspace);
    skySize=compteLesMarques(skyline, data_size);
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

}
#endif // NEGATIVE_H_INCLUDED
