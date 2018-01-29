#ifndef EXPERIMENTATIONS_H_INCLUDED
#define EXPERIMENTATIONS_H_INCLUDED

#include "declarations.h"

void displayResult(string dataName, DataType n, Space d, DataType k, string step, long structSize, double timeToPerform, int method){
    cout<<dataName<<" "<<n<<" "<<d<<" "<<k<<" "<<methodNamesDisplayed[method]<<" "<<NB_THREADS<<" "<<step<<" "<<structSize<<" "<<timeToPerform<<endl;
}


void experimentation_NAIF(string dataName, TableTuple &donnees, Space d, DataType k, vector<vector<Space>> &listNTabSpace, vector<vector<Space>> &listAllTabSpace){
    double timeToPerform;
    long structSize=0;
    DataType i;
    DataType N=listNTabSpace.size();
    DataType All=listAllTabSpace.size();
    stringstream ss;ss<<"N="<<N;

    displayResult(dataName, donnees.size(), d, k, "BUILD", 0, 0, NAIF);

    structSize=0;
    timeToPerform=debut();
        for (i=0;i<N;i++){
            structSize+=subspaceSkylineSize_NAIF(donnees, listNTabSpace[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, ss.str(), structSize, timeToPerform, NAIF);

    structSize=0;
    timeToPerform=debut();
        for (i=0;i<All;i++){
            structSize+=subspaceSkylineSize_NAIF(donnees, listAllTabSpace[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "SKYCUBE", structSize, timeToPerform, NAIF);
}

void experimentation_TREE(string dataName, TableTuple &donnees, Space d, DataType k, vector<vector<Space>> &vectSpaceN, vector<vector<Space>> &vectSpaceAll){
    double timeToPerform;
    long structSize=0;
    DataType i;
    DataType N=vectSpaceN.size();
    DataType All=vectSpaceAll.size();
    stringstream ss;ss<<"N="<<N;

    displayResult(dataName, donnees.size(), d, k, "BUILD", 0, 0, TREE);

    structSize=0;
    timeToPerform=debut();
        for (i=0;i<N;i++){
            structSize+=subspaceSkylineSize_TREE(vectSpaceN[i], donnees);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, ss.str(), structSize, timeToPerform, TREE);

    structSize=0;
    timeToPerform=debut();
        for (i=0;i<All;i++){
            structSize+=subspaceSkylineSize_TREE(vectSpaceAll[i], donnees);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "SKYCUBE", structSize, timeToPerform, TREE);
}

void experimentation_CSC(string dataName, TableTuple &donnees, Space d, DataType k, vector<Space> &subspaceN, vector<Space> &subspaceAll){
    double timeToPerform;
    long structSize=0;
    DataType i;
    DataType N=subspaceN.size();
    DataType All=subspaceAll.size();
    stringstream ss;ss<<"N="<<N;
    DataType n=donnees.size();

    VectorUSetId structureCOMP(All+1);

    structSize=0;
    timeToPerform=debut();
        structSize+=compressedSkycube2(structureCOMP, donnees, d);
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "BUILD", structSize, timeToPerform, CSC);

    structSize=0;
    timeToPerform=debut();
        for (i=0;i<N;i++){
            structSize+=subspaceSkylineSize_CSC(donnees, structureCOMP, n, d, subspaceN[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, ss.str(), structSize, timeToPerform, CSC);

    structSize=0;
    timeToPerform=debut();
        for (i=0;i<All;i++){
            structSize+=subspaceSkylineSize_CSC(donnees, structureCOMP, n, d, subspaceAll[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "SKYCUBE", structSize, timeToPerform, CSC);
}

void experimentation_NSC(string dataName, TableTuple &donnees, Space d, DataType k, vector<Space> &subspaceN, vector<Space> &subspaceAll, bool nSC, bool kDom){

    double timeToPerform;
    double timeToPerform2;
    long structSize=0;
    long response_size=0;
    long structSize2;
    DataType i;
    DataType N=1;//subspaceN.size();
    DataType All=subspaceAll.size();
    DataType all=(1<<d)-1;
    stringstream ss;ss<<"N="<<N;
    //DataType n=donnees.size();

    NegSkyStr structureNSC;
    map<DataType, DataType> newIndexes;
    map<DataType, DataType> prvIndexes;

    structSize=0;
    timeToPerform=debut();
    DataType n=donnees.size();
    vector<USetDualSpace> listUSetDualSpace(n);
    TableTuple topmost;
    NEG::negativeSkycubeAux(listUSetDualSpace, donnees, topmost, d);
    structSize=NEG::negativeSkycube(structureNSC, newIndexes, prvIndexes, listUSetDualSpace, d);
    timeToPerform=duree(timeToPerform);
    structSize2=structSize;
    timeToPerform2=timeToPerform;
    displayResult(dataName, donnees.size(), d, k, "BUILD", structSize, timeToPerform, NSC);

    timeToPerform=debut();
    for (i=0;i<N;i++){
            
        response_size+=NEG::subspaceSkylineSize_NSC(structureNSC, newIndexes, prvIndexes, subspaceN[i]);

    }
    
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, ss.str(), response_size, timeToPerform, NSC);
    response_size=0;
    timeToPerform=debut();
        for (i=0;i<All;i++){
            response_size+=NEG::subspaceSkylineSize_NSC(structureNSC, newIndexes, prvIndexes, subspaceAll[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "SKYCUBE", response_size, timeToPerform, NSC);



    if (kDom){
        structSize=0;
        timeToPerform=debut();
            for (auto j=0;j<All;j++){
                for (i=1;i<=spaceSize(subspaceAll[j]);i++){
                    structSize+=NEG::subspaceSkylineSize_NSC_kDom(structureNSC, newIndexes, prvIndexes, subspaceAll[j], i);
                }
            }
        timeToPerform=duree(timeToPerform);
        if (nSC) cout<<endl;
        displayResult(dataName, donnees.size(), d, k, "KDOMSKYCUBE", structSize, timeToPerform, NSCk);
        displayResult(dataName, donnees.size(), d, k, "+BUILD", structSize2, timeToPerform+timeToPerform2, NSCk);
    }
    

}

void experimentation_NSCwM(string dataName, TableTuple &donnees, Space d, DataType k, vector<Space> &subspaceN, vector<Space> &subspaceAll, bool nSC, bool kDom){
    
    

    double timeToPerform;
    long structSize=0;
    long response_size=0;

    DataType i;
    DataType N=1;//subspaceN.size();
    DataType All=subspaceAll.size();
    DataType all=(1<<d)-1;
    stringstream ss;ss<<"N="<<N;
    

    NegSkyStr structureNSC;
    map<DataType, DataType> newIndexes;
    map<DataType, DataType> prvIndexes;

    
    timeToPerform=debut();
    DataType n=donnees.size();
    vector<mapDualSpace> listMapDualSpace(n); // Structure that stores pairs and counters: NSC
    TableTuple topmost;
    values_map topmost_map(d);

    //Build pairs
    NEG_wM::negativeSkycubeAux(listMapDualSpace, donnees, topmost, topmost_map, d);
    structSize=NEG_wM::negativeSkycube(structureNSC, newIndexes, prvIndexes, listMapDualSpace, d);
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "BUILD", structSize, timeToPerform, NSCwM);
    
    // identify IDs that are not in TM
    vector<int> notInTopmost(n), ids_donnees(n), ids_topmost(topmost.size());
    for (int i=0; i<n; i++) ids_donnees[i]=donnees[i][0];
    for (int i=0; i<topmost.size(); i++) ids_topmost[i]=topmost[i][0];
    sort(ids_topmost.begin(),ids_topmost.end());
    auto it1=std::set_difference (ids_donnees.begin(), ids_donnees.end(), ids_topmost.begin(), ids_topmost.end(),notInTopmost.begin());
    notInTopmost.resize(it1-notInTopmost.begin());
 

    timeToPerform=debut();
        for (i=0;i<N;i++){
            
            response_size+=NEG::subspaceSkylineSize_NSC(structureNSC, newIndexes, prvIndexes, subspaceN[i]);

        }
    
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, ss.str(), response_size, timeToPerform, NSCwM);
    response_size=0;
    timeToPerform=debut();
        for (i=0;i<All;i++){
            response_size+=NEG::subspaceSkylineSize_NSC(structureNSC, newIndexes, prvIndexes, subspaceAll[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "SKYCUBE", response_size, timeToPerform, NSCwM);

    //the application is interactive by the menu below 

    cerr<<mendl(3)<<"*********************NB*********************"<<mendl(2);
    cerr<<"*****Some options are not handled by this version and could lead to a misbehavior of the software. So, please: "<<endl;
    cerr<<"1/ Experiment Deletions and Insertions separately."<<endl;
    cerr<<"2/ Delete tuples in descending order of ids."<<endl;

    bool interactive_menu=true;

    while (interactive_menu){
        cerr << endl<< endl<< "****************choose an option****************"<< endl;

        cerr << "1 -> Deletion of one tuple"<< endl;

        cerr << "2 -> Insertion of one tuple" << endl;

        cerr << "3 -> Query on a subspace and skycube" << endl;

        cerr << "4 -> Deletion of multiple tuples" << endl;

        cerr << "5 -> Insertion of multiple tuples" << endl;

        cerr << "6 -> Tuple impact" << endl;

        cerr << "other -> out" << endl;

        int choix;
        cerr <<endl<<"Please, enter your choice: ";
        cin >> choix;
        cerr <<endl;

        switch(choix){

        case 1: cerr <<"++++++++ Deletion of one tuple"<<endl;
                NEG_wM::deletion_option(donnees, topmost, notInTopmost, listMapDualSpace, d, structureNSC, newIndexes, prvIndexes);
                NEG_wM::negativeSkycube(structureNSC, newIndexes, prvIndexes, listMapDualSpace, d);
                break;

        case 2: cerr <<"++++++++ Insertion of one tuple"<<endl;
                NEG_wM::insertion_option(dataName, k, donnees, topmost, listMapDualSpace, d, structureNSC, newIndexes, prvIndexes);
                NEG_wM::negativeSkycube(structureNSC, newIndexes, prvIndexes, listMapDualSpace, d);
                break;

        case 3: cerr <<"++++++++ Query on a subspace and skycube"<<endl;
                NEG::skylinequery(dataName, donnees, structureNSC, newIndexes, prvIndexes, d, k, subspaceN, subspaceAll);
                break;

        case 4: cerr <<"++++++++ Deletion of multiple tuples"<<endl;
                NEG_wM::multiple_deletion_option(dataName ,donnees, topmost, notInTopmost, listMapDualSpace, d, structureNSC, newIndexes, prvIndexes, topmost_map);
                NEG_wM::negativeSkycube(structureNSC, newIndexes, prvIndexes, listMapDualSpace, d);
                break;

        case 5: cerr <<"++++++++ Insertion of multiple tuples"<<endl;
                NEG_wM::multiple_insertion_option(dataName, k, donnees, topmost, listMapDualSpace, d, structureNSC, newIndexes, prvIndexes);
                NEG_wM::negativeSkycube(structureNSC, newIndexes, prvIndexes, listMapDualSpace, d);
                break;

        case 6: cerr <<"++++++++ Tuple impact"<<endl;
                NEG_wM::tuple_impact(dataName, donnees, topmost, listMapDualSpace, d);
                break;        

        default: exit(0);break;

        }

    }

    //cout<<endl;display_NSC(structureNSC, newIndexes, prvIndexes, d);cout<<endl;
}




void experimentation_BUA(string dataName, TableTuple &donnees, Space d, DataType k, vector<vector<Space>> &vectSpaceAll){
    double timeToPerform;
    long structSize=0;
    DataType All=vectSpaceAll.size();

    timeToPerform=debut();
        for (auto i=0;i<All;i++){
            structSize += BUA_SkylineCubeSize(donnees, vectSpaceAll[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "KDOMSKYCUBE", structSize, timeToPerform, BUA);
}

void experimentation_UBA(string dataName, TableTuple &donnees, Space d, DataType k, vector<vector<Space>> &vectSpaceAll){
    double timeToPerform;
    long structSize=0;
    DataType All=vectSpaceAll.size();

    timeToPerform=debut();
        for (auto i=0;i<All;i++){
            structSize += UBA_SkylineCubeSize(donnees, vectSpaceAll[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "KDOMSKYCUBE", structSize, timeToPerform, UBA);
}

void experimentation_DPI(string dataName, TableTuple &donnees, Space d, DataType k, vector<vector<Space>> &vectSpaceAll){
    double timeToPerform;
    long structSize=0;
    DataType All=vectSpaceAll.size();

    timeToPerform=debut();
        for (auto i=0;i<All;i++){
            structSize += DPI_SkylineCubeSize(donnees, vectSpaceAll[i]);
        }
    timeToPerform=duree(timeToPerform);
    displayResult(dataName, donnees.size(), d, k, "KDOMSKYCUBE", structSize, timeToPerform, DPI);
}

void experimentSkycube(string dataName, string path, DataType k, DataType n, Space d, const bool* selectedMethod){

    TableTuple donnees;
    //DataType nbSpace=1<<d;
    Space spAux;
    DataType i;

    Space N=1; //number of random queries
    Space All=(1<<d)-1; //number of queries in the skycube

    vector<Space> subspaceN;
    vector<vector<Space>> listNTabSpace(N);
    for (i=1;i<=N;i++){
        spAux=rand() % All + 1;
        subspaceN.push_back(spAux);
        listeAttributsPresents(spAux, d, listNTabSpace[i-1]);
        //displaySubspace(spAux, d);cout<<endl;
    }

    vector<Space> subspaceAll;
    vector<vector<Space>> listAllTabSpace(All);
    for (i=1;i<=All;i++){
        subspaceAll.push_back(i);
        listeAttributsPresents(i, d, listAllTabSpace[i-1]);
    }

    loadData(dataName, path, n, d, k, donnees);
    cerr << "Random query subspace: " << subspaceN[0]<<endl<<endl;


    if (selectedMethod[NAIF]) {experimentation_NAIF(dataName, donnees, d, k, listNTabSpace, listAllTabSpace);cout<<endl;}
    if (selectedMethod[TREE]) {experimentation_TREE(dataName, donnees, d, k, listNTabSpace, listAllTabSpace);cout<<endl;}
    if (selectedMethod[CSC]) {experimentation_CSC(dataName, donnees, d, k, subspaceN, subspaceAll);cout<<endl;}
    if (selectedMethod[NSC]) {experimentation_NSC(dataName, donnees, d, k, subspaceN, subspaceAll, true, selectedMethod[NSCk]);cout<<endl;}
    if (selectedMethod[NSCwM]) {experimentation_NSCwM(dataName, donnees, d, k, subspaceN, subspaceAll, true, selectedMethod[NSCk]);cout<<endl;}
    if (selectedMethod[BUA]) {experimentation_BUA(dataName, donnees, d, k, listAllTabSpace);cout<<endl;}
    if (selectedMethod[UBA]) {experimentation_UBA(dataName, donnees, d, k, listAllTabSpace);cout<<endl;}
    if (selectedMethod[DPI]) {experimentation_DPI(dataName, donnees, d, k, listAllTabSpace);cout<<endl;}
    if (selectedMethod[NSCk] && !selectedMethod[NSC]) {experimentation_NSC(dataName, donnees, d, k, subspaceN, subspaceAll, false, true);cout<<endl;}


    for (i=0;i<n;i++) delete[] donnees[i];
}

#endif // EXPERIMENTATIONS_H_INCLUDED
