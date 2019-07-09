#ifndef NEGATIVEWM_H_INCLUDED
#define NEGATIVEWM_H_INCLUDED
#include "sys/types.h"
#include "sys/sysinfo.h"
namespace NEG_wM {


int parseLine(char* line){
    // This assumes that a digit will be found, the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getMemoryValue(){ //this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}



void visualisation_pairs(vector<mapDualSpace> listMapDualSpace){

    cout <<"*****visualisation_pairs*****"<<endl;
    for(int i=0; i<listMapDualSpace.size(); i++){
        cout <<"t"<<i<<": ";
        for (auto it_uset = listMapDualSpace[i].begin(); it_uset!=listMapDualSpace[i].end(); it_uset++){
            cout <<it_uset->first.dom <<" "<<it_uset->first.equ <<" "<<it_uset->second <<" ; ";
        }
        cout <<endl;
    }

}



void build_pairs(mapDualSpace &mDS, TableTuple& donnees, TableTuple& topmost, Space d, int i){
    
    
    USetDualSpace usDS_local;
    mapDualSpace mDS_local;

    
    DualSpace ds;
    
    for (int l=0;l<topmost.size();l++){
        ds=NEG::domDualSubspace_1(topmost[l], donnees[i], d);

        auto it = mDS_local.find(ds);
        if (it != mDS_local.end()){
            it->second++;
        }
        else{
            mDS_local.insert(pair<DualSpace, int>(ds, 1));
        }
    }
    

    //  CompressionByInclusion

    list<DualSpace> maListe;
    for(auto it=mDS_local.begin(); it!=mDS_local.end(); ++it)maListe.push_back(it->first);
     
    maListe.sort(NEG::pet_pair);

   
    bool trouve;
    auto it=maListe.begin();
    while (it!=maListe.end()){//pour chaque dual space, s'il est couvert alors on le supprime
        trouve=false;
        Space spdom=(*it).dom, spequ=(*it).equ;
        auto it1=it;
        it1++;
        if(it1!=maListe.end()){
        while(it1!=maListe.end()){
            //if(it!=it1){
                
            trouve=false;
                if (estInclusDans(it1->dom,spdom) && (estInclusDans(it1->equ, spequ+ spdom))){
                             it1=maListe.erase(it1);
                            trouve=true;
                }
                if(!trouve)it1++;
            //}
        }
    }

        it++;
    }
    for(auto it=maListe.begin();it!=maListe.end();it++)usDS_local.insert(*it);


    // CompressionByGreedy
    // size_t t1 = usDS_local.size();
    // int nb=0;
    // for(auto it=usDS_local.begin(); it!=usDS_local.end() ;++it) if(it->equ!=0) {nb++; break;}
    // if((t1>1) && (nb!=0)){
    //     fusionGloutonne(usDS_local, d);
    // }
       
    mapDualSpace tmp_mDS;
    for (auto it_usDS_local=usDS_local.begin(); it_usDS_local != usDS_local.end(); it_usDS_local++ ){
        tmp_mDS[*it_usDS_local]=mDS_local[*it_usDS_local];
    }
    mDS_local.swap(tmp_mDS);

    mDS.swap(mDS_local);
}

void negativeSkycubeAux(vector<mapDualSpace> &listMapDualSpace, TableTuple& donnees, TableTuple& topmost,values_map& topmost_map, Space d){


    long structSize;
    DataType i;
    DataType n=donnees.size();
    Space j;
    vector<Space> attList;
    for (j=1;j<=d;j++) attList.push_back(j);
    ExecuteBSkyTree(attList, donnees, topmost);
    NEG::triScore(topmost, d);
    double debut=omp_get_wtime();


    #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
    for (i=0;i<n;++i){

        build_pairs(listMapDualSpace[i], donnees, topmost, d, i);   

    }

}



long creationStructureNSC(NegSkyStrAux &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, vector<mapDualSpace> &listMapDualSpace, Space d){
    Space all=(1<<d)-1;
    long structSize=0;
    DataType i;
    DataType nbTuples=0;

    for (i=0;i<(DataType)listMapDualSpace.size();++i){ // on boucle sur tous les tuples taille n
        if (listMapDualSpace[i].size()!=1 || (listMapDualSpace[i].begin())->first.dom<all){  
            DataType idTuple=nbTuples;
            for (auto it=listMapDualSpace[i].begin();it!=listMapDualSpace[i].end();++it){ // on boucle sur tous les paris (X|Y) de ce tuple
                Space spaceXY=it->first.dom+it->first.equ;
                Space spaceY=it->first.equ;
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

long negativeSkycube(NegSkyStr &structure, map<DataType,DataType> &newIndexes, map<DataType,DataType> &prvIndexes, vector<mapDualSpace> listMapDualSpace, Space d){

    long structSize;
    Space spXY, spY;

    NegSkyStrAux structure0;
    structSize = creationStructureNSC(structure0, newIndexes, prvIndexes, listMapDualSpace, d);
    
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

// Returns impact of every tuple in the topmost

void tuple_impact(string dataName, TableTuple& donnees, TableTuple& topmost, vector<mapDualSpace> &listMapDualSpace, Space d){

    double debut = omp_get_wtime();

    int n = donnees.size();

    cerr<< endl<<"This operation can take a while" <<mendl(2);
    cerr<< "Results will be printed in "<<"te-"<<dataName<<"-d-"<<std::to_string(d)<<"-n-"<<std::to_string(n)<<endl;
    string const nomFichier1("./te-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    ofstream monFlux1(nomFichier1.c_str());
    
    bool inclus;

    int total1;

    omp_lock_t writelock;

    omp_init_lock(&writelock);

    vector<int> total(topmost.size(),0);

    if(monFlux1){
        #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
        for (int i=0; i <topmost.size() ; i++){
            total1=0;

            
            for (int j=0; j <n; j++){
                if (i!=j){
                    DualSpace ds;
                    ds = NEG::domDualSubspace_1(topmost[i], donnees[j], d);

                    auto it = listMapDualSpace[j].find(ds);
                    if (it !=listMapDualSpace[j].end()){
                        if (it->second==1)
                        {
                            total1++; 
                        }
                    }
                }
            }

            total[i]=total1;
            
            
        }
        for (int i=0; i <topmost.size() ; i++){
            monFlux1 << "id: "<<topmost[i][0]<<", impact: "<< total[i]<<endl;
        }

        monFlux1.close();
    }
    else
    {
        cout << "ERROR: Couldn't open the file." << endl;
    }

    omp_destroy_lock(&writelock);

    cerr << "*done* "<< endl;
}




void identify_topmost_pairs(TableTuple& donnees, TableTuple& topmost,  TableTuple& tmp_topmost, int id_tuple_to_delete, vector<int> &notInTopmost, Space d){

    //1- identify descendants of t-
    TableTuple descendants;
    vector<int> descendants_id(notInTopmost.size(),0);

    Space all=(1<<d)-1;


    #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
    for (int i = 0; i < notInTopmost.size(); i++){
        DualSpace ds;
        ds = NEG::domDualSubspace_1(donnees[id_tuple_to_delete], donnees[notInTopmost[i]], d);
        if (ds.dom+ds.equ==all){
            descendants_id[i]=1;
        }
    }

    for (int i=0; i<descendants_id.size();i++){
        if (descendants_id[i]==1){
            descendants.push_back(donnees[notInTopmost[i]]);
        }
    }

    //2- identify the topmost of the descendants,

    TableTuple TM_descendants;
    vector<Space> attList;
    for (int l=1;l<=d;l++) attList.push_back(l);
    if (descendants.size()>0) ExecuteBSkyTree(attList, descendants, TM_descendants);
    
    //3- identify from descendants, those that are not dominated by topmost\t-

    vector<int> TM_descendants_id(TM_descendants.size(),0);
    #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
    for(int i = 0; i < TM_descendants.size(); i++){
        bool dominated = false;
        DualSpace ds;
        for (int j = 0; j< topmost.size(); j++ ){
            if (topmost[j][0]!=id_tuple_to_delete){
                ds = NEG::domDualSubspace_1(topmost[j], TM_descendants[i], d);
                if (ds.dom+ds.equ==all){
                    dominated=true;
                    break;
                }
            }
        }
        if (!dominated){

            TM_descendants_id[i]=1;

        }
    }
    for (int i =0; i<TM_descendants_id.size(); i++){
        if (TM_descendants_id[i]==1){
            tmp_topmost.push_back(TM_descendants[i]);        
        }
    
    }

}



void DeleteTuple(TableTuple& donnees, TableTuple& topmost, int id_tuple_to_delete ,vector<mapDualSpace> &listMapDualSpace, Space d, vector<int> notInTopmost){

    int n = donnees.size();


    TableTuple new_topmost;
    TableTuple desc_in_topmost;
    TableTuple tmp_donnees(donnees);

    bool inTopmost=false;

   //1-check if the tuple to delete is in topmost

    for (int i = 0; i < topmost.size(); i++){
        if (id_tuple_to_delete==topmost[i][0]){
            inTopmost=true;
            break;
        }
    }  

    //2- if yes, recalculate pairs of i if id_tuple_to_delete affects i

    if (inTopmost){        
        
        //2.1 identify descendants

        identify_topmost_pairs( donnees, topmost, desc_in_topmost,  id_tuple_to_delete, notInTopmost,  d);

        new_topmost=desc_in_topmost;
        for (int i = 0; i < topmost.size(); i++){
            
            if (topmost[i][0]!=id_tuple_to_delete){
             new_topmost.push_back(topmost[i]);
            }
        } 
        
        
        // 2.2 identify the impacted tuples
        std::vector<int> tuples_impacted(n);
        int total_tuples_impacted=0;
        #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic) //reduction (+:total_tuples_impacted)
        for (int i=0; i <n; i++){
        DualSpace ds;
            if (id_tuple_to_delete!=i){

                ds = NEG::domDualSubspace_1(donnees[id_tuple_to_delete], donnees[i], d);;

                auto it = listMapDualSpace[i].find(ds);
                if (it !=listMapDualSpace[i].end()){
                    if (it->second==1)
                    {   

                        listMapDualSpace[i].clear();
                        build_pairs(listMapDualSpace[i], donnees, new_topmost, d, i);
                    }else{
                        it->second--;
                    }

                } 
            }   
        }


        //updating data structures of NSCwM 

        topmost.swap(new_topmost);
        
        vector<int> vc_ids_desc(desc_in_topmost.size()),tmp_NTM(notInTopmost.size());
        for (int j=0; j < desc_in_topmost.size(); j++) vc_ids_desc[j]=desc_in_topmost[j][0];
        sort(desc_in_topmost.begin(),desc_in_topmost.end());
        auto it1=std::set_difference (notInTopmost.begin(), notInTopmost.end(), vc_ids_desc.begin(), vc_ids_desc.end(),tmp_NTM.begin());
        tmp_NTM.resize(it1-tmp_NTM.begin());
        notInTopmost.swap(tmp_NTM);
        
        donnees.erase(donnees.begin()+id_tuple_to_delete);
        listMapDualSpace.erase(listMapDualSpace.begin()+id_tuple_to_delete);

    }
    else{
        // if id_tuple_to_delete is not in TM
        donnees.erase(donnees.begin()+id_tuple_to_delete);

        listMapDualSpace.erase(listMapDualSpace.begin()+id_tuple_to_delete);


    }

}

void BatchDeleteSetOfTuples(TableTuple& donnees, TableTuple& topmost, vector<int> id_tuples_to_delete, vector<mapDualSpace> &listMapDualSpace, Space d){

    TableTuple tmp_topmost;
    TableTuple tmp_donnees(donnees);
    
    //1-calculate new topmost
    
    for (int i= id_tuples_to_delete.size()-1;i>=0; i-- ) tmp_donnees.erase(tmp_donnees.begin()+id_tuples_to_delete[i]);
     
    //cerr << "1"<<endl;    
    
    //2-detect impacted tuples
    int n = tmp_donnees.size();
    vector<int> vc_tuples_impacted(n,0);
    
    DualSpace ds;

    unordered_set<int> IdsDansTopmost;//contient les ids de l'ancien topmost
    for(int i=0; i<topmost.size(); i++ ){
        IdsDansTopmost.insert(topmost[i][0]);
    }

    vector <int> TopDel;//contient les supprimés qui sont dans topmost

    for(int i=0; i<id_tuples_to_delete.size(); i++ ){
        if(IdsDansTopmost.find(id_tuples_to_delete[i])!=IdsDansTopmost.end()) TopDel.push_back(id_tuples_to_delete[i]);
    }

    //cerr<<"la taille de topdel est "<<TopDel.size()<<endl;

    if(TopDel.size()>0){
        vector<Space> attList;
        vector<int> new_tuples_id;
        for (int l=1;l<=d;l++) attList.push_back(l);
        //on calcule le nouveau topmost    
        ExecuteBSkyTree(attList, tmp_donnees, tmp_topmost);


    
        for (int i=0; i<TopDel.size(); i++ ){
            #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
            for (int j=0 ; j < n; j++){
                ds = NEG::domDualSubspace_1(donnees[TopDel[i]], donnees[tmp_donnees[j][0]], d);
                //cout << ds.dom<<" " <<ds.equ<<endl;
                auto it = listMapDualSpace[tmp_donnees[j][0]].find(ds);
                if (it !=listMapDualSpace[tmp_donnees[j][0]].end()){
                    if (it->second==1)
                    {                        
                        vc_tuples_impacted[j]=1;                       
                       
                    }else{
                        it->second--;
                    }
                } 
            }
        }
    }

    //cerr << "2"<<endl;   
    //3- rebuild pairs of impaceted tuples 


    #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic) 
    for (int i =0; i < vc_tuples_impacted.size(); i++){
        if (vc_tuples_impacted[i]==1){

            listMapDualSpace[i].clear();
            build_pairs(listMapDualSpace[i], donnees, tmp_topmost, d, i);
        }
    }

    for (int i=id_tuples_to_delete.size()-1;i>=0; i-- ) {
        donnees.erase(donnees.begin()+id_tuples_to_delete[i]);
        listMapDualSpace.erase(listMapDualSpace.begin()+id_tuples_to_delete[i]);
    }
}











void InsertTuple(TableTuple &tuple, TableTuple &donnees, TableTuple &topmost, Space d, vector<mapDualSpace> &listMapDualSpace){

    int n=donnees.size();
    std::vector<int> notTopmostAnymore;
    Space all=(1<<d)-1;
    bool dominated = false;
    USetDualSpace usDS;
    mapDualSpace mDS;

    for (int i=0; i < topmost.size(); i++){
        DualSpace ds;
        ds=NEG::domDualSubspace_1(topmost[i], tuple[0], d);

        auto it = mDS.find(ds);
        if (it != mDS.end()){
            it->second++;
        }
        else{
            mDS.insert(pair<DualSpace, int>(ds, 1));
        } 
        if ((ds.dom+ds.equ)==all & ds.equ!=all){
            dominated=true;
        }
        else if(ds.dom==0 & ds.equ!=all) {
            notTopmostAnymore.push_back(i);
        }
    }
    // compression

    //  methode liste pour l'inclusion

    list<DualSpace> maListe;
    for(auto it=mDS.begin(); it!=mDS.end(); ++it)maListe.push_back(it->first);
     
    maListe.sort(NEG::pet_pair);
   
    bool trouve;
    auto it=maListe.begin();
    while (it!=maListe.end()){//pour chaque dual space, s'il est couvert alors on le supprime
        trouve=false;
        Space spdom=(*it).dom, spequ=(*it).equ;
        auto it1=it;
        it1++;
        if(it1!=maListe.end()){
        while(it1!=maListe.end()){
            //if(it!=it1){
                
            trouve=false;
                if (estInclusDans(it1->dom,spdom) && (estInclusDans(it1->equ, spequ+ spdom))){
                        //auto it_msDS = mDS_local.find(*it1);
                        //if(it_msDS!=mDS_local.end() && it_msDS->second>1 ){
                            //trouve=true;//(*it) est entièrement couvert par (*it1)
                             it1=maListe.erase(it1);
                            trouve=true;
                            //it=it1;
                        //}
                }
                if(!trouve)it1++;
            //}
        }
    }
        //if(!trouve) 
        it++;
    }

    //usDS_local.clear();
 
    for(auto it=maListe.begin();it!=maListe.end();it++)usDS.insert(*it);
    //  fin methode liste pour l'inclusion

    // size_t t1 = usDS.size();
    // int nb=0;
    // for(auto it=usDS.begin(); it!=usDS.end() ;++it) if(it->equ!=0) {nb++; break;}
    // if((t1>1) && (nb!=0)){
    //     fusionGloutonne(usDS, d);
    // }
       
    mapDualSpace tmp_mDS;
    for (auto it_usDS=usDS.begin(); it_usDS != usDS.end(); it_usDS++ ){
        tmp_mDS[*it_usDS]=mDS[*it_usDS];
    }
    mDS.swap(tmp_mDS);

    donnees.push_back(tuple[0]);
    listMapDualSpace.push_back(mDS);

    // cout << "usds"<< usDS.size()<<endl;
    // cout << "mds"<< mDS.size()<<endl;
    // cout <<"NIT: "<< notTopmostAnymore.size()<<endl;


    //if t+ not dominated

    if (!dominated){
        #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
        for (int i=0; i<donnees.size(); i++){

            bool new_p=false;
            DualSpace ds;
            ds=NEG::domDualSubspace_1(tuple[0], donnees[i], d);

            auto it_map = listMapDualSpace[i].find(ds);
            if (it_map != listMapDualSpace[i].end()){
                it_map->second++;
            }
            else{
                list<DualSpace> maListe;
                for(auto it=listMapDualSpace[i].begin(); it!=listMapDualSpace[i].end(); ++it)maListe.push_back(it->first);
                 
                maListe.sort(NEG::pet_pair);
                
                bool inclus=false;
                auto it=maListe.begin();
                while (it!=maListe.end()&&!inclus){//pour chaque dual space, s'il est couvert alors on le supprime
                    Space spdom=(*it).dom, spequ=(*it).equ;
                    if (estInclusDans(ds.dom,spdom) && (estInclusDans(ds.equ, spdom + spequ))){
                        inclus=true;
                    }
                    it++;
                }
                if (!inclus) {   
                    listMapDualSpace[i].insert(pair<DualSpace, int>(ds, 1));
                }    
            }
        }

        //update topmost

        for(int i=notTopmostAnymore.size()-1;i>=0;i--){
            topmost.erase(topmost.begin()+notTopmostAnymore[i]);
        }
        topmost.push_back(tuple[0]);

    }

}



void BatchInsertSetOfTuples(TableTuple &new_tuples, TableTuple &donnees, TableTuple &topmost, Space d, vector<mapDualSpace> &listMapDualSpace){

    int n=donnees.size();

    //calculate topmost of new tuples
    TableTuple topmost_new_tuples;
    Space j;
    vector<Space> attList;
    for (j=1;j<=d;j++) attList.push_back(j);
    ExecuteBSkyTree(attList, new_tuples, topmost_new_tuples);





    TableTuple union_TM = topmost ;
    for (int i =0;i < topmost_new_tuples.size();i++){
        union_TM.push_back(topmost_new_tuples[i]);
    }
    TableTuple union_TM_2;

    ExecuteBSkyTree(attList, union_TM, union_TM_2);



  
    // intersection entre TM_union et TM+

    std::vector<int> union_TM_2_id;
    for (int i =0;i < union_TM_2.size();i++) union_TM_2_id.push_back(union_TM_2[i][0]);
    sort(union_TM_2_id.begin(), union_TM_2_id.end());
    int pointeur;
    for ( pointeur =0;pointeur < union_TM_2_id.size();pointeur++) if (union_TM_2_id[pointeur]>=n) break;



    #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
    for (int i=0; i<donnees.size(); i++){
        USetDualSpace usDS;
        mapDualSpace mDS;
        bool new_p=false;
       
        for (int j=pointeur; j <union_TM_2_id.size();j++){
        //for (int j=0; j <inter_TMnew_union_TM.size();j++){
            DualSpace ds;
            ds=NEG::domDualSubspace_1(new_tuples[union_TM_2_id[j]-n], donnees[i], d);
            //ds=NEG::domDualSubspace_1(new_tuples[inter_TMnew_union_TM[j]-n], donnees[i], d);


            auto it = listMapDualSpace[i].find(ds);
            if (it != listMapDualSpace[i].end()){
                it->second++;
            }
            else{
                //new_p=true;
                //listMapDualSpace[i].insert(pair<DualSpace, int>(ds, 1));
                auto it1 = mDS.find(ds);
                if (it1!=mDS.end()){
                    it1->second++;
                }
                else{
                    mDS.insert(pair<DualSpace, int>(ds, 1));
                }
            }
        }

        //  methode liste pour l'inclusion
    
        list<DualSpace> maListe;
        for(auto it=listMapDualSpace[i].begin(); it!=listMapDualSpace[i].end(); ++it)maListe.push_back(it->first);
         
        maListe.sort(NEG::pet_pair);
    
        for ( auto it2 =mDS.begin(); it2!=mDS.end(); it2++){

            bool inclus=false;
            auto it=maListe.begin();
            while (it!=maListe.end()&&!inclus){//pour chaque dual space, s'il est couvert alors on le supprime
                Space spdom=(*it).dom, spequ=(*it).equ;
                if (estInclusDans(it2->first.dom,spdom) && (estInclusDans(it2->first.equ, spdom + spequ))){
                    inclus=true;
                }
                it++;
            }
            if (!inclus) {   
                listMapDualSpace[i].insert(pair<DualSpace, int>(it2->first, it2->second));
                maListe.push_back(it2->first);
            } 


        }

     

    //  fin methode liste pour l'inclusion



           

    }


    listMapDualSpace.resize(n+new_tuples.size());
    donnees.resize(n+new_tuples.size());
    #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
    for (int i=0; i < new_tuples.size();i++){
        USetDualSpace usDS;
        for (int j=0; j <union_TM_2.size();j++){
            DualSpace ds;
            ds=NEG::domDualSubspace_1(union_TM_2[j], new_tuples[i], d);

            auto it = listMapDualSpace[i+n].find(ds);
            if (it != listMapDualSpace[i+n].end()){
                it->second++;
            }
            else{
                listMapDualSpace[i+n].insert(pair<DualSpace, int>(ds, 1));
            }
        }


        //  CompressByInclusion
    
        list<DualSpace> maListe;
        for(auto it=listMapDualSpace[i+n].begin(); it!=listMapDualSpace[i+n].end(); ++it)maListe.push_back(it->first);
         
        maListe.sort(NEG::pet_pair);
    
       
        bool trouve;
        auto it=maListe.begin();
        while (it!=maListe.end()){//
            trouve=false;
            Space spdom=(*it).dom, spequ=(*it).equ;
            auto it1=it;
            it1++;
            if(it1!=maListe.end()){
            while(it1!=maListe.end()){
                //if(it!=it1){
                    
                trouve=false;
                    if (estInclusDans(it1->dom,spdom) && (estInclusDans(it1->equ, spequ+ spdom))){
                            //auto it_msDS = mDS_local.find(*it1);
                            //if(it_msDS!=mDS_local.end() && it_msDS->second>1 ){
                                //trouve=true;//(*it) est entièrement couvert par (*it1)
                                 it1=maListe.erase(it1);
                                trouve=true;
                                //it=it1;
                            //}
                    }
                    if(!trouve)it1++;
                //}
            }
        }
            //if(!trouve) 
            it++;
        }
    
        //usDS_local.clear();
     
        for(auto it=maListe.begin();it!=maListe.end();it++)usDS.insert(*it);

        // CompressionByGreedy
        // size_t t1 = usDS.size();
        // int nb=0;
        // for(auto it=usDS.begin(); it!=usDS.end() ;++it) if(it->equ!=0) {nb++; break;}
        // if((t1>1) && (nb!=0)){
        //     fusionGloutonne(usDS, d);
        // }
           
        mapDualSpace tmp_mDS;
        for (auto it_usDS_local=usDS.begin(); it_usDS_local != usDS.end(); it_usDS_local++ ){
            tmp_mDS[*it_usDS_local]=listMapDualSpace[i+n][*it_usDS_local];
        }
        listMapDualSpace[i+n].swap(tmp_mDS);
        donnees[i+n]=new_tuples[i];
    }

}












void multiple_deletion_option(string dataName, TableTuple& donnees, TableTuple& topmost, std::vector<int> notInTopmost, vector<mapDualSpace> &listMapDualSpace,Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0, values_map& topmost_map){

    int n =donnees.size();

    int size;

    vector<int> ids_to_delete;

    cerr <<"Please, enter the size of the subset to delete: ";

    cin >> size;

    structure0.clear();

    newIndexes0.clear();

    prvIndexes0.clear();

    ids_to_delete.clear();
    for ( int j = n; j>n-size; j-- ) {
        ids_to_delete.push_back(j-1);
    }
    sort(ids_to_delete.begin(),ids_to_delete.end());

    double timeToPerform=debut();
    BatchDeleteSetOfTuples(donnees, topmost, ids_to_delete, listMapDualSpace, d);
    timeToPerform=duree(timeToPerform);

    cerr <<endl<<"BatchDeleteSetOfTuples, size: " << size << ", time: "<< timeToPerform<<std::endl;

}

void deletion_option(TableTuple& donnees, TableTuple& topmost, std::vector<int> notInTopmost, vector<mapDualSpace> &listMapDualSpace, Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0){

    int id_tuple_to_delete;

    cerr << "Please, enter the id of the tuple to delete: ";

    cin >> id_tuple_to_delete;

    structure0.clear();

    newIndexes0.clear();

    prvIndexes0.clear();

    double timeToPerform=debut();

    DeleteTuple(donnees, topmost, id_tuple_to_delete, listMapDualSpace, d, notInTopmost);

    timeToPerform=duree(timeToPerform);

    std::cerr << "DeleteTuple, time: "<< timeToPerform<<std::endl; 

}







void multiple_insertion_option(string dataName,  DataType k, TableTuple& donnees, TableTuple& topmost, std::vector<mapDualSpace> &listMapDualSpace, Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0){  

    int n=donnees.size();

    int size;

    TableTuple new_tuples;

    cerr <<"Please, enter the size of the subset to insert: ";

    cin >> size;
    
    // generate new tuples
    loadData(dataName, "", size, d, k, new_tuples);

    for (int k=0 ; k< new_tuples.size();k++){
        new_tuples[k][0]=n+k;
    }

    structure0.clear();

    newIndexes0.clear();

    prvIndexes0.clear();

    double timeToPerform=debut();

    BatchInsertSetOfTuples(new_tuples, donnees, topmost, d, listMapDualSpace);

    timeToPerform=duree(timeToPerform);

    cerr << "BatchInsertSetOfTuples, size: "<< size<<", time: "<< timeToPerform<<std::endl;    

}

void insertion_option(string dataName,  DataType k, TableTuple& donnees, TableTuple& topmost, vector<mapDualSpace> &listMapDualSpace,Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0){

    int n = donnees.size();

    TableTuple new_tuples;

    loadData(dataName, "", 1, d, k, new_tuples);

    for (int i=0; i<new_tuples.size(); i++){

    TableTuple new_tuple;

    new_tuple.push_back(new_tuples[i]);

    new_tuple[0][0]=donnees.size();

    structure0.clear();

    newIndexes0.clear();

    prvIndexes0.clear();

    double timeToPerform=debut();

    InsertTuple(new_tuple, donnees, topmost, d, listMapDualSpace);

    timeToPerform=duree(timeToPerform);

    cerr << "InsertTuple, time: "<< timeToPerform<<std::endl;

    }
}


}
#endif // NEGATIVEWM_H_INCLUDED