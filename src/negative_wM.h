#ifndef NEGATIVEWM_H_INCLUDED
#define NEGATIVEWM_H_INCLUDED
#include "sys/types.h"
#include "sys/sysinfo.h"
namespace NEG_wM {


int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getMemoryValue(){ //Note: this value is in KB!
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



void build_pairs(USetDualSpace &usDS ,mapDualSpace &mDS, TableTuple& donnees, TableTuple& topmost, Space d, int i){
    
    
    USetDualSpace usDS_local;
    mapDualSpace mDS_local;

    
    DualSpace ds;
    
    for (int l=0;l<topmost.size();l++){
        ds=NEG::domDualSubspace_1(topmost[l], donnees[i], d);

       // usDS_local.insert(ds);
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

    usDS.swap(usDS_local);
    mDS.swap(mDS_local);
}

void negativeSkycubeAux(vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace, TableTuple& donnees, TableTuple& topmost,values_map& topmost_map, Space d){


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

        USetDualSpace USDS;

        build_pairs(listUSetDualSpace[i], listMapDualSpace[i], donnees, topmost, d, i);
      

        //listUSetDualSpace[i].clear();
        //cout << "uset  " <<  listUSetDualSpace[i].size() << "  map  " << listMapDualSpace[i].size() << endl;
        //if (((i+1) % 1000) ==0)cout << i<<" "<<getMemoryValue() << endl;
    }

    // visualisation pairs

    //visualisation_pairs(listUSetDualSpace);

}
















// tuple_effect: calcul l'impact d'un tuple sur T

void tuple_effect(string dataName, TableTuple& donnees, TableTuple& topmost, vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace, Space d){

    double debut = omp_get_wtime();

    int n = donnees.size();

    

    string const nomFichier1("./experiments/tuple_effect/testmap-te-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
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
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }

    omp_destroy_lock(&writelock);

    cout << "time_tuple_effect: " << omp_get_wtime() - debut << endl;
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



void DeleteTuple(TableTuple& donnees, TableTuple& topmost, int id_tuple_to_delete ,vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace, Space d, vector<int> notInTopmost){

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

        // tmp_donnees.erase(tmp_donnees.begin()+id_tuple_to_delete);
        // vector<Space> attList;
        // for (int l=1;l<=d;l++) attList.push_back(l);
        // ExecuteBSkyTree(attList, tmp_donnees, new_topmost);    
        //cout << "TM skytree"<< new_topmost.size()<<endl;    
        
        //2.1 identify descendants

        identify_topmost_pairs( donnees, topmost, desc_in_topmost,  id_tuple_to_delete, notInTopmost,  d);
        new_topmost=desc_in_topmost;
        for (int i = 0; i < topmost.size(); i++){
            
            if (topmost[i][0]!=id_tuple_to_delete){
             new_topmost.push_back(topmost[i]);
            }
        } 
        //cerr<<"-- desc size: "<<desc_in_topmost.size()<<endl;
        
        
        // 2.2 identify the impacted tuples
        std::vector<int> tuples_impacted(n);
        int total_tuples_impacted=0;
        #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic) //reduction (+:total_tuples_impacted)
        for (int i=0; i <n; i++){
        DualSpace ds;
            if (id_tuple_to_delete!=i){

                ds = NEG::domDualSubspace_1(donnees[id_tuple_to_delete], donnees[i], d);
                //cout << ds.dom<<" " <<ds.equ<<endl;

                auto it = listMapDualSpace[i].find(ds);
                if (it !=listMapDualSpace[i].end()){
                    if (it->second==1)
                    {   
                        //tuples_impacted[i]=1;
                        //total_tuples_impacted++;
                        USetDualSpace USDS;
                        //listUSetDualSpace[i].clear();
                        listMapDualSpace[i].clear();
                        build_pairs(USDS, listMapDualSpace[i], donnees, new_topmost, d, i);
                        //listMapDualSpace[i].erase(it);
                    }else{
                        it->second--;
                    }

                } 
            }   
        }





        // updating data structures of NSCwM 

        // topmost.swap(new_topmost);
        
        // vector<int> vc_ids_desc(desc_in_topmost.size()),tmp_NTM(notInTopmost.size());
        // for (int j=0; j < desc_in_topmost.size(); j++) vc_ids_desc[j]=desc_in_topmost[j][0];
        // sort(desc_in_topmost.begin(),desc_in_topmost.end());
        // auto it1=std::set_difference (notInTopmost.begin(), notInTopmost.end(), vc_ids_desc.begin(), vc_ids_desc.end(),tmp_NTM.begin());
        // tmp_NTM.resize(it1-tmp_NTM.begin());
        // notInTopmost.swap(tmp_NTM);
        
        // donnees.erase(donnees.begin()+id_tuple_to_delete);
        // listMapDualSpace.erase(listMapDualSpace.begin()+id_tuple_to_delete);
        // listUSetDualSpace.erase(listUSetDualSpace.begin()+id_tuple_to_delete);
    }
    else{
        // if id_tuple_to_delete is not in TM
        donnees.erase(donnees.begin()+id_tuple_to_delete);
        listMapDualSpace.erase(listMapDualSpace.begin()+id_tuple_to_delete);
        listUSetDualSpace.erase(listUSetDualSpace.begin()+id_tuple_to_delete);
        //notInTopmost.erase();

    }

}

void BatchDeleteSetOfTuples(TableTuple& donnees, TableTuple& topmost, vector<int> id_tuples_to_delete ,vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace, Space d){

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

            USetDualSpace USDS;
            listMapDualSpace[i].clear();
            build_pairs(USDS, listMapDualSpace[i], donnees, tmp_topmost, d, i);
        }
    }
}











void InsertTuple(TableTuple &tuple, TableTuple &donnees, TableTuple &topmost, Space d, vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace){

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

    donnees.push_back(tuple[0]);
    listUSetDualSpace.resize(n+1);
    //listUSetDualSpace[n]=usDS;
    listMapDualSpace.push_back(mDS);

}



void BatchInsertSetOfTuples(TableTuple &new_tuples, TableTuple &donnees, TableTuple &topmost, Space d, vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace){

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

   /* std::vector<int> topmost_new_tuples_id;
    for (int i =0;i < topmost_new_tuples.size();i++) topmost_new_tuples_id.push_back(topmost_new_tuples[i][0]);
    sort(topmost_new_tuples_id.begin(), topmost_new_tuples_id.end());
    vector<int> inter_TMnew_union_TM(topmost_new_tuples_id.size());
    auto it = set_intersection(union_TM_2_id.begin(),union_TM_2_id.end(),topmost_new_tuples_id.begin(), topmost_new_tuples_id.end(),inter_TMnew_union_TM.begin());
    inter_TMnew_union_TM.resize(it-inter_TMnew_union_TM.begin());
    cout << "inter TM+ union TM: "<< inter_TMnew_union_TM.size()<<endl;
*/

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



    listUSetDualSpace.resize(n+new_tuples.size());
    listMapDualSpace.resize(n+new_tuples.size());
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


        //  methode liste pour l'inclusion
    
        list<DualSpace> maListe;
        for(auto it=listMapDualSpace[i+n].begin(); it!=listMapDualSpace[i+n].end(); ++it)maListe.push_back(it->first);
         
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
        for (auto it_usDS_local=usDS.begin(); it_usDS_local != usDS.end(); it_usDS_local++ ){
            tmp_mDS[*it_usDS_local]=listMapDualSpace[i+n][*it_usDS_local];
        }
        listMapDualSpace[i+n].swap(tmp_mDS);
        //listUSetDualSpace[i+n]=usDS;

    }


}






















void multiple_deletion_option(string dataName, TableTuple& donnees, TableTuple& topmost, std::vector<int> notInTopmost, vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace,Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0, values_map& topmost_map){

    int n =donnees.size();

    int id_to_delete;

    vector<int> ids_to_delete;

    // generate ids_to_delete with same proportion of topmost and notTopmost

    // double proportionTM = float(topmost.size())/n;

    // for (int i =0; i< int (n*(0.01*proportionTM)); i++){
    //     //generate ids
    //     id_to_delete=rand()%topmost.size();

    //     ids_to_delete.push_back(topmost[id_to_delete][0]);
               
    // }

    // for (int i =0; i < ((n/100)- (n*(0.01*proportionTM))) ; i++){

    //     id_to_delete=rand()%notInTopmost.size();

    //     ids_to_delete.push_back(notInTopmost[id_to_delete]);

    // }

   // deletion NSCwM x1

    cerr <<"**** Begin: Deletion NSCwM x1"<<endl;

    string const nomFichier1("./experiments/deletion_NSCwM/del-x1-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    ofstream monFlux1(nomFichier1.c_str());
    if (!monFlux1){
        cerr << "Impossible d'ouvrir le fichier" << nomFichier1 << endl;
    }

    for(int i = 0; i < topmost.size(); i++){

        ids_to_delete.push_back(topmost[i][0]);
    
    }

    for (int i=0; i< ids_to_delete.size(); i++){

        TableTuple copy_donnees(donnees);

        TableTuple copy_topmost(topmost);

        std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

        std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

        std::vector<int> copy_notInTopmost(notInTopmost);

        double timeToPerform=debut();
    
        DeleteTuple( copy_donnees, copy_topmost, ids_to_delete[i], copy_listUSetDualSpace, copy_listMapDualSpace, d, copy_notInTopmost);
    
        timeToPerform=duree(timeToPerform);

        monFlux1 <<"id: " << ids_to_delete[i] << ", time: "<< timeToPerform<<std::endl;       

    }



    // Comparaison results between add_pairs and rebuild_pairs

    // cerr <<"**** Begin: Comparaison between add_pairs and rebuild_pairs, by Skycube result"<<endl;

    // Space N=1; //number of random queries
    // Space All=(1<<d)-1; //number of queries in the skycube
    // Space spAux;
    // vector<Space> subspaceN;
    // vector<vector<Space>> listNTabSpace(N);
    // for (int i=1;i<=N;i++){
    //     spAux=rand() % All + 1;
    //     subspaceN.push_back(spAux);

    //     //displaySubspace(spAux, d);cout<<endl;
    // }

    // vector<Space> subspaceAll;
    // vector<vector<Space>> listAllTabSpace(All);
    // for (int i=1;i<=All;i++){
    //     subspaceAll.push_back(i);

    // }

    // // run delete_tuple with copies of data
    // cerr << "query with add_pairs"<< endl;
    // for (int i=5859; i< 5860; i++){

    //     TableTuple copy_donnees(donnees);

    //     TableTuple copy_topmost(topmost);

    //     std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

    //     std::vector<int> copy_notInTopmost(notInTopmost);

    //     double timeToPerform=debut();
    
    //     delete_tuple( copy_donnees, copy_topmost, i, copy_listUSetDualSpace, copy_listMapDualSpace, d, copy_notInTopmost, topmost_map);
    
    //     timeToPerform=duree(timeToPerform);

    //     structure0.clear();

    //     newIndexes0.clear();

    //     prvIndexes0.clear();

    //     NEG::negativeSkycube(structure0, newIndexes0, prvIndexes0, copy_listUSetDualSpace, d);

    //     NEG::skylinequery_bySH_option("INDE", copy_donnees, structure0, newIndexes0, prvIndexes0, d, 100, subspaceN, subspaceAll);        

    // }
    // cerr << "query with rebuild_pairs"<< endl;
    // for (int i=5859; i< 5860; i++){

    //     TableTuple copy_donnees(donnees);

    //     TableTuple copy_topmost(topmost);

    //     std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

    //     std::vector<int> copy_notInTopmost(notInTopmost);

    //     double timeToPerform=debut();
    
    //     DeleteTuple( copy_donnees, copy_topmost, i, copy_listUSetDualSpace, copy_listMapDualSpace, d, copy_notInTopmost);
    
    //     timeToPerform=duree(timeToPerform);

    //     structure0.clear();

    //     newIndexes0.clear();

    //     prvIndexes0.clear();

    //     NEG::negativeSkycube(structure0, newIndexes0, prvIndexes0, copy_listUSetDualSpace, d);

    //     NEG::skylinequery_bySH_option("INDE", copy_donnees, structure0, newIndexes0, prvIndexes0, d, 100, subspaceN, subspaceAll);        

    // }


    //deletion NSCwM x21

    // cerr <<"**** Begin: Deletion NSCwM x21"<<endl;

    // string const nomFichier1("./experiments/deletion_NSCwM/del-x21-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    // ofstream monFlux1(nomFichier1.c_str());
    // if (!monFlux1){
    //     cerr << "Impossible d'ouvrir le fichier" << nomFichier1 << endl;
    // }

    
    // int fin=90;
    // int pas=20;
    // int deb=90;

    // for (int i=deb ; i<=fin; i=i+pas){

    //     cout << "-- sup: "<< i<<endl;

    //     //   sequential deletion of a subset x21

    //     // cout << "/seq"<<endl;

    //     TableTuple copy_donnees(donnees);

    //     TableTuple copy_topmost(topmost);

    //     std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

    //     std::vector<int> copy_notInTopmost(notInTopmost);

    //     double timeToPerform=debut();

    //     // for ( int j = n-1; j>n-i-1; j-- ) {

    //     //     DeleteTuple( copy_donnees, copy_topmost, j, copy_listUSetDualSpace, copy_listMapDualSpace, d, copy_notInTopmost);
            
    //     // }

    //     // timeToPerform=duree(timeToPerform);

    //     // monFlux1 <<"s-j: " << i << ", time: "<< timeToPerform<<std::endl;

    //     //grouped deletion of a subset x21

    //     cout <<"/batch"<<endl;

    //     ids_to_delete.clear();

    //     for ( int j = n; j>n-i; j-- ) {
    //         ids_to_delete.push_back(j-1);
    //     }

    //     sort(ids_to_delete.begin(),ids_to_delete.end());

    //     copy_donnees=donnees;

    //     copy_topmost=topmost;

    //     copy_listUSetDualSpace=listUSetDualSpace;

    //     copy_listMapDualSpace=listMapDualSpace;

    //     timeToPerform=debut();

    //     BatchDeleteSetOfTuples(copy_donnees, copy_topmost, ids_to_delete ,copy_listUSetDualSpace, copy_listMapDualSpace, d);

    //     timeToPerform=duree(timeToPerform);

    //     monFlux1 <<"b-j: " << i << ", time: "<< timeToPerform<<std::endl;

    // }

        //deletion NSCwM x22

    // cerr <<"**** Begin: Deletion NSCwM x22"<<endl;

    // string const nomFichier1("./experiments/deletion_NSCwM/del-x22-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    // ofstream monFlux1(nomFichier1.c_str());
    // if (!monFlux1){
    //     cerr << "Impossible d'ouvrir le fichier" << nomFichier1 << endl;
    // }

    
    // int fin=100000;
    // int pas=10;
    // int deb=1000;

    // for (int i=deb ; i<=fin; i=i+pas){

    //     // cout << "-- sup: "<< i<<endl;

    //     // //   sequential deletion of a subset x21

    //     // cout << "/seq"<<endl;

    //     // TableTuple copy_donnees(donnees);

    //     // TableTuple copy_topmost(topmost);

    //     // std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     // std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

    //     // std::vector<int> copy_notInTopmost(notInTopmost);

    //     // double timeToPerform=debut();

    //     // for ( int j = n-1; j>n-i-1; j-- ) {

    //     //     DeleteTuple( copy_donnees, copy_topmost, j, copy_listUSetDualSpace, copy_listMapDualSpace, d, copy_notInTopmost);
            
    //     // }

    //     // timeToPerform=duree(timeToPerform);

    //     // monFlux1 <<"s-j: " << i << ", time: "<< timeToPerform<<std::endl;

    //     //grouped deletion of a subset x22

    //     cout <<"/batch"<<endl;

    //     ids_to_delete.clear();

    //     for ( int j = n; j>n-i; j-- ) {
    //         ids_to_delete.push_back(j-1);
    //     }

    //     sort(ids_to_delete.begin(),ids_to_delete.end());

    //     TableTuple copy_donnees(donnees);

    //     TableTuple copy_topmost(topmost);

    //     std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

    //     double timeToPerform=debut();

    //     BatchDeleteSetOfTuples(copy_donnees, copy_topmost, ids_to_delete ,copy_listUSetDualSpace, copy_listMapDualSpace, d);

    //     timeToPerform=duree(timeToPerform);

    //     monFlux1 <<"b-j: " << i << ", time: "<< timeToPerform<<std::endl;

    // }
}

void deletion_option(TableTuple& donnees, TableTuple& topmost, std::vector<int> notInTopmost, vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace, Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0){

    int id_tuple_to_delete;

   // cerr << "enter the id of the tuple to delete: ";

    //cin >> id_tuple_to_delete;

    structure0.clear();

    newIndexes0.clear();

    prvIndexes0.clear();

    double timeToPerform=debut();
    
    //tuple_effect(donnees, id_tuple_to_delete, listUSetDualSpace, d);

    //delete_tuple(donnees, topmost, 32, listUSetDualSpace, listMapDualSpace, d, notInTopmost);

    timeToPerform=duree(timeToPerform);

    std::cerr << "delete tuple, time: "<< timeToPerform<<std::endl; 

    timeToPerform=debut();

    NEG::negativeSkycube(structure0, newIndexes0, prvIndexes0, listUSetDualSpace, d);

    timeToPerform=duree(timeToPerform);

    std::cerr << "fin SH, time: "<< timeToPerform<<std::endl;

}







void multiple_insertion_option(string dataName,  DataType k, TableTuple& donnees, TableTuple& topmost, vector<USetDualSpace> &listUSetDualSpace, std::vector<mapDualSpace> &listMapDualSpace, Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0){  

    int n=donnees.size();
    
    cerr <<"**** Begin: Multiple Insertion NSCwM"<<endl;

    // // x11

    // string const nomFichier1("./experiments/insertion_NSCwM/ins-x11-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    // ofstream monFlux1(nomFichier1.c_str());
    // if (!monFlux1){
    //     cerr << "Impossible d'ouvrir le fichier" << nomFichier1 << endl;
    // }

    // cerr <<"-x11"<<endl;

    // int fin=110;

    // int pas=20;

    // int deb=10;
    
    // for (int i=deb;i<=fin;i=i+pas){
        
    //     cerr <<"= insertion: "<<i<<endl;

    //     TableTuple new_tuples;
        
    //     // generate new tuples
    //     loadData(dataName, "", i, d, k, new_tuples);

    //     for (int k=0 ; k< new_tuples.size();k++){
    //         new_tuples[k][0]=n+k;
    //     }

    //     // sequential

    //     cerr <<"/sequential"<<endl;

    //     TableTuple copy_donnees(donnees);

    //     TableTuple copy_topmost(topmost);

    //     std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

    //     double timeToPerform=debut();

    //     for (int j = 0; j< i;j++){

    //         TableTuple new_tuple;
    //         new_tuple.push_back(new_tuples[j]);
    //         new_tuple[0][0]=n+j;
    //         InsertTuple(new_tuple, copy_donnees, copy_topmost, d, copy_listUSetDualSpace, copy_listMapDualSpace);
    //     }

    //     timeToPerform=duree(timeToPerform);

    //     monFlux1 << "s: "<< i<<", time: "<< timeToPerform<<std::endl;

    //     // batch

    //     cerr <<"/batch"<<endl;

    //     copy_donnees=donnees;   

    //     copy_topmost=topmost;

    //     copy_listUSetDualSpace=listUSetDualSpace;

    //     copy_listMapDualSpace=listMapDualSpace;
        
    //     timeToPerform=debut();

    //     BatchInsertSetOfTuples(new_tuples, copy_donnees, copy_topmost, d, copy_listUSetDualSpace, copy_listMapDualSpace);

    //     timeToPerform=duree(timeToPerform);
    
    //     monFlux1 << "b: "<< i<<", time: "<< timeToPerform<<std::endl; 

    //     // rebuild

    //     cerr <<"/rebuild"<<endl;

    //     copy_donnees=donnees;  

    //     for (int j = 0; j< i;j++){            
    //         new_tuples[j][0]=n+j;
    //         copy_donnees.push_back(new_tuples[j]);
    //     } 

    //     copy_topmost.clear();

    //     copy_listUSetDualSpace.clear();
    //     copy_listUSetDualSpace.resize(n+i);

    //     copy_listMapDualSpace.clear();
    //     copy_listMapDualSpace.resize(n+i);

    //     values_map topmost_map(d);

    //     timeToPerform=debut();

    //     NEG_wM::negativeSkycubeAux(copy_listUSetDualSpace, copy_listMapDualSpace, copy_donnees, copy_topmost, topmost_map, d);

    //     timeToPerform=duree(timeToPerform);

    //     monFlux1 << "r: "<< i<<", time: "<< timeToPerform<<std::endl;
    // }

    // // x12

    // string const nomFichier1("./experiments/insertion_NSCwM/ins-x12-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    // ofstream monFlux1(nomFichier1.c_str());
    // if (!monFlux1){
    //     cerr << "Impossible d'ouvrir le fichier" << nomFichier1 << endl;
    // }

    // cerr <<"-x12"<<endl;

    // int fin=100000;

    // int pas=10;

    // int deb=1000;
    
    // for (int i=deb;i<=fin;i=i*pas){
        
    //     cerr <<"= insertion: "<<i<<endl;

    //     TableTuple new_tuples;
        
    //     // generate new tuples
    //     loadData(dataName, "", i, d, k, new_tuples);

    //     for (int k=0 ; k< new_tuples.size();k++){
    //         new_tuples[k][0]=n+k;
    //     }

    //     // batch

    //     cerr <<"/batch"<<endl;

    //     TableTuple copy_donnees(donnees);

    //     TableTuple copy_topmost(topmost);

    //     std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);
        
    //     double timeToPerform=debut();

    //     BatchInsertSetOfTuples(new_tuples, copy_donnees, copy_topmost, d, copy_listUSetDualSpace, copy_listMapDualSpace);

    //     timeToPerform=duree(timeToPerform);
    
    //     monFlux1 << "b: "<< i<<", time: "<< timeToPerform<<std::endl; 

    //     // rebuild

    //     cerr <<"/rebuild"<<endl;

    //     copy_donnees=donnees;  

    //     for (int j = 0; j< i;j++){            
    //         new_tuples[j][0]=n+j;
    //         copy_donnees.push_back(new_tuples[j]);
    //     } 

    //     copy_topmost.clear();

    //     copy_listUSetDualSpace.clear();
    //     copy_listUSetDualSpace.resize(n+i);

    //     copy_listMapDualSpace.clear();
    //     copy_listMapDualSpace.resize(n+i);

    //     values_map topmost_map(d);

    //     timeToPerform=debut();

    //     NEG_wM::negativeSkycubeAux(copy_listUSetDualSpace, copy_listMapDualSpace, copy_donnees, copy_topmost, topmost_map, d);

    //     timeToPerform=duree(timeToPerform);

    //     monFlux1 << "r: "<< i<<", time: "<< timeToPerform<<std::endl;
    // }

    // // x3

    // string const nomFichier1("./experiments/insertion_NSCwM/ins-x3-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    // ofstream monFlux1(nomFichier1.c_str());
    // if (!monFlux1){
    //     cerr << "Impossible d'ouvrir le fichier" << nomFichier1 << endl;
    // }

    // cerr <<"-x3"<<endl;

    // int NSC_size=0;
    // for (int i=0; i<listMapDualSpace.size();i++) NSC_size+=listMapDualSpace[i].size();
    // monFlux1 << "begin, NSC_size: "<< NSC_size<<std::endl;    

    // int fin=1000;

    // int pas=10;

    // int deb=10;

    // TableTuple new_tuples_gen;
    
    // // generate new tuples

    // loadData(dataName, "", fin, d, k, new_tuples_gen);
    
    // for (int i=deb;i<=fin;i=i*pas){
        
    //     cerr <<"= insertion: "<<i<<endl;

    //     TableTuple new_tuples;

    //     for (int j=0; j<i; j++) new_tuples.push_back(new_tuples_gen[j]);

    //     for (int k=0 ; k< new_tuples.size();k++){
    //         new_tuples[k][0]=n+k;
    //     }

    //     // sequential

    //     cerr <<"/sequential"<<endl;

    //     TableTuple copy_donnees(donnees);

    //     TableTuple copy_topmost(topmost);

    //     std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);

    //     std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);

    //     double timeToPerform=debut();

    //     add_multiple_tuple_nothing(new_tuples, copy_donnees, copy_topmost, d, copy_listUSetDualSpace, copy_listMapDualSpace);

    //     timeToPerform=duree(timeToPerform);
        
    //     NSC_size=0;
    //     for (int i=0; i<copy_listMapDualSpace.size();i++) NSC_size+=copy_listMapDualSpace[i].size();

    //     monFlux1 << "s: "<< i<<", NSC_size: "<< NSC_size<<std::endl;
    //     monFlux1 << "s: "<< i<<", time: "<< timeToPerform<<std::endl;

    //     // batch with inclusion

    //     cerr <<"/batch i"<<endl;

    //     copy_donnees=donnees;   

    //     copy_topmost=topmost;

    //     copy_listUSetDualSpace=listUSetDualSpace;

    //     copy_listMapDualSpace=listMapDualSpace;
        
    //     timeToPerform=debut();

    //     add_multiple_tuple_inclusion(new_tuples, copy_donnees, copy_topmost, d, copy_listUSetDualSpace, copy_listMapDualSpace);

    //     timeToPerform=duree(timeToPerform);

    //     NSC_size=0;
    //     for (int i=0; i<copy_listMapDualSpace.size();i++) NSC_size+=copy_listMapDualSpace[i].size();
    
    //     monFlux1 << "bi: "<< i<<", NSC_size: "<< NSC_size<<std::endl; 
    //     monFlux1 << "bi: "<< i<<", time: "<< timeToPerform<<std::endl;

    //     // batch with fg

    //     cerr <<"/batch fg"<<endl;

    //     copy_donnees=donnees;   

    //     copy_topmost=topmost;

    //     copy_listUSetDualSpace=listUSetDualSpace;

    //     copy_listMapDualSpace=listMapDualSpace;
        
    //     timeToPerform=debut();

    //     BatchInsertSetOfTuples(new_tuples, copy_donnees, copy_topmost, d, copy_listUSetDualSpace, copy_listMapDualSpace);

    //     timeToPerform=duree(timeToPerform);

    //     NSC_size=0;
    //     for (int i=0; i<copy_listMapDualSpace.size();i++) NSC_size+=copy_listMapDualSpace[i].size();
    
    //     monFlux1 << "bfg: "<< i<<", NSC_size: "<< NSC_size<<std::endl;
    //     monFlux1 << "bfg: "<< i<<", time: "<< timeToPerform<<std::endl; 

    // fin x3

   // }

}

void insertion_option(string dataName,  DataType k, TableTuple& donnees, TableTuple& topmost, vector<USetDualSpace> &listUSetDualSpace, vector<mapDualSpace> &listMapDualSpace,Space d, NegSkyStr &structure0, map<DataType,DataType> &newIndexes0, map<DataType,DataType> &prvIndexes0){

    cerr <<"**** Begin: Simple insertion NSCwM"<<endl;

    int n = donnees.size();

    TableTuple new_tuples;

    loadData(dataName, "", 2, d, k, new_tuples);

    //int id_new_tuple;

    //int tuple[15]={1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000};

    //cout << "chooose id of the new tuple: ";

    //cin >> id_new_tuple;

    //new_tuple[0]=tuple; // on lui donne un identifiant inexistant auparavant

    // new_tuple[0][0]=donnees.size();

    // structure0.clear();

    // newIndexes0.clear();

    // prvIndexes0.clear();

    // double timeToPerform=debut();

    // InsertTuple(new_tuple, donnees, topmost, d, listUSetDualSpace, listMapDualSpace);

    // timeToPerform=duree(timeToPerform);

    // std::cout << "Adding tuple, time: "<< timeToPerform<<std::endl;

    // timeToPerform=debut();

    // NEG::negativeSkycube(structure0, newIndexes0, prvIndexes0, listUSetDualSpace, d);

    // timeToPerform=duree(timeToPerform);

    // std::cout << "fin SH, time: "<< timeToPerform<<std::endl; 

    //matrixQuery(pairMAP);

    //x2

    string const nomFichier1("./experiments/insertion_NSCwM/ins-x2-"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n));
    ofstream monFlux1(nomFichier1.c_str());
    if (!monFlux1){
        cerr << "Impossible d'ouvrir le fichier" << nomFichier1 << endl;
    }

    cerr <<"-x2"<<endl;

    for (int i=0; i<new_tuples.size(); i++){

    TableTuple new_tuple;

    new_tuple.push_back(new_tuples[i]);

    new_tuple[0][0]=donnees.size();

    // with add tuple

    cerr <<"/InsertTuple"<<endl;

    TableTuple copy_donnees(donnees);
    TableTuple copy_topmost(topmost);
    std::vector<USetDualSpace> copy_listUSetDualSpace(listUSetDualSpace);
    std::vector<mapDualSpace> copy_listMapDualSpace(listMapDualSpace);    

    double timeToPerform=debut();

    InsertTuple(new_tuple, copy_donnees, copy_topmost, d, copy_listUSetDualSpace, copy_listMapDualSpace);

    timeToPerform=duree(timeToPerform);

    monFlux1 << "a: "<<", time: "<< timeToPerform<<std::endl;

    // with rebuild

    cerr <<"/rebuild"<<endl;

    copy_donnees=donnees;  
    new_tuple[0][0]=n;
    copy_donnees.push_back(new_tuple[0]);
    copy_topmost.clear();
    copy_listUSetDualSpace.clear();
    copy_listUSetDualSpace.resize(n+1);
    copy_listMapDualSpace.clear();
    copy_listMapDualSpace.resize(n+1);
    values_map topmost_map(d);
    
    timeToPerform=debut();
    
    NEG_wM::negativeSkycubeAux(copy_listUSetDualSpace, copy_listMapDualSpace, copy_donnees, copy_topmost, topmost_map, d);
    
    timeToPerform=duree(timeToPerform);
    
    monFlux1 << "r: "<<", time: "<< timeToPerform<<std::endl;

    }
}


}
#endif // NEGATIVEWM_H_INCLUDED