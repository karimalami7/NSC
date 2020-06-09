#ifndef QUERIES_H_INCLUDED
#define QUERIES_H_INCLUDED


void spaces_dominated_for_tuples(vector<unordered_set<int>> &spaces_dominated, NegSkyStr &structureNSC, Space d)
{
    for (auto itXY=structureNSC.rbegin();itXY!=structureNSC.rend();++itXY){
        for (auto itY=(itXY->second).begin();itY!=(itXY->second).end();++itY) {
            vector<Space> listCouv;
            listCouverts(itXY->first-itY->first, itY->first, d, listCouv);
            for (int id : itY->second){
                spaces_dominated[id].insert(listCouv.begin(), listCouv.end());
            }
        }
    }
}

// Rank the tuples by their skyline frequency, i.e, the number of skyline to which they belong 
void rank_F(string dataName, TableTuple &donnees, NegSkyStr &structureNSC, int data_size, DataType k, Space d, int set_size){

    Space full_space=(1<<d)-1;

    //***********************************************//
    // rank tuples

    double timeToPerform=debut();
    vector<unordered_set<int>> spaces_dominated(data_size);
    spaces_dominated_for_tuples(spaces_dominated, structureNSC, d);


    vector<pair <int,int>> score(data_size);
    for (int i=0; i<data_size; i++){ // for each tuple
        score[i]=make_pair(full_space-spaces_dominated[i].size(),i);
    }
    sort(score.begin(), score.end());
    timeToPerform=duration(timeToPerform);

    point_set_t* result_set=alloc_point_set(set_size);
    int n=donnees.size();
    for (int i = n-1; i >= n-set_size; i--)
    {   
        point_t* p = alloc_point(d, n-1-i);
        for (int j = 0; j < d; j++)
        {   
            p->coord[j]=1-((double)donnees[score[i].second][j+1]/k);
        }
        result_set->points[n-1-i] = p;
    }
    
    //***********************************************//


    //***********************************************//
    // compute the regret of the result set
    
    //1 get topmost ids
    bool* topmost=NEG::subspaceSkyline_NSC(structureNSC, donnees.size(), full_space);

    //2 create topmost dataset
    int topmost_size=0;
    for (int i=0; i< donnees.size(); i++) if (topmost[i]) topmost_size++;
    point_set_t* topmost_set=alloc_point_set(topmost_size);
    
    int l=0;
    for (int i = 0; i < donnees.size(); i++)
    {   
        if (topmost[i])
        {
            point_t* p = alloc_point(d, l);
            for (int j = 0; j < d; j++) p->coord[j]=1-((double)donnees[i][j+1]/k);
            topmost_set->points[l] = p;
            l++;
        }
    }

    //3 evaluate the regret
    double mrr = evaluateLP(topmost_set, result_set, 0);

    printf("%15s%15f%15f\n", "Frequency", mrr, timeToPerform);    

    // int n=donnees.size();
    // string const nomFichier1("../datasets/"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n)+"-k-"+std::to_string(k)+"-F");
    // ofstream monFlux1(nomFichier1.c_str());
    // if (monFlux1){
    //     for(int i=n-1; i>=0; --i){
    //         // monFlux1 <<"id:"<<score[i].second<<" ";
    //         monFlux1 <<"score:"<<score[i].first<<" ";
    //         for (Space j=1;j<=d;j++){
    //             monFlux1<<donnees[score[i].second][j]<<" ";
    //         }
    //         //for (auto subspace: spaces_dominated[score[i].second]) monFlux1 << subspace<<",";
    //         monFlux1 <<endl;
    //     }
    //     monFlux1.close();
    // }
    // else{
    //     cout << "ERROR: Couldn't open the file." << endl;
    // }

    // displayResult(dataName, valid_data_size, d, k, "topK" , score[0], duree(timeToPerform), "topK");

}

// Rank the tuples by their skyline width, i.e, the sum of the width of the skyline  
void rank_D(string dataName, TableTuple &donnees, NegSkyStr &structureNSC, int data_size, DataType k, Space d, int set_size){


    //***********************************************//
    // rank tuples

    double timeToPerform=debut();
    vector<unordered_set<int>> spaces_dominated(data_size);
    spaces_dominated_for_tuples(spaces_dominated, structureNSC, d);

    vector<pair <int,int>> score(data_size);
    for (int i=0; i<data_size; i++){ // for each tuple
        int tuple_score=0;
        for(auto subspace: spaces_dominated[i]){
            tuple_score+=spaceSize(subspace)+1;
        }
        score[i]=make_pair(tuple_score,i);
    }
    sort(score.begin(), score.end());
    timeToPerform=duration(timeToPerform);

    point_set_t* result_set=alloc_point_set(set_size);
    int n=donnees.size();
    for (int i = 0; i < set_size; i++)
    {   
        point_t* p = alloc_point(d, i);
        for (int j = 0; j < d; j++)
        {   
            p->coord[j]=1-((double)donnees[score[i].second][j+1]/k);
        }
        result_set->points[i] = p;
    }

    //***********************************************//

    //***********************************************//
    // compute the regret of the result set
    
    //1 get topmost ids
    Space full_space=(1<<d)-1;
    bool* topmost=NEG::subspaceSkyline_NSC(structureNSC, donnees.size(), full_space);

    //2 create topmost dataset
    int topmost_size=0;
    for (int i=0; i< donnees.size(); i++) if (topmost[i]) topmost_size++;
    point_set_t* topmost_set=alloc_point_set(topmost_size);
    
    int l=0;
    for (int i = 0; i < donnees.size(); i++)
    {   
        if (topmost[i])
        {
            point_t* p = alloc_point(d, l);
            for (int j = 0; j < d; j++) p->coord[j]=1-((double)donnees[i][j+1]/k);
            topmost_set->points[l] = p;
            l++;
        }
    }

    //3 evaluate the regret
    double mrr = evaluateLP(topmost_set, result_set, 0);

    printf("%15s%15f%15f\n", "Dimensionality", mrr, timeToPerform); 


    // int n=donnees.size();
    // string const nomFichier1("../datasets/"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n)+"-k-"+std::to_string(k)+"-D");
    // ofstream monFlux1(nomFichier1.c_str());
    // if (monFlux1){
    //     for(int i=0; i<n; ++i){
    //         // monFlux1 <<"id:"<<score[i].second<<" ";
    //         monFlux1 <<"score:"<<score[i].first<<" ";
    //         for (Space j=1;j<=d;j++){
    //             monFlux1<<donnees[score[i].second][j]<<" ";
    //         }
    //         monFlux1 <<endl;
    //     }
    //     monFlux1.close();
    // }
    // else{
    //     cout << "ERROR: Couldn't open the file." << endl;
    // }

    // displayResult(dataName, valid_data_size, d, k, "topK" , score[0], duree(timeToPerform), "topK");

}

void rank_P(string dataName, TableTuple &donnees, NegSkyStr &structureNSC, int data_size, DataType k, Space d){


    Space All=(1<<d)-1; 

    // foreach subspace 

    vector<unordered_set<int>> spaces_dominated(data_size);
    for (auto itXY=structureNSC.rbegin();itXY!=structureNSC.rend();++itXY){
        for (auto itY=(itXY->second).begin();itY!=(itXY->second).end();++itY) {
            vector<Space> listCouv;
            listCouverts(itXY->first-itY->first, itY->first, d, listCouv);
            for (int id : itY->second){
                spaces_dominated[id].insert(listCouv.begin(), listCouv.end());
            }
        }
    }


    vector <int> smallest_skyline_space(data_size, d+1);

    for (int subspace=1; subspace<=All ; subspace++ ){
        for (int i=0; i<data_size; i++){
            if (spaces_dominated[i].find(subspace)==spaces_dominated[i].end()){
                // *it is skyline on subspace
                int size=spaceSize(subspace);
                if (smallest_skyline_space[i]>size){
                    smallest_skyline_space[i]=size;
                }
            }
        }
    }


    vector<pair <int,int>> score(data_size);
    for (int i=0; i<data_size; i++){
        score[i]=make_pair(smallest_skyline_space[i],i);
    }

    sort(score.begin(), score.end());

    int n=donnees.size();
    string const nomFichier1("../datasets/"+dataName+"-d-"+std::to_string(d)+"-n-"+std::to_string(n)+"-k-"+std::to_string(k)+"-P");
    ofstream monFlux1(nomFichier1.c_str());
    if (monFlux1){
        for(int i=0; i<n; ++i){
            monFlux1 <<"score:"<<score[i].first<<" ";
            for (Space j=1;j<=d;j++){
                monFlux1<<donnees[score[i].second][j]<<" ";
            }
            monFlux1 <<endl;
        }
        monFlux1.close();
    }
    else{
        cout << "ERROR: Couldn't open the file." << endl;
    }

    // displayResult(dataName, valid_data_size, d, k, "topK" , score[0], duree(timeToPerform), "topK");

}

void representative_by_sphere(string dataName, TableTuple &donnees, NegSkyStr &structureNSC, DataType k, Space d, int output_size){

    
    //1 define input as topmost or the whole dataset
    Space full_space=(1<<d)-1;
    bool* input;

    cout << "Which input you want: "<<endl;
    cout << "1- The whole dataset: "<<endl;
    cout << "2- Skyline: "<<endl;
    cout << "3- K frequent: "<<endl;
    cout << "4- K width: "<<endl;
    cout << "5- K priority: "<<endl;
    int option;
    cin >> option;
    int input_size;
    if (option==3 || option==4 || option==5){
        cout << "give input size: "<<endl;
        cin >> input_size;
    }
    point_set_t* input_set;
    double timeToPrepareInput=debut();
    switch(option){
        case 1:
            //*********************************************//
            {input=new bool[donnees.size()];for (int i=0;i<donnees.size();i++) input[i]=true;
            input_size=0;
            for (int i=0; i< donnees.size(); i++) if (input[i]) input_size++;
            input_set=alloc_point_set(input_size);
            
            int l=0;
            for (int i = 0; i < donnees.size(); i++)
            {   
                if (input[i])
                {   
                    point_t* p = alloc_point(d, l);
                    for (int j = 0; j < d; j++) p->coord[j]=1-((double)donnees[i][j+1]/k);
                    input_set->points[l]=p;
                    l++;
                }
            }
            break;}
            //*********************************************//
        case 2:
            //*********************************************//
            {double timeToPrepareInput=debut();
            input=NEG::subspaceSkyline_NSC(structureNSC, donnees.size(), full_space);
            timeToPrepareInput=duration(timeToPrepareInput);
            printf("%15s%15f\n", "Sky", timeToPrepareInput);

            input_size=0;
            for (int i=0; i< donnees.size(); i++) if (input[i]) input_size++;
            input_set=alloc_point_set(input_size);
            
            int l=0;
            for (int i = 0; i < donnees.size(); i++)
            {   
                if (input[i])
                {   
                    point_t* p = alloc_point(d, l);
                    for (int j = 0; j < d; j++) p->coord[j]=1-((double)donnees[i][j+1]/k);
                    input_set->points[l]=p;
                    l++;
                }
            }        
            break;}
            //*********************************************//
        case 3:
            //*************************************//
            {double timeToPrepareInput=debut();
            vector<unordered_set<int>> spaces_dominated(donnees.size());
            spaces_dominated_for_tuples(spaces_dominated, structureNSC, d);

            vector<pair <int,int>> score(donnees.size());
            for (int i=0; i<donnees.size(); i++){ // for each tuple
                score[i]=make_pair(spaces_dominated[i].size(),i);
            }
            sort(score.begin(), score.end());
            timeToPrepareInput=duration(timeToPrepareInput);
            printf("%15s%15f\n", "Rank by F", timeToPrepareInput);

            input_set=alloc_point_set(input_size);
            int n=donnees.size();
            for (int i = 0; i < input_size; i++)
            {   
                point_t* p = alloc_point(d, i);
                for (int j = 0; j < d; j++)
                {   
                    p->coord[j]=1-((double)donnees[score[i].second][j+1]/k);
                }
                input_set->points[i] = p;
            }
            break;}
            //*************************************//
        case 4:
            //*************************************//
            {double timeToPrepareInput=debut();
            vector<unordered_set<int>> spaces_dominated(donnees.size());
            spaces_dominated_for_tuples(spaces_dominated, structureNSC, d);

            vector<pair <int,int>> score(donnees.size());
            for (int i=0; i<donnees.size(); i++){ // for each tuple
                int tuple_score=0;
                for(auto subspace: spaces_dominated[i]){
                    tuple_score+=spaceSize(subspace)+1;
                }
                score[i]=make_pair(tuple_score,i);
            }
            sort(score.begin(), score.end());
            timeToPrepareInput=duration(timeToPrepareInput);
            printf("%15s%15f\n", "Rank by D", timeToPrepareInput);

            input_set=alloc_point_set(input_size);
            int n=donnees.size();
            for (int i = 0; i < input_size; i++)
            {   
                point_t* p = alloc_point(d, i);
                for (int j = 0; j < d; j++)
                {   
                    p->coord[j]=1-((double)donnees[score[i].second][j+1]/k);
                }
                input_set->points[i] = p;
            }
            break;}
            //*************************************//
        case 5:
            //*************************************//
            {double timeToPrepareInput=debut();
            vector<unordered_set<int>> spaces_dominated(donnees.size());
            spaces_dominated_for_tuples(spaces_dominated, structureNSC, d);
            vector <int> smallest_skyline_space(donnees.size(), d+1);
            for (int subspace=1; subspace<=full_space ; subspace++ ){
                for (int i=0; i<donnees.size(); i++){
                    if (spaces_dominated[i].find(subspace)==spaces_dominated[i].end()){
                        // *it is skyline on subspace
                        int size=spaceSize(subspace);
                        if (smallest_skyline_space[i]>size){
                            smallest_skyline_space[i]=size;
                        }
                    }
                }
            }
            vector<pair <int,int>> score(donnees.size());
            for (int i=0; i<donnees.size(); i++){
                score[i]=make_pair(smallest_skyline_space[i],i);
            }
            sort(score.begin(), score.end());
            timeToPrepareInput=duration(timeToPrepareInput);
            printf("%15s%15f\n", "Rank by P", timeToPrepareInput);

            input_set=alloc_point_set(input_size);
            int n=donnees.size();
            for (int i = 0; i < input_size; i++)
            {   
                point_t* p = alloc_point(d, i);
                for (int j = 0; j < d; j++)
                {   
                    p->coord[j]=1-((double)donnees[score[i].second][j+1]/k);
                }
                input_set->points[i] = p;
            }
            break;}
            //*************************************//
        default: exit(0);break;
    }


    //3 run sphere
    point_set_t* S;
    double timeToPerformSphere=debut();
    S = sphereWSImpLP(input_set, output_size);
    timeToPerformSphere=duration(timeToPerformSphere);


    //***********************************************//
    // compute the regret of the result set

    //1 get topmost ids
    bool* topmost=NEG::subspaceSkyline_NSC(structureNSC, donnees.size(), full_space);

    //2 create topmost dataset
    int topmost_size=0;
    for (int i=0; i< donnees.size(); i++) if (topmost[i]) topmost_size++;
    point_set_t* topmost_set=alloc_point_set(topmost_size);
    
    int l=0;
    for (int i = 0; i < donnees.size(); i++)
    {   
        if (topmost[i])
        {
            point_t* p = alloc_point(d, l);
            for (int j = 0; j < d; j++) p->coord[j]=1-((double)donnees[i][j+1]/k);
            topmost_set->points[l] = p;
            l++;
        }
    }

    //double mrr = evaluateLP(topmost_set, S, 0);
    double mrr = 0;
    printf("%15s%15f%15f\n", "Sphere", mrr, timeToPerformSphere);
}

void skylinequery(string dataName, TableTuple &donnees, NegSkyStr structure0, Space d, DataType k,vector<Space> &subspaceN, vector<Space> &subspaceAll){

    Space subspace;

    cerr << "Choose a subspace for the skyline query: ";

    cin >> subspace;

    //query on 1 subspace

    int structSize=0;
    double timeToPerform=debut();
    structSize+=NEG::subspaceSkylineSize_NSC(structure0, donnees.size(), subspace);
    timeToPerform=duree(timeToPerform);
    NEG::displayResultv2(dataName, donnees.size(), d, k, "N=1", structSize, timeToPerform, NSC); 

    //query on all subspaces

    Space All=(1<<d)-1;
    structSize=0;
    timeToPerform=debut();
    for (int i=0;i<All;i++){
        structSize+=NEG::subspaceSkylineSize_NSC(structure0, donnees.size(), subspaceAll[i]);
    }
    timeToPerform=duree(timeToPerform);
    NEG::displayResultv2(dataName, donnees.size(), d, k, "SKYCUBE", structSize, timeToPerform, NSC);
}

#endif