#ifndef COMPRESSED_H_INCLUDED
#define COMPRESSED_H_INCLUDED

#include "../common/declarations.h"

Space COMPRESSED_DIM;

void sortIndexComp(const vector<Point> &d, vector<DataType> &result) {
    vector<DataType> idx(d.size());
    DataType i;
    for (i = 0; i != (DataType)idx.size(); i++) idx[i] = i;
    sort(idx.begin(), idx.end(), [&d](DataType i1, DataType i2) {return d[i1][COMPRESSED_DIM] < d[i2][COMPRESSED_DIM];});
    result.swap(idx);
}

void listSpaceLevelR(set<Space>& setSp, Space d, Space level, Space subs, Space dd, Space jeton){
    if (jeton==0){
        setSp.insert(subs<<dd);
    }else{
        if (jeton<dd) listSpaceLevelR(setSp, d, level, subs*2, dd-1, jeton);
        else{
            setSp.insert((subs<<dd)+(1<<dd)-1);
        }
        if (jeton>0) listSpaceLevelR(setSp, d, level, subs*2+1, dd-1, jeton-1);
    }
}

void listSpaceLevel(Space d, Space level, vector<Space> &result){
    set<Space> setSp;
    vector<Space> vect;
    listSpaceLevelR(setSp,d,level,0,d,level);
    for (auto it=setSp.begin();it!=setSp.end();it++) vect.push_back(*it);
    result.swap(vect);
}

void collecteTuples_TD(USetId &tabId, TableTuple &donnees, VectorUSetId &structure, vector<Space> vectSpace, Space realSubspace, Space originSpace, Space actualSpace, vector<Space> &missingAttrib, Space nbDiscarded, Space d){
/*  Cette procédure marque à 'false' tous les tuples retrouvées dans tous les sur-espaces de 'actualSpace'
    missingAttribPow est la liste de attribut absents dans 'actualSpace' et qui seront utilisés pour accéder aux sur-espaces */
    Space i;
    TableTuple skyline, listTuple;
    USetId::iterator itId;
    if (structure[actualSpace].size()>0){
        for (itId=structure[actualSpace].begin();itId!=structure[actualSpace].end();itId++){
            listTuple.push_back(donnees[*itId]);
        }
        ExecuteBSkyTree(vectSpace, listTuple, skyline);
    }

    TableTuple::iterator it=skyline.begin();
    while (it!=skyline.end()){
        tabId.insert((*it)[0]);
        it++;
    }
    for (i=nbDiscarded;i<(Space)missingAttrib.size();i++) collecteTuples_TD(tabId, donnees, structure, vectSpace, realSubspace, originSpace, actualSpace+(1<<(missingAttrib[i]-1)), missingAttrib, i+1, d);
}

void queryCSC(TableTuple &donnees, VectorUSetId &structure, DataType n, Space d, Space subspace, USetId &resultat){
    if (subspace==((1<<d)-1)){
        resultat=structure[subspace];
        return;
    }
    USetId result;
    vector<Space> listAttrib;
    vector<Space> vectSpace;
    listeAttributsPresents(subspace, d, listAttrib);
    vectSpace=listAttrib;
    collecteTuples_TD(result, donnees, structure, vectSpace, subspace, 0, 0, listAttrib, 0, d);
    resultat.swap(result);
}

DataType subspaceSkylineSize_CSC(TableTuple &donnees, VectorUSetId &structure, DataType n, Space d, Space subspace){
    USetId result;
    queryCSC(donnees, structure, n, d, subspace, result);
    return result.size();
}

bool isDominated(TableTuple& donnees, DataType p, vector<Space> &listAttrib){
    Space j, nbI, nbE;
    bool dominated=false;
    DataType i=0;
    while (!dominated && i<(DataType)donnees.size()){
        nbI=0;
        nbE=0;
        for (j=0;j<(Space)listAttrib.size();j++)
            if (donnees[i][listAttrib[j]]<donnees[p][listAttrib[j]]) nbI++;
            else if (donnees[i][listAttrib[j]]==donnees[p][listAttrib[j]]) nbE++;
        if (nbI+nbE==(Space)listAttrib.size() && nbI>0) dominated=true;
        i++;
    }
    return dominated;
}

bool isDominated2(TableTuple& donnees, Point p, vector<Space> &listAttrib){
    Space j, nbI, nbE;
    bool dominated=false;
    DataType i=0;
    while (!dominated && i<(DataType)donnees.size()){
        nbI=0;
        nbE=0;
        for (j=0;j<(Space)listAttrib.size();j++)
            if (donnees[i][listAttrib[j]]<p[listAttrib[j]]) nbI++;
            else if (donnees[i][listAttrib[j]]==p[listAttrib[j]]) nbE++;
        if (nbI+nbE==(Space)listAttrib.size() && nbI>0) dominated=true;
        i++;
    }
    return dominated;
}

long compressedSkycube_PARALLELE(VectorUSetId &structure, TableTuple& donnees, Space d){//BUS_CSC
/*  Cette procédure construit à partir des données la structure de compression de skycube suivant l'idéologie COMPRESSED SKYCUBE*/
    long structSize=0;
    Space i;
    DataType n=donnees.size();
    vector<DataType>* sortedDim=new vector<DataType>[d];

    for (auto j=1;j<=d;j++){
        COMPRESSED_DIM=j;
        sortIndexComp(donnees, sortedDim[j-1]);
    }

    for (i=1;i<=d;i++){
        //vector<Space>::iterator it;
        vector<Space> listSpace;
        listSpaceLevel(d, i, listSpace);

        vector<vector<Space>> listAttrib(listSpace.size());
        vector<USetId> fs(listSpace.size());
        vector<vector<DataType>> CU(listSpace.size());

        #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic)
        for (auto ii=0;ii<(Space)listSpace.size();ii++){
            vector<DataType> f;
            listeAttributsPresents(listSpace[ii], d, listAttrib[ii]);
            queryCSC(donnees, structure, n, d, listSpace[ii], fs[ii]);
            for (auto itfs=fs[ii].begin(); itfs!=fs[ii].end(); itfs++){
                DataType aux=0;
                for (auto j=0;j<(Space)listAttrib[ii].size();j++) aux+=donnees[(*itfs)][listAttrib[ii][j]];
                f.push_back(aux);
            }
            for (auto k=0;k<n;k++){
                DataType t=sortedDim[listAttrib[ii][0]-1][k];
                //t=sortedDim[listAttrib[ii][0]-1][k];
                auto itfs=fs[ii].find(t);
                if (itfs!=fs[ii].end()) fs[ii].erase(itfs);
                else{
                    if (!isDominated(donnees, t, listAttrib[ii])){
                        fs[ii].insert(t);
                        CU[ii].push_back(t);
                    }
                }
            }
        }
        for (auto ii=0;ii<(Space)listSpace.size();ii++){structSize+=CU[ii].size();
            for (auto jj=0;jj<(Space)CU[ii].size();jj++) structure[listSpace[ii]].insert(CU[ii][jj]);
        }
    }
    delete []sortedDim;
    return structSize;
}

long compressedSkycube_SEQUENTIEL(VectorUSetId &structure, TableTuple& donnees, Space d){//BUS_CSC
/*  Cette procédure construit à partir des données la structure de compression de skycube suivant l'idéologie COMPRESSED SKYCUBE*/
    long structSize=0;
    Space i, j;
    DataType n=donnees.size(), aux, k, t;
    vector<DataType>* sortedDim=new vector<DataType>[d];
    vector<Space> listSpace;
    vector<Space>::iterator it;
    USetId::iterator itfs;
    USetId fs;
    vector<DataType> f;
    vector<Space> listAttrib;

    for (j=1;j<=d;j++){
        COMPRESSED_DIM=j;
        sortIndexComp(donnees, sortedDim[j-1]);
    }

    for (i=1;i<=d;i++){
        listSpace.clear();
        listSpaceLevel(d, i, listSpace);
        for (it=listSpace.begin();it!=listSpace.end();it++){
            listeAttributsPresents(*it, d, listAttrib);
            queryCSC(donnees, structure, n, d, *it, fs);
            for (itfs=fs.begin();itfs!=fs.end();itfs++){
                aux=0;
                for (j=0;j<(Space)listAttrib.size();j++) aux+=donnees[(*itfs)][listAttrib[j]];
                f.push_back(aux);
            }
            for (k=0;k<n;k++){
                t=sortedDim[listAttrib[0]-1][k];
                itfs=fs.find(t);
                if (itfs==fs.end()){//f.erase(itfs);  else
                    if (!isDominated(donnees, t, listAttrib)){
                        //fs.insert(t);
                        structure[*it].insert(t);
                        structSize++;
                    }
                }
            }
            f.clear();
        }
    }
    delete []sortedDim;
    return structSize;
}


void afficheVectorTableTuple(VectorUSetId &structure, Space d){
    long nbSubspace=1<<d;
    USetId::iterator it;
    Space i;
    for (i=1;i<nbSubspace;i++){
        displaySubspace(i, d);
        cout<<": ";
        for (it=structure[i].begin();it!=structure[i].end();it++){
            cout<<"t"<<(*it)+1<<" ";
        }
        cout<<endl;
    }
}

/**************************************************DEBUT NOUVELLE VERSION DE COMPRESSED SKYCUBE***********************************************************************/

void collecteTuples_TD(unordered_set<DataType> &uSetId, VectorUSetId &structure, vector<Space> vectSpace, Space realSubspace, Space originSpace, Space actualSpace, vector<Space> &missingAttrib, Space nbDiscarded, Space d){
/*  Cette procédure collecte et stocke dans uSetId tous les tuples retrouvées dans tous les sur-espaces de 'actualSpace'
    missingAttribPow est la liste de attribut absents dans 'actualSpace' et qui seront utilisés pour accéder aux sur-espaces
    en partant de originSpace=0 et spécifiant missing attrib comme une liste d'attributs on a la liste des tuples des sous-espaces de l'espace défini par la liste spécifiée à missingAttrib
*/
    Space i;
    TableTuple skyline, listTuple;
    USetId::iterator itId;
    for (itId=structure[actualSpace].begin();itId!=structure[actualSpace].end();itId++){
        uSetId.insert(*itId);
    }
    for (i=nbDiscarded;i<(Space)missingAttrib.size();i++) collecteTuples_TD(uSetId, structure, vectSpace, realSubspace, originSpace, actualSpace+(1<<(missingAttrib[i]-1)), missingAttrib, i+1, d);
}

long compressedSkycube2(VectorUSetId &structure, TableTuple& donnees, Space d){
    Space All=(1<<d)-1;
    vector<Space> subspaceAll;
    vector<vector<Space>> listAllTabSpace(All);
    for (auto i=1;i<=All;i++){
        subspaceAll.push_back(i);
        listeAttributsPresents(i, d, listAllTabSpace[i-1]);
    }
    long structSize=0;

    for (auto i=1;i<=d;i++){
        vector<Space> listSpace;
        listSpaceLevel(d, i, listSpace);

        #pragma omp parallel for num_threads(NB_THREADS) schedule(dynamic) reduction(+:structSize)
        for (auto ii=0;ii<listSpace.size();ii++){
            Space space=listSpace[ii];
            TableTuple skySpace;
            ExecuteBSkyTree(listAllTabSpace[space-1], donnees, skySpace);

            vector<DataType> idSkySpace(skySpace.size());
            for (auto iii=0;iii<skySpace.size();iii++) idSkySpace[iii]=skySpace[iii][0];

            unordered_set<DataType> sASupprimer;
            vector<DataType> vASupprimer;

            if (space!=All){/*On élimine rien pour le topmost*/
                collecteTuples_TD(sASupprimer, structure, listAllTabSpace[space-1], space, 0, 0, listAllTabSpace[space-1], 0, d);

                for (auto it=sASupprimer.begin();it!=sASupprimer.end();it++) vASupprimer.push_back(*it);
            }
            sort(idSkySpace.begin(),idSkySpace.end());
            sort(vASupprimer.begin(), vASupprimer.end());

            vector<DataType> difference(idSkySpace.size());
            vector<DataType>::iterator it = set_difference(idSkySpace.begin(),idSkySpace.end(), vASupprimer.begin(), vASupprimer.end(), difference.begin());
            difference.resize(it-difference.begin());
            for (auto j=0;j<difference.size();j++) structure[space].insert(difference[j]);

            structSize+=difference.size();
        }
    }
    //afficheVectorTableTuple(structure, d);
    return structSize;
}


/**************************************************FIN NOUVELLE VERSION DE COMPRESSED SKYCUBE***********************************************************************/


long compressedSkycube(VectorUSetId &structure, TableTuple& donnees, Space d){
    if (NB_THREADS>1) return compressedSkycube_PARALLELE(structure, donnees, d);
    else return compressedSkycube_SEQUENTIEL(structure, donnees, d);
}

#endif // COMPRESSED_H_INCLUDED
