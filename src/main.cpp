//définit le nombre de threads à exécuter en parallèle (si NB_THREADS==1 alors exécution séquentielle)
int NB_THREADS=22;
int K_MIN=1;

#include "generator/DataGenerator.h"
#include "generator/generateur.h"
#include "generator/generateKossmann.h"
#include "common/affichage.h"
#include "common/skyline.h"
#include "bskytree/bskytree.h"
#include "csc/compressed.h"
#include "nsc/algoglouton.h"
#include "nsc/negative.h"
#include "nsc/negative_wM.h"
#include "kdominance/kdomskycube.h"
#include "sphere/sphere.h"
#include "nsc/queries.h"
#include "experimentations.h"



int main(int argc,char **argv){

    //syntaxe d'appel      ./cmd "ANTI"|"INDE"|"CORR"|"PERS"|DataName   k|Path    n   d    NB_THREADS    [NAIF|TREE|NSC|CSC]*
    //DataName est le nom qui apparaitra dans le résultat tandis que Path est le chemin vers le fichier de données
    //lorsque l'une des option "ANTI"|"INDE"|"CORR" est choisie alors k doit être spécifié
    //lorsque DataName est spécifié alors Path doit également être spécifié
    //NB_THREADS==1 -> Sequentiel ; NB_THREADS==m -> m threads créés (uniquement dans la phase de pré-calcul)
    //si aucune méthode n'est spécifiée alors toutes les méthodes sont exécutées

    //PERS est une version de données ANTICORRELEES avec deux groupes de variables, les variables sont corrélés intra-groupe et anticorrélés inter-groupes
    //les groupes sont respectivement d'indice pair et impair

    /*Space nbDim=4;
    for (Space level=1; level<nbDim;level++){
        vector<Space> result;
        listSpaceLevel(nbDim, level, result);
        cout<<level<<"("<<result.size()<<"): ";
        for (int i=0;i<result.size();i++) {cout<<" ";displaySubspace(result[i], nbDim);}
        cout<<endl;
    }
    return 0;*/

    time_t val=1448639227;
    srand (val);
//srand (time(NULL));
//    string dataName="ANTI";//"baseball";
//    string path="";//"data/Baseball";
//    DataType k=100;
//    DataType n=10000;
//    Space d=10;
//    bool selectedMethod[]={
//      false // naif (Non Parallelisable)
//    , false // tree (Non Parallelisable)
//    , false // CSC (Parallelisable)
//    , false // BUA
//    , false // UBA
//    , false // DPI
//    , true // NSC (Parallelisable)  voir option NB_THREADS au début de ce fichier
//    , true // NSC2 (Parallelisable)  voir option NB_THREADS au début de ce fichier
//    , false // NSCk
//    };

    string dataName="INDE";//"baseball";
    string path="";//"data/Baseball";
    DataType k=100;
    DataType n=1000;
    Space d=14;
    bool selectedMethod[]={
      false // naif (Non Parallelisable)
    , true // tree (Non Parallelisable)
    , false // CSC (Parallelisable)
    , false // BUA
    , false // UBA
    , false // DPI
    , true // NSC (Parallelisable)  voir option NB_THREADS au début de ce fichier
    , true // NSC2 (Parallelisable)  voir option NB_THREADS au début de ce fichier
    , false // NSCk
    , true // NSC3 calcul du skycube sans passer par la structure de données intermédiaire
    , false // NSCv2
    , false // NSCwM
    };

    int i, j;
    NB_THREADS=22;//C'est le nombre de threads qu'on veut lancer en parallèle.
    if (argc==1){
        experimentSkycube(dataName, path, k, n, d, selectedMethod);
    }else{
        dataName=argv[1];
        if (dataName=="ANTI"||dataName=="INDE"||dataName=="CORR"||dataName=="PERS"){
            k=atoi(argv[2]);
            path="";
        }else{
            k=0;
            path=argv[2];
        }
        n=atoi(argv[3]);
        d=atoi(argv[4]);
        NB_THREADS=atoi(argv[5]);
        //NB_THREADS=2;
        if (argc>6) for (j=0;j<NB_METHOD;j++) selectedMethod[j]=false;
        for (i=6;i<argc;i++) for (j=0;j<NB_METHOD;j++) if (argv[i]==methodNames[j]) selectedMethod[j]=true;
        experimentSkycube(dataName, path, k, n, d, selectedMethod);
    }
    //cout<<endl<<"seed: "<<val<<endl;
    return 0;
}
