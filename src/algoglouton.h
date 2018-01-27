#ifndef ALGOGLOUTON_H_INCLUDED
#define ALGOGLOUTON_H_INCLUDED

#include "declarations.h"

void set_of_bits(int i, 	vector<int> &V){
    //on récupère les position des bits à 1 dans i qu'on met dans V
    //Par ex: si i = 101 (en binaire) alors V=[1, 3], i.e., les bits 1 et 3 sont à 1
    int j=__builtin_ffs(i);//on récupère la position du bit le plus faible positionné à 1	
    while(i!=0){
        if (j==1) {
            i=i-1;
            V.push_back(j);
        }
        else{
            i = (i & (~(1<<(j-1))));//on met à 0 le bit j
            V.push_back(j);
        }
        j=__builtin_ffs(i);//on récupère la position du prochain bit à 1
    }
}

void mis(vector<int> &f, int i, vector<vector<int> > &V1)
{
	//on génère les sous-ensemble de f qu'on met dans V1
	// Ex: si f=[1, 3], alors V1=[ [1], [3], [1,3] ]
	// Noter qu'on ne met pas l'ensemble vide dans V1
	vector<int> v;
	if (i >f.size()){
		return ;
	}
	else if (i == f.size()){
		for(int j = 0; j < f.size(); ++j){
			if (f[j]!=0) {
				v.push_back(f[j]);
			}
		}
		if(v.size()>0) V1.push_back(v);
		return ;
	}
	int x=f[i];
	f[i] = 0;
	mis(f, i + 1, V1);
	f[i] = x;
	mis(f, i + 1, V1);
	return ;
}

void vect_to_int(vector<int> &V, int &i){
	//transforme un vecteur V de positions en un entier i
	//Ex: si V=[1, 3, 4] alors i = 1101 (en binair)
	i=0;
	for(int j=0; j<V.size();j++){
		if(V[j]==1){
			i++;
		}
		else{
			i+=(1<<(V[j]-1));
		}
	}
}

void vect_of_vect_to_vect_of_int (vector<vector<int> > &v1, vector<int> &v2){
	for(int i=0; i<v1.size();i++){
		int k;
		vect_to_int(v1[i],k);
		v2.push_back(k);
	}
}

void deleteEmptyEntries(vector<Space> &poidsDualSpace, vector<DataType> &indiceValide){
    DataType i;
    vector<DataType> result;
    for (i=0;i<(DataType)indiceValide.size();i++) if (poidsDualSpace[indiceValide[i]]>0) result.push_back(indiceValide[i]);
    indiceValide.swap(result);
}

void fusionGloutonne(USetDualSpace &uSetDualSpace, Space d){
	//fonction qui retourne un sous-ensemble minimal de USetDualSpace et qui lui est équivalent
	//le principe de la fonction est un algorithme glouton
    DataType i;
    USetDualSpace result;
    Space all=1<<d;

    map<DualSpace, vector<Space>> M;

    auto itu=uSetDualSpace.begin();
    //pour chaque paire, on génère les espaces qu'elle couvre 
    while (itu!=uSetDualSpace.end()){
        vector<Space> listCouv;
        listCouverts(itu->dom, itu->equ, d, listCouv);//les espaces couverts par la paire itu
        M.insert(pair<DualSpace, vector<Space>>(*itu, listCouv));
        ++itu;
    }

/* On crée un vecteur espacesCouverts de bool du nombre total d'espaces. 
espaceCouverts[i]=true ssi i est un espace couvert par au moins une paire */
/*
    bool* espacesCouverts= new bool[all];//quand d=30, ceci devient un bottleneck: 1G de booléens !!
    DataType nbCouverts=0;
    for (i=0;i<all;i++) espacesCouverts[i]=false;//on intialise
    for (auto it1=M.begin();it1!=M.end();++it1){//pour chaque paire
        for (auto j=0;j<(it1->second).size();++j){ //pour chaque espace couvert par la paire
            if(!espacesCouverts[(it1->second)[j]]){
                espacesCouverts[(it1->second)[j]]=true;
                nbCouverts++;//le nombre de couverts est incrémenté
            }
        }
    }
	delete [] espacesCouverts;

    

    if (nbCouverts==all-1){
        USetDualSpace aux;
        DualSpace ds;ds.dom=all-1;ds.equ=0;ds.poids=all-1;
        aux.insert(ds);
        uSetDualSpace.swap(aux);
       // delete [] espacesCouverts;
        return;
    }
*/

    DataType maxValue;
    map<DualSpace, vector<Space>>::iterator itX, itm, iti;
    DualSpace newDs, inter, ds;
    DataType nbDelete;
    vector<Space>::iterator itj, itk;

    while(M.size()>0){
        //on cherche la paire qui couvre le max d'espaces non encore couverts
        itX=M.begin();
        maxValue=(itX->second).size();
        itm=M.begin();
        while (itm!=M.end()){
            if ((Space)(itm->second).size()>maxValue){
                itX=itm;
                maxValue=(itm->second).size();
            }
            itm++;
        }
	   //on met les espaces couverts par la paire choisie dans un vecteur intermédiaire
        newDs = itX->first;
	   //on ajoute cette paire au résultat
        result.insert(newDs);
	   //on efface cette paire
        M.erase(itX);
        
        nbDelete=0;
        for (iti=M.begin(); iti!=M.end(); iti++){
            ds=iti->first;
		  //on construit une paire inter qui résume les espaces communs à newDS et iti->first
            //inter.dom = newDs.dom & ds.dom;
            //inter.equ = soustraction((newDs.dom+newDs.equ) & (ds.dom+ds.equ), inter.dom);
            vector<Space> usp;
            //chaque espace de iti->second qui n'est pas couvert par inter, on le met dans usp
            for (auto itk=(iti->second).begin();itk!=(iti->second).end();itk++){
                //if (!estCouvertPar(*itk,inter)) usp.push_back(*itk);
				if (!estCouvertPar(*itk,newDs)) usp.push_back(*itk);
            }
            //les espaces qui restent associés à iti sont ceux non couverts par inter, cad ceux dans usp
		  //on fait donc l'échange entre les deux
            (iti->second).swap(usp);
            if ((iti->second).size()==0) nbDelete++;
        }
        //on enlève de M toutes les paires dont le nombre d'espaces qui restent à couvrir est nul.
        if (nbDelete>0){
            iti=M.begin();
            while (iti!=M.end()){
                itX=iti;iti++;
                if ((itX->second).size()==0) M.erase(itX);
            }
        }
    }
    uSetDualSpace.swap(result);
}


/*resultat minimal Version 1*/
void fusionGloutonne2(USetDualSpace &uSetDualSpace, Space d){
    long i;
    USetDualSpace result;

    /*Construction de Dom(t)*/
    USetDualSpace::iterator itu=uSetDualSpace.begin();
    unordered_set<Space> domt;
    while (itu!=uSetDualSpace.end()){
        vector<Space> listCouv;
        listCouverts(itu->dom, itu->equ, d, listCouv);
        for (auto it=listCouv.begin();it!=listCouv.end();it++) domt.insert(*it);
        itu++;
    }

    vector<vector<Space>> vdomt(d+1);
    vector<vector<char>> vcovt(d+1);
    i=0;
    for (auto it=domt.begin();it!=domt.end();it++){
        vdomt[__builtin_popcount(*it)].push_back(*it);
        vcovt[__builtin_popcount(*it)].push_back(0);//non couvert
        i++;
    }


    long nb_a_couvrir=domt.size();
    for (auto dim=d; dim>=1 && nb_a_couvrir>0; dim--){
        for (i=0;i<vdomt[dim].size() && nb_a_couvrir>0;i++){
            if (vcovt[dim][i]==0){
                //vector<Space> lesFils;
			 int nbFils=0;
                Space inter=vdomt[dim][i];
                for (auto j=0;j<vdomt[dim-1].size();j++){
                    if (estInclusDans(vdomt[dim-1][j], vdomt[dim][i])){
                        inter=(inter & vdomt[dim-1][j]);
                        //lesFils.push_back(vdomt[dim-1][j]);
					nbFils++;
                    }
                }
                DualSpace ds;
                //if (lesFils.size()==dim){//si le nombre de fils où t est dominé est égal à dim
			 if(nbFils==dim){
                    ds.dom=vdomt[dim][i];
                    ds.equ=0;
				vector<Space> listc;
                    long nb=0;
                    listCouverts(ds.dom, ds.equ, d, listc);
                    Space Z=0;
                    for (auto it=listc.begin(); it!=listc.end(); it++) if (domt.find(*it)==domt.end()){
                        Z = Z | (*it);
                    }
                    ds.dom=ds.dom-Z;
                    ds.equ=Z;

                }else{
                    ds.dom=inter;
                    ds.equ=vdomt[dim][i]-inter;
                }

                result.insert(ds);

                /********************/

                vcovt[dim][i]=1;nb_a_couvrir--;

                for (auto dim1=dim-1; dim1>=1; dim1--){
                    for (auto j=0;j<vdomt[dim1].size();j++){
                        //cerr<<"je suis la\n";cerr<<" et dim="<<dim<<endl;
                        if (estCouvertPar(vdomt[dim1][j], ds)){
                            if(vcovt[dim1][j]==0){vcovt[dim1][j]=1;nb_a_couvrir--;}
                        }
                    }
                }

            }

        }
    }
    uSetDualSpace.swap(result);
  //  cerr<<"---------------------------------\n";
}

/*resultat minimal Version 2*/
void fusionGloutonne3(USetDualSpace &uSetDualSpace, Space d){
    long i;
    USetDualSpace result;
    Space all=1<<d;
    bool* domt=new bool[all];
    for (i=0;i<all;i++) domt[i]=false;

    /*Construction de Dom(t)*/
    auto  itu=uSetDualSpace.begin();
    while (itu!=uSetDualSpace.end()){
        vector<Space> listCouv;
        listCouverts(itu->dom, itu->equ, d, listCouv);
        for (auto it=listCouv.begin();it!=listCouv.end();it++){
            domt[*it]=true;
        }
        itu++;
    }
	/*for (auto i=1; i<all-1;i++){
		for(auto it=uSetDualSpace.begin(); it!=uSetDualSpace.end();it++){
			if((i&it->equ)!=i && ((i & (it->dom|it->equ))==i)){
				domt[i]=true;
				break;
			}
		}
	}
	*/
    vector<vector<Space>> vdomt(d+1);
    vector<vector<Space>> ndomt(d+1);
    vector<vector<char>> vcovt(d+1);
    Space sizeDomt=0;
    for (i=1;i<all;i++){
        if (domt[i]){
            vdomt[__builtin_popcount(i)].push_back(i);
            vcovt[__builtin_popcount(i)].push_back(0);//non couvert
            sizeDomt++;
        }
        else{
            ndomt[__builtin_popcount(i)].push_back(i);
        }
    }

    long nb_a_couvrir=sizeDomt;
    for (auto dim=d; dim>=1 && nb_a_couvrir>0; dim--){
        for (i=0;i<vdomt[dim].size() && nb_a_couvrir>0;i++){
            if (vcovt[dim][i]==0){
                DualSpace ds;
			//vector<Space> listc;
                ds.dom=vdomt[dim][i];
                ds.equ=0;

                Space ddim=dim-1;
                Space Y=0;
                while (ddim>0){
                    for (auto j=0;!Y && j<ndomt[ddim].size();j++) if (estInclusDans(ndomt[ddim][j], ds.dom)){
                        Y=ndomt[ddim][j];
                        ddim=0;
                    }
                    ddim--;
                }

                ds.dom-=Y;
                ds.equ=Y;

                result.insert(ds);

                vcovt[dim][i]=1;nb_a_couvrir--;

                for (auto dim1=dim-1; dim1>=1; dim1--){
                    for (auto j=0;j<vdomt[dim1].size();j++){
                        if (estCouvertPar(vdomt[dim1][j], ds)){
                            if(vcovt[dim1][j]==0){vcovt[dim1][j]=1;nb_a_couvrir--;}
                        }
                    }
                }

            }

        }
    }
    delete []domt;
    uSetDualSpace.swap(result);
  //  cerr<<"---------------------------------\n";
}

#endif // ALGOGLOUTON_H_INCLUDED
