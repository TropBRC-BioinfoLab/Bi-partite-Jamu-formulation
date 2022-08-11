#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>
#include <string>
#include <cctype>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <vector>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <list>
#include <bitset>
#include <functional>
#include <ctime>
#include <omp.h>

using namespace std;

#define mp make_pair
#define pb push_back
#define kandidat 1000
#define T 4 //thread
#define K 4 //kombinasi tanaman

typedef long long int ll;
typedef pair<double,double> dd;
typedef map<ll,pair<string,dd> > kamus; //gi : gene,(bobotProt, bobotEdge)
typedef pair<double,pair<vector<string>,kamus> > DT; // score, (namaLatin, kamus)

map<string,kamus > db; // namaLatin: kamus
priority_queue<DT> pq[T],pqAns;
map<string,kamus >::iterator it;
vector<pair<string,kamus> > mydata[T];
kamus value[T];
vector<string> nama[T];

ll counter = 0;
ll previ = 0;
ll alls;
ll part;

int cnt;

ll nChoosek( ll n, ll k ){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    ll result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (ll)(n-i+1);
        result /= (ll)i;
    }
    return result;
}

void rek(int depth, int loop){
	int id = omp_get_thread_num();
	//printf("thread: %d, depth: %d\n", id,depth);
	if(depth==K){
		double sc = 0.0;

		for(kamus::iterator itk = value[id].begin();itk!=value[id].end();++itk){
			double prot = itk->second.second.first;
			double edge = itk->second.second.second;
			//printf("%lld : %.10lf, %.10lf_", itk->first, prot,edge);
			//cout<<itk->second.first<<endl;

			sc += (prot*edge);
		}

		pq[id].push(mp(-sc,mp(nama[id],value[id])));

		int saiz = pq[id].size();

		if(saiz>kandidat){
			pq[id].pop();
		}

		if(id==0){
			counter++;
			if(counter/part>previ){
//				cout<<(double)100.0*counter/alls<<endl<<flush;
				previ = counter/part;
			}
		}
		//printf("thread: %d, counter: %lld\n", id,counter);
	}else if(depth==0){
		for(int i=loop;i<=cnt-(K-depth);i+=T){
			kamus previV = value[id];
			vector<string> previN = nama[id];

			string namaLatin = mydata[id][i].first;
			kamus value1 = mydata[id][i].second;

			for(kamus::iterator itk = value1.begin();itk!=value1.end();++itk){
				ll gi = itk->first;
				if(value[id].find(gi)!=value[id].end()){
					double sc1 = value[id][gi].second.second;
					double sc2 = value1[gi].second.second;
					if(sc2>sc1){
						value[id][gi].second.second = sc2;
					}
				}else{
					value[id][gi] = value1[gi];
				}
			}
			nama[id].pb(namaLatin);

			rek(depth+1,i+1);

			value[id] = previV;
			nama[id] = previN;
		}
	}else{
		for(int i=loop;i<=cnt-(K-depth);i++){
			kamus previV = value[id];
			vector<string> previN = nama[id];

			string namaLatin = mydata[id][i].first;
			kamus value1 = mydata[id][i].second;

			for(kamus::iterator itk = value1.begin();itk!=value1.end();++itk){
				ll gi = itk->first;
				if(value[id].find(gi)!=value[id].end()){
					double sc1 = value[id][gi].second.second;
					double sc2 = value1[gi].second.second;
					if(sc2>sc1){
						value[id][gi].second.second = sc2;
					}
				}else{
					value[id][gi] = value1[gi];
				}
			}
			nama[id].pb(namaLatin);

			rek(depth+1,i+1);

			value[id] = previV;
			nama[id] = previN;
		}
	}
}

void postKomb(){
	ll ctr = 0;
	for(int id = 0; id<T; id++){
		while(!pq[id].empty()){
			pqAns.push(pq[id].top());
			pq[id].pop();
			ctr++;
			if(ctr>kandidat){
				pqAns.pop();
			}
		}
	}
}

int main(){
	freopen("tanaman_protein.in","r",stdin);

	int nTan;
	scanf("%d ",&nTan);
	//getchar();
	while(nTan--){
		//char latin[100];
		//gets(latin);
		string latin;
		getline(cin,latin);
		//puts(latin);
		//getchar();
		int nProt;
		scanf("%d",&nProt);
		kamus tmp;
		for(int i=0;i<nProt;i++){
			string gen;
			double b,e;
			ll gi;
			//scanf("%lld %s %lf %lf ",&gi,&gen,&b,&e);
			cin>>gi>>gen>>b>>e;
			getchar();
			//printf("%lld %s %.10lf %.10lf\n", gi,gen,b,e);

			tmp[gi] = mp(gen,mp(b,e));
		}
		db[latin] = tmp;
		//printf("nTan : %d\n", nTan);
	}


	for(it= db.begin();it!=db.end();++it){
		for(int id = 0;id<T;id++)
			mydata[id].push_back(mp(it->first,it->second));
	}

	cnt = mydata[0].size();
	alls = nChoosek(cnt,K)/T;
	part = alls/10000;

	double waktu = omp_get_wtime();

	#pragma omp parallel num_threads(T)
	rek(0,omp_get_thread_num());

	waktu = omp_get_wtime() - waktu;
	printf("waktu: %.10lf\n", waktu);

	postKomb();

	while(!pqAns.empty()){
		DT atas = pqAns.top();
		pqAns.pop();

		pair<vector<string>,kamus> out = atas.second;
		cout<<-atas.first<<" : "<<out.second.size()<<", [";
		for(int i=0;i<out.first.size();i++){
			cout<<out.first[i]<<", ";
		}
		cout<<"], {";

		for(kamus::iterator itk=out.second.begin();itk!=out.second.end();++itk){
			cout<<"("<<itk->first<<" : "<<itk->second.first<<", "<<itk->second.second.first<<", "<<itk->second.second.second<<"),";
		}
		cout<<"}\n";
		break ;
	}

	return 0;
}
