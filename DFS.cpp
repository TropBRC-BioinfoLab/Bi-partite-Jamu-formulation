#include <bits/stdc++.h>
using namespace std;
#define nl "\n"
#define ios ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL)

struct Item 
{ 
    int weight;
    unordered_map <string, double> value ;
    string name ;
    double totalValue;
}; 
  

struct Node 
{
    unordered_map <string, int> idxTan ;
    int level;
    double profit, bound;
    int weight; 
};
  
bool cmp(Item a, Item b) 
{ 
    double r1 = (double)a.totalValue / a.weight; 
    double r2 = (double)b.totalValue / b.weight; 
    return r1 > r2;
}
long long sumz ;
string proteinName[14] = {"Akt1", "Ep300", "Foxo1", "Gcgr", "Ins2", "Insr", "Kcnj11", "Mtnr1b", "Ppara", "Pparg", "Prkaca", "Sod3", "Stat3", "Tcf7l2"};
unordered_map <string,double> protein ;
Item arr[462] ;
set<int> sTan ;

  
double bound(Node u, int n, int W) 
{
    if (u.weight >= W)
        return 0;
  
    double profit_bound = u.profit;   
    int j = u.level + 1;
    int totweight = u.weight;
  
    while ((j < n) && (totweight + arr[j].weight <= W)) 
    { 
        totweight    += arr[j].weight; 
        profit_bound += arr[j].totalValue; 
        j++; 
    } 
  
    return profit_bound; 
}
  
double knapsack(int W, int n) 
{ 
    stack <Node> Q;
    Node u, v; 
    double temp ; 
    // create root node 
    u.level = -1; 
    u.profit = 0.0 ;
    u.weight = 0;
    for(int i=0 ; i < 14 ; i++)
    {
      u.idxTan[proteinName[i]] = -1 ;
    }
    Q.push(u); 
  
    double maxProfit = 0.0; 
    sumz = 0 ;
    while (!Q.empty())
    { 
        u = Q.top(); 
        Q.pop();
  
        if (u.level == -1) 
            v.level = 0; 
  
        if (u.level == n-1)
          continue;

        // generate first child
        v.level = u.level + 1;
        v.weight = u.weight; 
        temp = 0.0 ;
        for(int k = 0 ; k < 14 ;k++)
        {
          v.idxTan[proteinName[k]] = u.idxTan[proteinName[k]];
        }
        v.profit = u.profit ;
        v.bound = bound(v, n, W);
        if (v.bound > maxProfit) 
            Q.push(v);

        // generate second child
        v.weight = u.weight + arr[v.level].weight;
        temp = 0.0 ;
        for(int k = 0 ; k < 14 ;k++)
        {
          if(u.idxTan[proteinName[k]] == -1)
          {
              if(arr[v.level].value[proteinName[k]] > 0.0)
                v.idxTan[proteinName[k]] = v.level ;
              else
                v.idxTan[proteinName[k]] = -1 ;
          }
          else
          {
            if(arr[v.level].value[proteinName[k]] > arr[u.idxTan[proteinName[k]]].value[proteinName[k]])
                v.idxTan[proteinName[k]] = v.level ;
            else
                v.idxTan[proteinName[k]] = u.idxTan[proteinName[k]] ;
          }
          if(v.idxTan[proteinName[k]] >= 0.0)
            temp += (arr[v.idxTan[proteinName[k]]].value[proteinName[k]]*protein[proteinName[k]]) ;
        }
        v.profit = temp ;

        // change best solution  
        if (v.weight <= W && v.profit > maxProfit) 
        {
          maxProfit = v.profit;
          sTan.clear() ;
          for(int k = 0 ; k < 14 ; k++)
          {
            sTan.insert(v.idxTan[proteinName[k]]);
          }
        }
        v.bound = bound(v, n, W); 
        if (v.bound > maxProfit) 
            Q.push(v);

        // sum total of leaf node
        if(v.weight == W)
          sumz++; 
    } 
  
    return maxProfit; 
} 
  
int main() 
{
    clock_t tStart = clock();
    int W = 3 ;
    int nTan = 460 ;
    cin >> nTan ;
    for(int k = 0 ; k < nTan ; k++)
    {
      string latin ;
      string gen ;
      double bobot,edge ;
      long long gi ;
      int nProt ;
      arr[k].totalValue = 0.0 ;
      
      // scanning data
      getchar();
      getline(cin,latin);
      scanf("%d",&nProt);
      for(int i=0;i<nProt;i++){
        cin >> gi >> gen >> bobot >> edge ;
        arr[k].totalValue += (bobot*edge) ;
        arr[k].value[gen] = edge ;
        protein[gen] = bobot ;
      }
      for(int i=0 ; i < 14 ; i++)
      {
        if(!arr[k].value[proteinName[i]])
        {
          arr[k].value[proteinName[i]] = 0.0 ; 
        }
      }
      arr[k].name = latin ;
      arr[k].weight = 1 ;
    }
    // sorting
    sort(arr, arr + nTan, cmp);

    // BnB function
    cout << "Maximum possible profit = " << knapsack(W, nTan)<< endl ; 
    printf("Time taken: %.16fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    cout << "Possible solution :"<< sumz << " leaf node"<< endl ;
    for(auto i : sTan)
    {
      if(i >= 0)
        cout << "\t- " << arr[i].name << endl ;
    }    return 0; 
}