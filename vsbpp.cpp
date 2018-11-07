#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <float.h>

#include <vector>
#include <bits/stdc++.h>
#include <random>
using namespace std;

#define DEBUG 0


typedef struct move_info
{
    double value;
    int i;
    int j;
} move_tabu;


/**********************************************************************
 *                                                                    *
 * SUBSET SUM USING PERFECT SUM ALGORITHM (Dynamic Programming)       *
 *                                                                    *
 **********************************************************************/
// dp[i][j] is going to store true if sum j is
// possible with array elements from 0 to i.
bool** dp;
 
void display(const vector<int>& v)
{
    for (int i = 0; i < v.size(); ++i)
        printf("%d ", v[i]);
    printf("\n");
}
 
// A recursive function to print all subsets with the
// help of dp[][]. Vector p[] stores current subset.
void printSubsetsRec(vector<int> &  arr, int i, int sum, vector<int>& p, vector< vector<int> > & items_in_bins)
{
    // If we reached end and sum is non-zero. We print
    // p[] only if arr[0] is equal to sun OR dp[0][sum]
    // is true.
    if (i == 0 && sum != 0 && dp[0][sum])
    {
        p.push_back(arr[i]);
        //display(p);
        items_in_bins.push_back(p);
        return;
    }
 
    // If sum becomes 0
    if (i == 0 && sum == 0)
    {
        //display(p);
        items_in_bins.push_back(p);
        return;
    }
 
    // If given sum can be achieved after ignoring
    // current element.
    if (dp[i-1][sum])
    {
        // Create a new vector to store path
        vector<int> b = p;
        printSubsetsRec(arr, i-1, sum, b, items_in_bins);
    }
 
    // If given sum can be achieved after considering
    // current element.
    if (sum >= arr[i] && dp[i-1][sum-arr[i]])
    {
        p.push_back(arr[i]);
        printSubsetsRec(arr, i-1, sum-arr[i], p, items_in_bins);
    }
}
 
// Prints all subsets of arr[0..n-1] with sum 0.
void printAllSubsets(vector<int> & arr, int n, int sum, vector< vector<int> > & items_in_bins)
{
    if (n == 0 || sum < 0)
       return;
 
    // Sum 0 can always be achieved with 0 elements
    dp = new bool*[n];
    for (int i=0; i<n; ++i)
    {
        dp[i] = new bool[sum + 1];
        dp[i][0] = true;
    }
 
    // Sum arr[0] can be achieved with single element
    if (arr[0] <= sum)
       dp[0][arr[0]] = true;
 
    // Fill rest of the entries in dp[][]
    for (int i = 1; i < n; ++i)
        for (int j = 0; j < sum + 1; ++j)
            dp[i][j] = (arr[i] <= j) ? dp[i-1][j] ||
                                       dp[i-1][j-arr[i]]
                                     : dp[i - 1][j];
    if (dp[n-1][sum] == false)
    {
        printf("There are no subsets with sum %d\n", sum);
        return;
    }
 
    // Now recursively traverse dp[][] to find all
    // paths from dp[n-1][sum]
    vector<int> p;
    printSubsetsRec(arr, n-1, sum, p, items_in_bins);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/**********************************************************************
 *                                                                    *
 *                MEMORY ALLOCATION FUNCTIONS                         *
 *                                                                    *
 **********************************************************************/
void nrerror(char *error_text)
{
    fprintf(stderr,"\n%s\n",error_text);
    fprintf(stderr,"... now exiting to system ...\n");
    exit(1);
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    int **m;
    m = (int **)malloc((unsigned)(nrh - nrl +1)*sizeof(int*));
    if (!m) nrerror("allocation failure 1 in imatrix()");
    m -= nrl;
    for(i = nrl; i <= nrh ; i++)
    {
        m[i]= (int*)malloc((unsigned)(nch - ncl +1)*sizeof(int));
        if (!m[i])
            nrerror("allocation failure 2 in imatrix()");
        m[i] -= ncl;
    }
    return m;
}
void free_imatrix(int **m, int nrl, int nrh, int ncl)
{
    int i;
    for (i = nrh; i >= nrl; i--)
        free((char*) (m[i] + ncl));
    free((char *)(m + nrl));
}
/******************************************************************************/
void read_instance(
    char *arquivo, 
    int bin_types, 
    vector<int> & binCap, 
    vector<int> & binCost,
    vector<int> & weight, 
    int items
    )
{
    FILE *fp;
    int i = 0;
    int j = 0;
    int tmp = 0;
    int *d;
    
    d = (int *) malloc(sizeof(int) * bin_types*2);
    
    fp = fopen(arquivo, "r");
    if(fp == NULL)
        exit(1);
    
    for(i=0; i<4; i++)
        fscanf(fp, "%d", &tmp);
    
    for(j = 0; j < bin_types*2; j++)
    {
        fscanf(fp, "%d", &tmp);
        d[j] = tmp;
    }
    
    i = 0;
    for(j = 0; j < bin_types; j++)
    {
        binCap[j] = d[i];
        binCost[j] = d[i+1];
        i+=2;
    }
    
    for(i=0; i<items; i++)
    {
        fscanf(fp, "%d", &tmp);
        weight[i] = tmp;
    }
    
    free(d);
    fclose(fp);
}
/******************************************************************************/
void read_prob_sizes(
    char * arquivo, 
    int *bin_types, 
    int *items
    )
{
    FILE *fp;
    fp = fopen(arquivo, "r");
    if(fp == NULL)
        exit(1);
    
    fscanf(fp, "%d", items);
    fscanf(fp, "%d", bin_types);
    fclose(fp);
}
/******************************************************************************/
double calcula_penalidade(
    vector<int>  & individuo, 
    vector<int>  & binCap, 
    vector<int>  & binCost,
    vector<int>  & weight, 
    int items, 
    int bin_types,
    int *etiqueta_pred, 
    int *etiqueta_custo
    )
{
    double penalidade = 0.0;
    double tmp = 0.0;
    double ociosidade = 0.0;
    double ociosidade_total = 0.0;
    int peso = 0;
    int peso_excesso = 0;
    int peso_excesso_min = INT_MAX;
    int i, j;

    i = items;
    do
    {
#if DEBUG
        printf("noh %d => predecessor %d; novo valor de i = %d; custo = %d\n", i, etiqueta_pred[i], etiqueta_pred[i], etiqueta_custo[i]);
#endif
        peso = 0;
        for(j = i-1; j >= etiqueta_pred[i]; j--)
        {
            if(j < items)
            {
#if DEBUG
                printf("%d) individuo[%d] = %d -  Peso = %d\n", j, j, individuo[j], weight[individuo[j]]);
                printf("%d;%d\n", individuo[j]+1, weight[individuo[j]]);
#endif                
                peso += weight[individuo[j]];
            }
        }
#if DEBUG
        printf("Peso acumulado;%d\n\n", peso);
#endif      
        if(peso <= binCap[0])
        {
            ociosidade = binCap[0] - peso;
        }
        else if(peso <= binCap[1])
        {
            ociosidade = binCap[1] - peso;
        }
        else if(peso <= binCap[2])
        {
            ociosidade = binCap[2] - peso;
        }

        ociosidade_total += ociosidade;
        i = etiqueta_pred[i];

    } while(i > 0);

#if DEBUG
    printf("ociosidade total = %.2f\n", ociosidade_total);
#endif    
    penalidade = ociosidade_total * ociosidade_total;
#if DEBUG    
    printf("Penalidade = %.2f\n", penalidade);
#endif
    return penalidade;
}

/******************************************************************************/

double calc_fitness(
    vector<int> & individuo, 
    vector<int>  & binCap, 
    vector<int>  & binCost,
    vector<int>  & weight, 
    int items, 
    int bin_types, 
    vector<int>  & bin_item, 
    int penalize
    )
{
    double f = 0.0;
    int *Q; /* Peso acumulado */
    int *etiqueta_custo;
    int *etiqueta_pred;
    int i, j;
    int custo_acumulado = 0;
    int peso = 0;
    double penalidade = 0.0;
    double fi_bi = DBL_MAX; /* minima razao entre fi / bi (custo / capacidade) */
    double tmp = 0.0;
    int peso_excesso = 0;
    int peso_excesso_min = INT_MAX;

    int binNumber = 1;
    
    Q = (int *) malloc(sizeof(int) * (items + 1));
    etiqueta_custo = (int *) malloc(sizeof(int) * (items + 1));
    etiqueta_pred = (int *) malloc(sizeof(int) * (items + 1));
    
    for(i = 0; i < items + 1; i++)
    {
        etiqueta_custo[i] = 0;
        etiqueta_pred[i] = -1;
    }
    
    Q[0] = 0;
    for(i = 1; i < items + 1; i++)
        Q[i] = Q[i - 1] + weight[individuo[i - 1]];

    for(i = 0; i < items; i++)
    {
        for(j = i + 1; j < items+1; j++)
        {
            int val = Q[j] - Q[i];
            int entrou = 0;
            custo_acumulado = etiqueta_custo[i];
            if(val <= binCap[0])
            {
                custo_acumulado += binCost[0];
                entrou = 1;
            }
            else if(val <= binCap[1])
            {
                custo_acumulado += binCost[1];
                entrou = 1;
            }
            else if(val <= binCap[2])
            {
                custo_acumulado += binCost[2];
                entrou = 1;
            }
            
            if(entrou)
            {
                if(etiqueta_pred[j] == -1 || custo_acumulado < etiqueta_custo[j])
                {
                    etiqueta_custo[j] = custo_acumulado /*+ etiqueta_custo[i]*/;
                    etiqueta_pred[j] = i;
                }
            } else break;
        }
    }

    f = etiqueta_custo[items];

    if(penalize)
    {
        penalidade = calcula_penalidade(individuo, binCap, binCost,
            weight, items, bin_types, 
            etiqueta_pred, etiqueta_custo);
#if DEBUG       
        printf("Penalidade = %.2f\n", penalidade);
#endif
        f += penalidade;
    }
    
    i = items;
    bin_item[i] = binNumber;
    do
    {
        for(j = i-1; j >= etiqueta_pred[i]; j--)
            if(j < items)
                bin_item[j] = binNumber;

        i = etiqueta_pred[i];
        binNumber++;
    } while(i > 0);

#if DEBUG
    if(penalize)
    {
        for(i=0; i < items+1; i++)
            printf("i = %d, individuo[%d] = %d, Peso[%d] = %d, bin_item[%d] = %d\n", i, i, individuo[i], individuo[i], weight[individuo[i]], i, bin_item[i]);
        printf("Custo = %.2f\n\n", f);
    }

    for(i=0; i < items+1; i++)
        printf("i = %d, individuo[%d] = %d, Peso[%d] = %d, bin_item[%d] = %d\n", i, i, individuo[i], individuo[i], weight[individuo[i]], i, bin_item[i]);
    printf("Custo = %.2f\n\n", f);
#endif
    
    free(Q);
    free(etiqueta_custo);
    free(etiqueta_pred);
    Q = NULL;
    etiqueta_custo = NULL;
    etiqueta_pred = NULL;
    return f;
}

/******************************************************************************/
void best_move(move_tabu *mov, int iter, double c_best, double c_curr,
               int **tabu_time, vector<int> & individuo, int items, vector<int>  & binCap,
               vector<int>  & binCost, vector<int>  & weight, int bin_types, vector<int>  & bin_item, int penalize)
{
    int i, j, k;
    k = 0;
    mov->value = (double) INT_MAX;
    for(i = 0; i < items; i++)
    {
        for(j = 0; j < items; j++)
        {
            if(bin_item[i] != bin_item[j] && weight[individuo[i]] != weight[individuo[j]])
            {
                int tmp = individuo[i];
                individuo[i] = individuo[j];
                individuo[j] = tmp;
                c_curr = calc_fitness(individuo, binCap, binCost, weight, items, bin_types, bin_item, penalize);
                if(tabu_time[i][j] < iter || c_curr < c_best)
                {
                    if(c_curr < mov->value)
                    {
                        mov->value = c_curr;
                        mov->i = i;
                        mov->j = j;
#if DEBUG
                        printf("Entrou: mov->value = %.2f\n", mov->value);
#endif
                    }
                }
                tmp = individuo[j];
                individuo[j] = individuo[i];
                individuo[i] = tmp;
            }
        }
    }
}

void execute_move(move_tabu *mov, vector<int>  & individuo)
{
    int tmp = individuo[mov->i];
    individuo[mov->i] = individuo[mov->j];
    individuo[mov->j] = tmp;
}

double tabu_search(vector<int>  & individuo, vector<int>  & binCap, vector<int>  & binCost,
                   vector<int>  & weight, int items, int bin_types, double cost, vector<int>  & bin_item,
                   int tabu_size, int max_iter, int penalize)
{
    double ts_cost = 0.0;
    int i, j;
    int iter;
    int **tabu_time;
    vector<int> tmp_sol;
    vector<int> best_sol;
    double c_best, c_curr;
    move_tabu mov;
    tabu_time = imatrix(0, items, 0, items);
    tmp_sol.resize(items+1);
    best_sol.resize(items+1);

    for(i = 0; i < items; i++)
    {
        tmp_sol[i] = individuo[i];
        best_sol[i] = individuo[i];
        for(j = 0; j < items; j++)
            tabu_time[i][j] = 0;
    }
    tmp_sol[items] = individuo[items];
    best_sol[items]= individuo[items];
    
    iter = 0;
    c_curr = cost;
    c_best = c_curr;
    
    while(iter < max_iter)
    {
        iter++;
        best_move(&mov, iter, c_best, c_curr, tabu_time, tmp_sol, items, binCap, binCost, weight, bin_types, bin_item, penalize);
        execute_move(&mov, tmp_sol);
        tabu_time[mov.i][mov.j] = iter + tabu_size;
        c_curr = mov.value;
        if(c_curr < c_best)
        {
            c_best = c_curr;
            for(i = 0; i < items; i++) individuo[i] = tmp_sol[i];
#if DEBUG           
            printf("Iter: %d) ", iter);
            printf("Cost = %.2f\n", c_best);
            printf("\n");
#endif          
        }
    }
    ts_cost = c_best;
#if DEBUG   
    if(penalize)
        printf("Final TS (penalized) cost = %.2f\n", ts_cost);
#endif

    ts_cost = calc_fitness(individuo, binCap, binCost, weight, items, bin_types, bin_item, 0);
    printf("Final TS (no penalized) cost = %.2f\n", ts_cost);
    free_imatrix(tabu_time, 0, items, 0);
    
    return ts_cost;
}

/******************************************************************************/
void shuffle(vector<int> &array, size_t n, int seed)
{
    srand(seed);
    //srand(time(NULL)); /* descomente se preferir uma semente aleatoria a cada execucao */
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}
/******************************************************************************/
void sort(vector<int>  & weight, vector<int>  & individuo, int items)
{
    int i, j, temp;
    for(i = 0; i < items; i++)
    {
        for(j = i + 1; j < items; j++)
        {
            if(weight[individuo[i]] < weight[individuo[j]])
            {
                temp = individuo[i];
                individuo[i] = individuo[j];
                individuo[j] = temp;
            }
        }
    }
}
/******************************************************************************/
int main(int argc, const char * argv[])
{
    clock_t start;
    double tempo = 0.0;
    char arquivo[80];
    int bin_types = 0;
    int items = 0;
    vector<int> binCap;
    vector<int> binCost;
    vector<int> weight;
    vector<int> individuo;
    vector<int> bin_item;
    double myFitness = 0.0;
    int i;
    int iseed; /* semente para geracao de numeros aleatorios */
    int tabu_size;
    int max_iter;
    int bin_for_subsetsum;
    
    if(argc == 6)
    {
        sscanf(argv[1], "%s", arquivo);
        tabu_size = atoi(argv[2]);
        max_iter = atoi(argv[3]);
        iseed = atoi(argv[4]);
        bin_for_subsetsum = atoi(argv[5]);
    } else {
        printf("Wrong number of arguments: a.out <instance> <tabu_size> <max_iter> <seed> <bin size for subset sum>\n");
        exit(1);
    }
    
    read_prob_sizes(arquivo, &bin_types, &items);
    binCap.resize(bin_types);
    binCost.resize(bin_types);
    weight.resize(items);
    individuo.resize(items+1);
    bin_item.resize(items+1);

    read_instance(arquivo, bin_types, binCap, binCost, weight, items);
    
#if DEBUG
    for(i = 0; i < items; i++)
        printf("%d) item %d - peso %d\n", i+1, individuo[i], weight[individuo[i]]);
    printf("\n");
#endif

    start = clock();
    for(i = 0; i < items; i++)
        individuo[i] = i;
    
    sort(weight, individuo, items);
    shuffle(individuo, items, iseed);
    /*
    // Solucao otima! Basta descomentar esse trecho de codigo.
    individuo[0] = 0;
    individuo[1] = 31;
    individuo[2] = 1;
    individuo[3] = 8;
    individuo[4] = 2;
    individuo[5] = 36;
    individuo[6] = 39;
    individuo[7] = 3;
    individuo[8] = 45;
    individuo[9] = 4;
    individuo[10] = 35;
    individuo[11] = 48;
    individuo[12] = 5;
    individuo[13] = 20;
    individuo[14] = 6;
    individuo[15] = 33;
    individuo[16] = 37;
    individuo[17] = 7;
    individuo[18] = 25;
    individuo[19] = 9;
    individuo[20] = 15;
    individuo[21] = 21;
    individuo[22] = 10;
    individuo[23] = 12;
    individuo[24] = 41;
    individuo[25] = 11;
    individuo[26] = 19;
    individuo[27] = 23;
    individuo[28] = 13;
    individuo[29] = 38;
    individuo[30] = 14;
    individuo[31] = 17;
    individuo[32] = 46;
    individuo[33] = 16;
    individuo[34] = 29;
    individuo[35] = 40;
    individuo[36] = 18;
    individuo[37] = 44;
    individuo[38] = 22;
    individuo[39] = 26;
    individuo[40] = 27;
    individuo[41] = 47;
    individuo[42] = 24;
    individuo[43] = 28;
    individuo[44] = 32;
    individuo[45] = 30;
    individuo[46] = 34;
    individuo[47] = 42;
    individuo[48] = 43;
    individuo[49] = 49;
    */

    // Guarda todos os pesos em um multimap tendo a chave o peso e o valor o número do item.
    multimap <int, int> mymap;
    for(int i = 0; i < items; i++)
        mymap.insert(pair <int, int> (weight[i], i));


    // Faz o subset-sum usando o perfect sum (dynamic programming)
    vector< vector<int> > items_in_bins;
    int val = 0;

    printAllSubsets(weight, items, bin_for_subsetsum, items_in_bins);

    // limpa o vetor individuo
    individuo.clear();

    // insere uma sequencia aleatoria do subsetsum no individuo.
    int qtd_combinacoes = items_in_bins.size();

    int min = 0;
    int max = qtd_combinacoes;
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(min,max); // guaranteed unbiased

    int random_integer = uni(rng);

    cout << "random_integer = " << random_integer << endl;
    for(int j = 0; j < items_in_bins[random_integer].size(); j++)
    {
        cout << items_in_bins[random_integer][j] << " ";
        val = items_in_bins[random_integer][j];
        // encontra o item com o referido peso.
        multimap<int,int>::iterator it = mymap.lower_bound(val);
        individuo.push_back(it->second);
        mymap.erase(mymap.lower_bound(val));
    }
    cout << endl;

    cout << "individuo\n";
    for(vector<int>::iterator it = individuo.begin(); it != individuo.end(); it++)
    {
        cout << "item = " << *it << " - peso = " << weight[*it] << endl;
    }

    /*
    for(int i = 0; i < items_in_bins.size(); i++)
    {
        cout << "i = " << i << " => ";
        val = 0;
        for(int j = 0; j < items_in_bins[i].size(); j++)
        {
            cout << items_in_bins[i][j] << " ";
            val += items_in_bins[i][j];
        }
        cout << " = " << val << endl;
    }
    */






    /*
    myFitness = calc_fitness(individuo, binCap, binCost, weight, items, bin_types, bin_item, 0);
    printf("Custo inicial sem penalizacao = %.2f\n", myFitness);

    myFitness = calc_fitness(individuo, binCap, binCost, weight, items, bin_types, bin_item, 1);
    printf("Custo inicial com a penalizacao = %.2f\n", myFitness);
    */
    /* Executa a busca tabu. O ultimo parametro = 1 se leva em conta a penalizacao, ou 0 c.c. */
    /*
    myFitness = tabu_search(individuo, binCap, binCost, weight, items, bin_types, myFitness, bin_item, tabu_size, max_iter, 1);
    printf("Custo tabu search = %.2f\n", myFitness);
    
    tempo = (double)(clock() - start)/CLOCKS_PER_SEC;
    
    printf("%s;%.2f;%.5f\n", arquivo, myFitness, tempo);
    */
    return 0;
}
