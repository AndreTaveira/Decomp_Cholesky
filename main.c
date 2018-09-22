#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//                         Lucas garavaglia cova
//                         NUSP 8956449

// Declaração de funções
double** coef(double **A, int tam);
double** alocarMatriz(int Linhas,int Colunas);
double* ResolveSis(double** L, int tam, double* b);
double* esper(double** mR,int tam, int lin);
double** matrR(double** P,int tam, int lin);
double* varian(double** mR, double* U,int tam, int lin);
double** Covarian(double** mR,double* U,int Natv, int lin);
double* multMatr(double** m1, double* m2, int tam);
double multVec(double* m1, double* m2, int tam);
double* mult_esc(double* vet, double a, int tam );
double* som_vec(double* vet1, double* vet2, int tam );

int main()
{

    FILE *f;

    char str[80]; // usado apenas para fazer os prints dos strings no txt
    int Natv; // armazena o numero de ativos
    int i; int j; int k; // variaveis auxiliares
    double dado; // variavel destinada a manusear os dados numericos do txt
    f = fopen("PRECOS.txt", "r");
    fscanf(f, "%d/n", &Natv); printf("%d\n",Natv); // Armazena o numero de ativos

    double **P = (double**) alocarMatriz(100,Natv); // aloca memória para a matriz de preços dos ativos

    for(i = 0; i < 56; i++) // não consegui lidar genericamente com os nomes dos ativos por isso fiz este loop especificamente para o txt PRECOS
    {
        fscanf(f, "%s", &str);
        printf("%s",str);
    }
    printf("\n");
    for(i = 0; i < 100; i++) // não consegui lidar genericamente com a quantidade de linhas da tabela de PRECOS então tive de fazer um for para setar um fim manualmente
    {
        fscanf(f, "%s", &str); // le a data
        printf("%s ",str); // imprime a data
        for(j = 0; j < Natv; j++) // le, printa e armazena em uma matriz de precos os dados do TXT
        {
            fscanf(f, "%lf", &dado);
            printf("%.4lf ", dado);
            P[i][j] = dado;
        }
        printf("\n");
    }

    double** mR = matrR( P, Natv, 100); // gera matriz com os Rij
    double* U = esper( mR, Natv, 100); // gera vetor R que armazena as esperanças
    double* Sigma = varian( mR, U, Natv, 100); // gera vetor de variancias
    double** Q = Covarian( mR, U, Natv, 100); // gera matriz de covariancia
    double** L = coef( Q, Natv); // gera matriz L relativa a matriz de covariancia Q

    double *u = (double*)malloc((Natv) * sizeof(double));
    u = ResolveSis(L, Natv, U); // Resolve sistema linear Q*u = R

    double *v = (double*)malloc((Natv) * sizeof(double));
    double *uns = (double*)malloc((Natv) * sizeof(double));
    for(k=0; k<Natv; k++) { uns[k] = 1; } // Cria vetor de uns
    v = ResolveSis(L, Natv, uns); // Resolve sistema linear Q*v = 1

    double e[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double **Qinv = (double**) alocarMatriz(Natv,Natv); // aloca memória para a matriz inversa
    for(k=0; k<Natv; k++)
    {
        e[k] = 1;
        Qinv[k] = ResolveSis(L, Natv, e);
        e[k] = 0;
    }

    double um[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    double *temp = (double*)malloc((Natv) * sizeof(double));
    temp = multMatr(Qinv,um,Natv);
    double A = multVec(um,temp,Natv);

    double *temp2 = (double*)malloc((Natv) * sizeof(double));
    temp2 = multMatr(Qinv,U,Natv);
    double B = multVec(um,temp2,Natv);

    double *temp3 = (double*)malloc((Natv) * sizeof(double));
    temp3 = multMatr(Qinv,U,Natv);
    double C = multVec(U,temp2,Natv);

    double max = U[0];
    for(k=0; k<Natv; k++)
    {
        if (U[k] > max)
        {
            max = U[k];
        }
    }
    // derivando a equação 7 e igualando a 0 temos o min de mi = B/A
    double min = B/A;
    double media = (max+min)/2;

    double lambda_min = (C - min*B)/(A*C-pow(B,2));
    double lambda_max = (C - max*B)/((A*C)-pow(B,2));
    double lambda_med = (C - media*B)/((A*C)-pow(B,2));

    double gama_min = (A*min-B)/(A*C-pow(B,2));
    double gama_max = (A*max-B)/(A*C-pow(B,2));
    double gama_med = (A*media-B)/(A*C-pow(B,2));

    double *v1 = (double*)malloc((Natv) * sizeof(double));
    double *v2 = (double*)malloc((Natv) * sizeof(double));
    double *w_min = (double*)malloc((Natv) * sizeof(double));
    v1 = mult_esc(v, lambda_min, Natv);
    v2 = mult_esc(u, gama_min, Natv);
    w_min = som_vec( v1, v2, Natv);

    double *v11 = (double*)malloc((Natv) * sizeof(double));
    double *v22 = (double*)malloc((Natv) * sizeof(double));
    double *w_max = (double*)malloc((Natv) * sizeof(double));
    v11 = mult_esc(v, lambda_max, Natv);
    v22 = mult_esc(u, gama_max, Natv);
    w_max = som_vec( v11, v22, Natv);

    double *v111 = (double*)malloc((Natv) * sizeof(double));
    double *v222 = (double*)malloc((Natv) * sizeof(double));
    double *w_med = (double*)malloc((Natv) * sizeof(double));
    v111 = mult_esc(v, lambda_med, Natv);
    v222 = mult_esc(u, gama_med, Natv);
    w_med = som_vec( v111, v222, Natv);


    // Faremos agora oa prints dos dados de interesse do problema. Cada passo será impresso no terminal

    // impressão de matriz de Retornos ( Contém os Rij que usados nos calculos )
    printf("------  Matriz de Retornos: \n");
    for(k=0; k<Natv-1; k++)
    {
        for(j=0; j<Natv; j++)
        {
            printf("%f ",mR[k][j]);
        }
        printf("\n");
    }

    // impressão de vetor de esperanças ( Vetor R )
    printf("------  (vetor R) Esperança: \n");
    for(k=0; k<Natv; k++)
    {
        printf("%f ",U[k]);
        printf("\n");
    }

    // impressão de vetor de Variancias ( Note que estes valores correspondem a diagonal da matriz de covariancia )
    printf("------  Variancia: \n");
    for(k=0; k<Natv; k++)
    {
        printf("%f ",Sigma[k]);
        printf("\n");
    }

    // impressão de matriz de covariancia
    printf("------ Covariancia: \n");
    for(k=0; k<Natv; k++)
    {
        for(j=0; j<Natv; j++)
        {
            printf("%f ",Q[k][j]);
        }
        printf("\n");
    }

    // impressão de matriz da decomposição L relativa a Q (matriz de covariancia)
    printf("------ L: \n");
    for(k=0; k<Natv; k++)
    {
        for(j=0; j<Natv; j++)
        {
            printf("%f ",L[k][j]);
        }
        printf("\n");
    }

    // impressão da solução do sistema Q*u = R
    printf("------  solução de u: \n");
    for(k=0; k<Natv; k++)
    {
        printf("%f ",u[k]);
        printf("\n");
    }

    // impressão da solução do sistema Q*v = 1
    printf("------  solução de v: \n");
    for(k=0; k<Natv; k++)
    {
        printf("%f ",v[k]);
        printf("\n");
    }

    // impressão de Qinv
    printf("------  Qinv: \n");
    for(k=0; k<Natv; k++)
    {
        for(j=0; j<Natv; j++)
        {
            printf("%f ",Qinv[k][j]);
        }
        printf("\n");
    }

    // impressão de A
    printf("\n------  A: \n A = %f \n",A);

    // impressão de B
    printf("\n------  B: \n B = %f \n",B);

    // impressão de C
    printf("\n------  C: \n C = %f \n\n",C);

    // portifolio min
    printf("------  portifolio min (w calculado com mi minimo): \n");
    for(k=0; k<Natv; k++)
    {
        printf("%f ",w_min[k]);
    }
    printf("\n");

    // portifolio max
    printf("------  portifolio max (w calculado com mi maximo): \n");
    for(k=0; k<Natv; k++)
    {
        printf("%f ",w_max[k]);
    }
    printf("\n");

    // portifolio media
    printf("------  portifolio media (w calculado com mi na media entre max e min): \n");
    for(k=0; k<Natv; k++)
    {
        printf("%f ",w_med[k]);
    }
    printf("\n");



/*
    int k;
    double a[4] = {1,2,5,3};
    double *w = (double*)malloc((4) * sizeof(double));
    w = mult_esc(a, 0.5, 4);
    for(k=0; k<4; k++)
    {
        printf("%f ",w[k]);
    }
    printf("\n");
*/
    return 0;
}

double* mult_esc(double* vet, double a, int tam )
{
    int j;
    double *temp = (double*)malloc((tam) * sizeof(double));
    for(j = 0; j < tam; j++)
    {
        temp[j] = a*vet[j];
    }
    return temp;
}

double* som_vec(double* vet1, double* vet2, int tam )
{
    int j;
    double *temp = (double*)malloc((tam) * sizeof(double));
    for(j = 0; j < tam; j++)
    {
        temp[j] = vet1[j]+vet2[j];
    }
    return temp;
}


double** Covarian(double** mR,double* U,int tam, int lin)
{
    int i; int j; int k;// contadores
    double som; // var aux
    int n = lin-1; // n é o numero de linhas de R
    double **Q = alocarMatriz(tam,tam);


    for(j = 0; j < tam; j++) // para cada coluna j
    {
        for(k = 0; k < tam; k++) // para nada coluna k
        {
            som = 0; // inicializa variavel que armazena somatorio
            for(i = 0; i < n; i++) // para cada linha i
            {
                som += (mR[i][j]-U[j])*(mR[i][k]-U[k]); // calcula somatorio
            }
            Q[j][k] = som/(n-1); // finaliza calculo da covariancia e armazena na matriz
        }
    }
    return Q; // retorna matriz de covarianca completa de tamanho tam x tam
}

double* multMatr(double** m1, double* m2, int tam)
{
    int i; int j; int k;
    double *m3 = (double*)malloc((tam) * sizeof(double));

    for (i=0;i<tam; i++)
    {
        for (j=0; j<tam; j++)
        {
            m3[i] = m3[i] + (m1[i][j] * m2[j]);

        }

    }
    return m3;
}

double multVec(double* m1, double* m2, int tam)
{
    int i; int j; int k;
    double m3;

    for (i=0;i<tam; i++)
    {
        m3 = m3 + m1[i]*m2[i];
    }
    return m3;
}


double* varian(double** mR, double* U,int tam, int lin)
{
    int i; int j; // contadores
    double som; // var aux
    int n = lin-1; // n é o numero de linhas de R
    double *Sigma = (double*)malloc((tam) * sizeof(double));
    for(j = 0; j < tam; j++) // para cada coluna
    {
        som = 0; // inicializa variavel auxiliar do somatorio
        for(i = 0; i < n; i++) // para cada linha
        {
            som += pow((mR[i][j]-U[j]),2); // calcula variavel auxiliar para somatorio
        }
        Sigma[j] = som/(n-1); // finaliza o calculo da variancia da coluna de mR e armazena
    }
    return Sigma; // retorna vetor de variancias de tamanho tam
}

double* esper(double** mR,int tam,int lin)
{
    int i; int j; // contadores
    int n = lin-1; // n é o numero de linhas de R
    double *U = (double*)malloc((tam) * sizeof(double));
    for(j = 0; j < tam; j++) // para cada coluna
    {
        U[j] = 0; // inicializa elemento do vetor de esperanças
        for(i = 0; i < n; i++) // para cada linha
        {
            U[j] += mR[i][j]; // calcula somatorio da esperança
        }
        U[j] = U[j]/(n); // finaliza o calculo do elemento do vetor de esperanças e armazena
    }

    return U; // retorna vetor de esperanças de tamanho tam
}

double** matrR(double** P,int tam, int lin)
{
    int i; int j; // contadores
    double **mR = alocarMatriz(lin-1,tam); // Alocamos a matiz (Tem linha a menos pois não tem R definido para o primeiro elemento)
    for(j = 0; j < tam; j++) // para cada coluna
    {
        for(i = 1; i < 100; i++) // para cada linha
        {
            mR[i-1][j] = log(P[i][j]/P[i-1][j]); // calcula e armazena elemento da matriz de Retornos
        }
    }

    return mR; // retorna matriz de Retornos de tamanho lin-1 x tam
}

double* ResolveSis(double** L, int tam, double* b) // Resolve sistema L*x = b para uma matriz L obtida pela decomposição
{
    int i = 0; int j = 0; // inicializa contadores
    double ly; double lx; // declara variaveis auxiliares
    double y[tam]; // declara vetor y

    y[0] = (b[0])/(L[0][0]);

    for(i=1; i<tam; i++)
    {
        ly = 0;
        for(j=0; j<i; j++)
        {
            ly += (L[i][j])*(y[j]);
        }
        y[i] = (1/(L[i][i]))*(b[i] - ly);
    }
    double *x = (double*)malloc(tam * sizeof(double));
    x[tam-1] = y[tam-1]/L[tam-1][tam-1];
    for(i=tam-2; i>=0; i--)
    {
        lx = 0;
        for(j=i+1; j<tam; j++)
        {
            lx = lx + L[j][i]*x[j];
        }
        x[i] = (1/(L[i][i]))*(y[i]-lx);
    }

    return x;
}

double** coef (double **A, int tam)
{
    double Alin[tam][tam];
    int k; int j;

    for(k=0; k<tam; k++)
    {
        for(j=0; j<tam; j++)
        {
            Alin[k][j] = A[k][j];
        }
    }
    double **L = (double**) alocarMatriz(tam,tam);
    j = 0;
    k = 0;
    int i;

    while (Alin[k][k] > 0 && k < tam)
    {
        L[k][k]= sqrt(Alin[k][k]);

        for(i=k+1; i<tam; i++)
        {
            L[i][k]= Alin[i][k]/L[k][k];
            for(j=k+1; j<=i; j++)
            {
                Alin[i][j]= Alin[i][j]-L[i][k]*L[j][k];
            }
        }
        k++;
    }

    if (k < tam)
    {
        printf("FALHA NO ALGORITMO");
        char c = getc(stdin);
        exit(0);
    }

    return L;
}

double** alocarMatriz(int Linhas,int Colunas){ //Recebe a quantidade de Linhas e Colunas como Parâmetro

  int i,j; //Variáveis Auxiliares

  double **m = (double**)malloc(Linhas * sizeof(double*)); //Aloca um Vetor de Ponteiros

  for (i = 0; i < Linhas; i++){ //Percorre as linhas do Vetor de Ponteiros
       m[i] = (double*) malloc(Colunas * sizeof(double)); //Aloca um Vetor de Inteiros para cada posição do Vetor de Ponteiros.
       for (j = 0; j < Colunas; j++){ //Percorre o Vetor de Inteiros atual.
            m[i][j] = 0; //Inicializa com 0.
       }
  }
return m; //Retorna o Ponteiro para a Matriz Alocada
}



