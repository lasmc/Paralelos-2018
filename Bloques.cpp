#include <iostream>
#include <ctime>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

using namespace std;

#define MAX 1000

void mult_matriz()
{
	int ** f = new int * [MAX];
	int ** c = new int * [MAX];
	int ** resultado = new int * [MAX];

	for(int i=0; i<MAX; i++)
    {
        f[i]=new int [MAX];
        c[i]=new int [MAX];
        resultado[i]=new int [MAX];
        for(int j=0; j<MAX; j++)
        {
            f[i][j] = rand()%1000;
            c[i][j] = rand()%1000;
            resultado[i][j]=0;
        }
    }

	for(int i=0; i<MAX; i++)
	{
		for(int j=0; j<MAX; j++)
		{
			for(int k=0; k<MAX; k++)
			{
				resultado[i][j] += f[i][k] * c[k][j];
			}
		}
	}
}

void mult_matriz_bloques()
{
	int tam_bloq = MAX/4;
	int **f = new int * [MAX];
	int **c = new int * [MAX];
	int ** resultado = new int * [MAX];

	for(int i=0; i<MAX; i++)
    {
        f[i]=new int [MAX];
        c[i]=new int [MAX];
        resultado[i]=new int [MAX];
        for(int j=0; j<MAX; j++)
        {
            f[i][j]=rand()%1000;
            c[i][j]=rand()%1000;
            resultado[i][j]=0;
        }
    }

    for(int x=0; x<MAX; x+=tam_bloq)
    {
    	for(int y=0; y<MAX; y+=tam_bloq)
    	{
    		for(int z=0; z<MAX; z+=tam_bloq)
    		{
    			for(int i=x; i<x+tam_bloq; i++)
    			{
    				for(int j=y; j<y+tam_bloq; j++)
    				{
    					for(int k=z; k<z+tam_bloq; k++)
    					{
    						resultado[i][j]+=f[i][k]*c[k][j];
    					}
    				}
    			}
    		}
    	}
    }
}

int main()
{
	cout<<"Tiempo de multiplicacion normal:\n";
	clock_t star,end;
    star=clock();
	mult_matriz();
	end=clock();
	cout<<"Con "<<MAX<<" datos es: "<<(end-star)/(double)CLOCKS_PER_SEC<<" s."<<endl;
	cout<<"\nTiempo de multiplicacion por bloques:\n";
    star = clock();
	mult_matriz_bloques();
	end = clock();
    cout<<"Con "<<MAX<<" datos es: "<<(end-star)/(double)CLOCKS_PER_SEC<<" s."<<endl;
	return 0;
}
