#include<omp.h>
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<iomanip>
#include<time.h>
#include<new>
using namespace std;
class KnapSolver
{
	int *a, *p, *w, *x, C, N;
public:
	void read(char* file_name);
	void solve();
	void backTrack();

};
void KnapSolver::read(char* file_name)
{
	ifstream g;
	int count = 0;
	g.open(file_name);
	if (!g)
	{
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}
	g >> N;
	g >> C;
	p = new (nothrow) int [N];
	if (p == nullptr)cout << "Error: memory could not be allocated for p.";
	w = new (nothrow) int [N];
	if (w == nullptr)cout << "Error: memory could not be allocated for w.";
	while (!g.eof())
	{
		g >> p[count];
		g >> w[count];
		count++;
		if(count > N)break;
	}
}
void KnapSolver::solve()
{
	int i, j;
	a = new (nothrow) int [N * (C+1)];
	if (a == nullptr)cout << "Error: memory could not be allocated for a.";
	x = new (nothrow) int [N];
	if (x == nullptr)cout << "Error: memory could not be allocated for x.";
	for (i = 0; i < N; i++)
	{
		#pragma omp parallel for
		for (j = 0; j < C + 1; j++)
		{
			if (j < w[i])
			{
				if (j == 0 || i == 0)
					a[i * (C+1) + j] = 0;
				else
					a[i * (C+1) + j] = a[(i-1) * (C+1) + j];
			}
			if (j >= w[i])
			{
				if (i == 0)
					a[i * (C+1) + j] = p[i];
				else
				{
					int k = j - w[i];
					if (a[(i-1) * (C+1) + j] > (a[(i-1) * (C+1) + k] + p[i]))
						a[i * (C+1) + j] = a[(i-1) * (C+1) + j];
					else
						a[i * (C+1) + j] = a[(i-1) * (C+1) + k] + p[i];
				}

			}
		}
	}
	//cout << "The maximum value is = " << a[C * N + N - 1] << endl;
}
void KnapSolver::backTrack()
{
	int k = C;
	#pragma omp parallel for
	for (int i = N - 1; i >= 0; i--)
	{
		if (i == 0)
		{
			if (a[i * (C+1) + k] == 0)
				x[i] = 0;
			else
				x[i] = 1;
		}
		else if (a[i * (C+1) + k] != a[(i-1) * (C+1) + k])
		{
			x[i] = 1;
			k = k - w[i];
		}
		else
			x[i] = 0;
	}
	//cout<<endl;
	//for(int i=0; i<N; i++)cout<<"x ="<<x[i]<<",";
	delete[] p;
	delete[] a;
	delete[] x;
}
int main(int argc, char* argv[])
{
	KnapSolver kp;
	char* str = NULL;
	int nt = 1;
	double start1=0, end1=0, start2=0,end2=0, time_all=0;
	if(argc >= 3) {
	  str = argv[1];
	  if(argc == 3){
	    nt = atoi(argv[2]);
	  }
	} else {
	  fprintf(stderr, "usage: %s <input_file> <number_of_threads>\n", argv[0]);
	  exit(-1);
	}
	kp.read(str);
	omp_set_num_threads(nt);
	//clock_t begin = clock();
	start1 = omp_get_wtime();
	kp.solve();
	end1 = omp_get_wtime();
	//clock_t end = clock();
	//double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//clock_t beginBT = clock();
	start2 = omp_get_wtime();
	kp.backTrack();
	end2 = omp_get_wtime();
	//clock_t endBT = clock();
	//double time_spent_BT = (double)(endBT - beginBT) / CLOCKS_PER_SEC;
	//cout << "The process took " << time_spent << " seconds to run.\n";
	cout << (end1-start1) + (end2-start2) << "\t";
}