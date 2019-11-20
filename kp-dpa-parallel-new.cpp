#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<iomanip>
#include<time.h>
#include<new>
#include<mpi.h>
using namespace std;
class KnapSolver
{
	int *a, *p, *w, *x, C, N;
public:
	void read(char* file_name);
	void solve();
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
	int i, j, P1, P2, size, rank;
	double start = 0, end = 0, startBT = 0, endBT = 0;
	P1 = (C+1) / 2;
	P2 = (C+1) - P1;
	a = new (nothrow) int [N * (C+1)];
	if (a == nullptr)cout << "Error: memory could not be allocated for a.";
	x = new (nothrow) int [N];
	if (x == nullptr)cout << "Error: memory could not be allocated for x.";
	MPI_Init(NULL, NULL); //initialize MPI operations
	start = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
	MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
	if(rank == 0){
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < P1; j++)
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
						a[i * (C+1) + j] = max(a[(i-1) * (C+1) + j], a[(i-1) * (C+1) + k] + p[i]);
					}
				}
			}
			if(i != N-1){
				MPI_Send(&a[i * (C+1) + (P1-min(w[i+1],P1))], min(w[i+1], P1), MPI_INT, 1, 1, MPI_COMM_WORLD);
			}
		}
	}
	if(rank == 1){
		for(i = 0; i < N; i++)
		{
			if(i != 0){
				MPI_Recv(&(a[(i-1) * (C+1) + (P1-min(w[i], P1))]), min(w[i], P1), MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				/*cout <<"This is rank 1 receiving from rank 0 with i = "<<i<<endl;
				for (int p = 0; p < P1; p++)cout << a[(i-1) * (C+1) + p] << ",";
				cout<<endl;*/
			}
			for(j = P1; j < C + 1; j++){
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
						a[i * (C+1) + j] = max(a[(i-1) * (C+1) + j], a[(i-1) * (C+1) + k] + p[i]);
					}
				}
			}
		}
		//cout << "The maximum value is = " << a[C * N + N - 1] << endl;
	}
	end = MPI_Wtime();
	if (rank == 1)cout << "The process took " << end - start << " seconds to run." << std::endl;
	//if (rank == 1)cout << end - start << "\n";
	MPI_Barrier(MPI_COMM_WORLD);
	/////Back tracking algorithm./////
	startBT = MPI_Wtime();
	int k = C, EndPointOfP1 = 0;
	if(rank == 1)
	{
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
				if(k-w[i]<P1)
				{
					EndPointOfP1 = i;
					i=i-1;
					MPI_Send(&i, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
					//MPI_Bcast(&x,N,MPI_INT,1,MPI_COMM_WORLD);
					break;
				}
				k = k - w[i];
			}
			else
				x[i] = 0;
		}
	}
	if(rank == 0)
	{
		int j;
		MPI_Recv(&j, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		k=k-w[j+1];
		for (int i = j; i >= 0; i--)
		{
			if(i == 0)
			{
				if (a[i * (C+1) + k] == 0)
					x[i] = 0;
				else
					x[i] = 1;
			}
			else 
			{
				if(a[i * (C+1) + k] != a[(i-1) * (C+1) + k])
				{
					x[i] = 1;
					k = k - w[i];
				}
				else
				{
					x[i] = 0;
				}
			}
		}
		cout<<"x from P0 = ";
		for(int i=0; i<=j; i++)cout<<x[i]<<",";
	}
	endBT = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==1)for(int i=EndPointOfP1; i<N; i++){
		cout<<"x from P1 = ";
		cout<<x[i]<<",";
	}
	cout<<endl;
	MPI_Finalize(); //finalize MPI operations
	delete[] p;
	delete[] a;
	delete[] x;
}
int main(int argc, char* argv[])
{
	KnapSolver kp;
	char* str = NULL;
	int nt = 1;
	if (argc >= 2) {
		str = argv[1];
	}
	else {
		fprintf(stderr, "usage: %s <input_file>\n", argv[0]);
		exit(-1);
	}
	kp.read(str);
	//clock_t begin = clock();
	kp.solve();
	//clock_t end = clock();
	//double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//cout << "The process took " << time_spent << " seconds to run.\n";
	//cout << time_spent << "\t";
}