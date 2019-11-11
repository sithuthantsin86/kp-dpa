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
	double start = 0, end = 0;
	P1 = C / 2;
	P2 = C - P1;
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
						if (a[(i-1) * (C+1) + j] > (a[(i-1) * (C+1) + k] + p[i]))
							a[i * (C+1) + j] = a[(i-1) * (C+1) + j];
						else
							a[i * (C+1) + j] = a[(i-1) * (C+1) + k] + p[i];
					}
				}
			}
			if(i != N-1)MPI_Send(&a[i * (C+1) + 0], P1, MPI_INT, 1, 1, MPI_COMM_WORLD);
		}
	}
	if(rank == 1){
		for(i = 0; i < N; i++)
		{
			if(i != 0){
				MPI_Recv(&(a[(i-1) * (C+1) + 0]), P1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
						if (a[(i-1) * (C+1) + j] > (a[(i-1) * (C+1) + k] + p[i]))
							a[i * (C+1) + j] = a[(i-1) * (C+1) + j];
						else
							a[i * (C+1) + j] = a[(i-1) * (C+1) + k] + p[i];
					}
				}
			}
		}
		cout << "The maximum value is = " << a[C * N + N - 1] << endl;
	}
	/*int k = max_w;
	for (int i = obj - 1; i >= 0; i--)
	{
		if (i == 0) {
			if (a[k * obj] == 0)
				x[i] = 0;
			else
				x[i] = 1;
		}
		else if (a[k * obj + i] != a[k * obj + i - 1])
		{
			x[i] = 1; k = k - w[i];
		}
		else
			x[i] = 0;
	}*/
	end = MPI_Wtime();
	if (rank == 1)cout << "The process took " << end - start << " seconds to run." << std::endl;
	MPI_Finalize(); //finalize MPI operations
	delete[] p;
	delete[] a;
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
	kp.solve();
}