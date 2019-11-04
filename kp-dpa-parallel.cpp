#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<iomanip>
#include<mpi.h>
using namespace std;
class KnapSolver
{
	int* a, * p, * w, * x, * send_arr, * recv_arr, max_w, obj;
public:
	void read(char* file_name);
	void solve();
};
void KnapSolver::read(char* file_name)
{
	ifstream g;
	int t = 0;
	g.open(file_name);
	if (!g)
	{
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}
	g >> obj;
	g >> max_w;
	w = (int*)malloc(obj * sizeof(int));
	if (w == NULL) { cerr << "Error : Your size is too much.\n"; exit(1); }
	p = (int*)malloc(obj * sizeof(int));
	if (p == NULL) { cerr << "Error : Your size is too much.\n"; exit(1); }
	while (!g.eof())
	{
		g >> p[t];
		g >> w[t];
		t++;
		if (t > obj)
			break;
	}
}
void KnapSolver::solve()
{
	int i, j, chunk = 10, size, rank, C1, C2;
	double start = 0, end = 0;
	C1 = max_w / 2;
	C2 = max_w - C1;
	a = (int*)malloc((max_w + 1) * obj * sizeof(int));
	if (a == NULL) { cerr << "Error : Your size is too much.\n"; exit(1); }
	x = (int*)malloc(obj * sizeof(int));
	if (x == NULL) { cerr << "Error : Your size is too much.\n"; exit(1); }
	MPI_Init(NULL, NULL); //initialize MPI operations
	start = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
	MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
	send_arr = (int*)malloc((max_w + 1) * sizeof(int));
	if (send_arr == NULL) { cerr << "Error : Your size is too much.\n"; exit(1); }
	recv_arr = (int*)malloc((max_w + 1) * sizeof(int));
	if (recv_arr == NULL) { cerr << "Error : Your size is too much.\n"; exit(1); }
	if (rank == 0) {
		for (i = 0; i < obj; i++)
		{
			for (j = 0; j < C1; j++)
			{
				if (j < w[i])
				{
					if (j == 0 || i == 0)
						send_arr[j] = a[j * obj + i] = 0;
					else
						send_arr[j] = a[j * obj + i] = a[j * obj + i - 1];
				}
				if (j >= w[i])
				{
					if (i == 0)
						send_arr[j] = a[j * obj + i] = p[i];
					else
					{
						int k = j - w[i];
						if (a[j * obj + i - 1] > (a[k * obj + i - 1] + p[i]))
							send_arr[j] = a[j * obj + i] = a[j * obj + i - 1];
						else
							send_arr[j] = a[j * obj + i] = a[k * obj + i - 1] + p[i];
					}

				}
			}
			if(i != obj)MPI_Send(send_arr, C1, MPI_INT, 1, 1, MPI_COMM_WORLD);
		}
	}
	if (rank == 1) {
		for (i = 0; i < obj; i++)
		{
			MPI_Recv(recv_arr, C1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			/*cout << endl<<"This is rank 1 receiving from rank 0 with i = "<<i<<endl;
			for (int p = 0; p < C1; p++)cout << recv_arr[p] << ",";
			cout << endl;*/
			for (int k = 0; k < C1; k++)a[k * obj + i] = recv_arr[k];
			for (j = C1; j < max_w + 1; j++)
			{
				if (j < w[i])
				{
					if (j == 0 || i == 0)
						a[j * obj + i] = 0;
					else
						a[j * obj + i] = a[j * obj + i - 1];
				}
				if (j >= w[i])
				{
					if (i == 0)
						a[j * obj + i] = p[i];
					else
					{
						int k = j - w[i];
						if (a[j * obj + i - 1] > (a[k * obj + i - 1] + p[i]))
							a[j * obj + i] = a[j * obj + i - 1];
						else
							a[j * obj + i] = a[k * obj + i - 1] + p[i];
					}

				}
			}
			/*if (i == 6) {
				cout << endl;
				for (j = 0; j < max_w+1; j++) {
					cout << a[j * obj + i] << ",";
				}
				cout << endl;
			}*/
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
		}
		cout << "\nThe Answer is = ";
		for (int i = 0; i < obj; i++)cout << x[i] << " ";*/
		//cout << "\n\nThe maximum value is = " << a[max_w * obj + obj - 1] << endl << endl;
	}
	end = MPI_Wtime();
	if (rank == 1)cout << "\nThe process took " << end - start << " seconds to run." << std::endl;
	MPI_Finalize(); //finalize MPI operations
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