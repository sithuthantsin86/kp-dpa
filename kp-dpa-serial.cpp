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
	void output();

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
	delete[] p;
}
void KnapSolver::output()
{
#if 0
	for (int j = 0; j < max_w + 1; j++) {
		for (int i = 0; i < obj; i++) {
			printf("%d ", a[j * obj + i]);
		}
		printf("\n");
	}
#endif
	//cout << "\nThe Answer is = ";
	/*for(int i=0;i<max_w+1;i++)
	{
		for(int j=0;j<obj;j++)
		{
			cout<<setw(4)<<a[i*obj+j]<<" ";
		}
		cout<<endl;
	}*/
	//for (int i = 0; i < obj; i++)cout << x[i] << " ";
	/*cout<<"\n\nThe table is;\n\n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < C + 1; j++)
		{
			cout<<a[i * (C+1) + j]<<"\t,";
		}
		cout<<endl;
	}*/
	//cout << "\n\nThe maximum value is = " << a[C * N + N - 1] << endl;
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
	clock_t begin = clock();
	kp.solve();
	kp.output();
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	cout << "\nThe process took " << time_spent << " seconds to run.\n" << std::endl;
}