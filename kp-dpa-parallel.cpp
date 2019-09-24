#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<iomanip>
#include<mpi.h>
//#include<conio.h>
using namespace std;
class KnapSolver
{
	int *a, *p, *w, *x, *send_arr, *recv_arr, max_w,obj;
	public:
        void read(char *file_name);
        void solve();
		void output();

};
void KnapSolver::read(char *file_name)
{
	ifstream g;
	int t=0;
	g.open(file_name);
	if(!g)
	{
      		cerr << "Error: file could not be opened" << endl;
      		exit(1);
        }
	g >> obj;
	g >> max_w;
	w=(int *)malloc(obj*sizeof(int));
	if(w == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
	p=(int *)malloc(obj*sizeof(int));
	if(p == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
	while(!g.eof())
	{
		g >> p[t];
		g >> w[t];
		t++;
		if(t > obj)
		  break;
	}
	printf("\nTotal number of object is %d.\n\n",obj);
    printf("The maximum weight is %d.\n\n",max_w);
}
void KnapSolver::solve()
{
	int i,j,chunk=10, size, rank;
	a=(int *)malloc((max_w+1)*obj*sizeof(int));
	if(a == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
	x=(int *)malloc(obj*sizeof(int));
	if(x == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
	MPI_Init(NULL, NULL); //initialize MPI operations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
    send_arr=(int *)malloc((max_w+1)*sizeof(int));
	if(send_arr == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
	recv_arr=(int *)malloc((max_w+1)*sizeof(int));
	if(recv_arr == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
    if (rank == 0){
		for(i=0;i<obj;i++)
		{  
			for(j=0;j<(max_w+1)/size;j++)
			{
				if(j<w[i])
				{
					if(j==0 || i == 0)
					  send_arr[j]=a[j*obj+i]=0;
					else 
					  send_arr[j]=a[j*obj+i]=a[j*obj+i-1];
				}
				if(j>=w[i])
				{
					if(i == 0)
					  send_arr[j]=a[j*obj+i]=p[i];
					else
					{
	                  int k=j-w[i];
					  if(a[j*obj+i-1]>(a[k*obj+i-1]+p[i]))
					    send_arr[j]=a[j*obj+i]=a[j*obj+i-1];
					  else 
					    send_arr[j]=a[j*obj+i]=a[k*obj+i-1]+p[i];
					}
				
				}
			}
		}
		MPI_Send(send_arr, (max_w+1)/size, MPI_INT, 1, 1, MPI_COMM_WORLD);
	}
	if (rank == 1){
		MPI_Recv(recv_arr, (max_w+1)/size, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//for(int k=0; k<(max_w+1)/size; k++)
		//	a[k*obj+i]=recv_arr[k];
		//cout<<"\n--------------------\n";
		//for(int k=0; k<(max_w+1); k++){
		//	cout<<recv_arr[k]<<",";
		//}
		//cout<<"\n--------------------\n";
		//getchar();
		//std::cin.get();
		//getch();
		for(i=0;i<obj;i++)
		{ 
			for(int k=0; k<(max_w+1)/size; k++)a[k*obj+i]=recv_arr[k]; 
			for(j=((max_w+1)/size);j<max_w+1;j++)
			{
				if(j<w[i])
				{
					if(j==0 || i == 0)
					  a[j*obj+i]=0;
					else 
					  a[j*obj+i]=a[j*obj+i-1];
				}
				if(j>=w[i])
				{
					if(i == 0)
					  a[j*obj+i]=p[i];
					else
					{
	                  int k=j-w[i];
					  if(a[j*obj+i-1]>(a[k*obj+i-1]+p[i]))
					    a[j*obj+i]=a[j*obj+i-1];
					  else 
					    a[j*obj+i]=a[k*obj+i-1]+p[i];
					}
				
				}
			}
		}
	}
	MPI_Finalize(); //finalize MPI operations
	int k=max_w;
	for(int i=obj-1;i>=0;i--)
	{
		if(i==0){
		  if(a[k*obj]==0) 
		    x[i]=0;
		  else
		    x[i] = 1;
		}else if(a[k*obj+i]!=a[k*obj+i-1])
		{
			x[i]=1;k=k-w[i];
		}
		else
		    x[i]=0;
	}
}
void KnapSolver::output()
{
#if 0
        for(int j = 0; j < max_w + 1; j ++) {
	  for(int i = 0; i < obj; i ++) {
	    printf("%d ", a[j * obj + i]);
	  }
	  printf("\n");
	}
#endif
	cout<<"\nThe Answer is = ";
	/*for(int i=0;i<max_w+1;i++)
	{
		for(int j=0;j<obj;j++)
		{
			cout<<setw(4)<<a[i*obj+j]<<" ";
		}
		cout<<endl;
	}*/
	for(int i=0;i<obj;i++)cout<<x[i]<<" ";
	cout<<"\n\nThe maximum value is = "<<a[max_w*obj+obj-1]<<endl;
}
int main(int argc, char* argv[])
{
	KnapSolver kp;
	char* str = NULL;
	int nt = 1;
	double start=0,end=0;
	if(argc >= 2) {
	  str = argv[1];
	} else {
	  fprintf(stderr, "usage: %s <input_file>\n", argv[0]);
	  exit(-1);
	}
	kp.read(str);
	kp.solve();
	kp.output();
}
