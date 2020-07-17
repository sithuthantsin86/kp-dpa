#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string>
#include<sstream>
#include<sys/stat.h>
using namespace std;
void printPair(FILE* f, int *p, int *w, int N, int C)
{
	fprintf(f, "%d\t%d\n", N, C);
	for(int i=0;i<N;i++)
		fprintf(f,"%d\t%d\n", p[i], w[i]);
}
string file_name(string method, int N, int C, int num_instance)
{
	ostringstream os;
	os<<method<<"/";
	os <<method<<"_"<<N<<"_"<<C<<"_";
	if(num_instance<10)os <<"0"<<"0"<<num_instance;
	if(num_instance>=10 && num_instance<100)os <<"0"<<num_instance;
	if(num_instance>=100)os<<num_instance;
	return os.str();
}
void generator(int number_of_items, int capacity, int num_instances)
{
	time_t t;
	struct stat st;
	srand (time(&t));
	FILE *a;
	int *p, *w, K=0;
	p=(int *)malloc(number_of_items*sizeof(int));
	if(p == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
	w=(int *)malloc(number_of_items*sizeof(int));
    if(w == NULL){cerr<<"Error : Your size is too much.\n";exit(1);}
	if(stat("NxC_test_MPI", &st) == 0)system("rm -r NxC_test_MPI");
	K = min(max(4 * capacity/number_of_items, 10), capacity);
	cout << "\nK = " << K << "\n";
	for(int j=0;j<num_instances;j++){
			string M="NxC_test_MPI";
			//if(stat("UC", &st) == 0)system("rm -r UC");
			mkdir("NxC_test_MPI", 0777);
			a=fopen(file_name(M, number_of_items, capacity, j+1).c_str(), "w");
			for(int i=0;i<number_of_items;i++){
				p[i] = (rand() % K + 1);
				w[i] = (rand() % K + 1);
			}
			printPair(a, p, w, number_of_items, capacity);
			fclose(a);
	}
}
int main(int argc, char* argv[])
{
	int num_instances;
	int number_of_items, R, capacity;
	if(argc >= 2){
		number_of_items = atoi(argv[1]);
		if(argc == 4){
			capacity = atoi(argv[2]);
			num_instances = atoi(argv[3]);
		}
	}
	else{
		fprintf(stderr,"Usage: %s <number of items> <capacity> <How many instances>\n", argv[0]);
		exit(-1);
	}
	generator(number_of_items, capacity, num_instances);
	return 0;
}
