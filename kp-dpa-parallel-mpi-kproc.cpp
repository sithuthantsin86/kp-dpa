#include<iostream>
#include<algorithm>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<iomanip>
#include<time.h>
#include<new>
#include<mpi.h>
using namespace std;

class KnapSolver {
    int *a, *p, *w, *x, C, N;
public:
    void read(int rank, char* file_name);
    void solve();
};

void KnapSolver::read(int rank, char* file_name) {
    ifstream g;
    int count = 0;
    if (rank == 0) {
        g.open(file_name);
        if (!g) {
            cerr << "Error: file could not be opened" << endl;
            exit(1);
        }
        g >> N;
        g >> C;
        p = new (nothrow) int [N];
        if (p == nullptr)cout << "Error: memory could not be allocated for p.";
        w = new (nothrow) int [N];
        if (w == nullptr)cout << "Error: memory could not be allocated for w.";
        while (!g.eof()) {
            g >> p[count];
            g >> w[count];
            count++;
            if (count > N)break;
        }
        g.close();
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&C, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        p = new (nothrow) int [N];
        if (p == nullptr)cout << "Error: memory could not be allocated for p.";
        w = new (nothrow) int [N];
        if (w == nullptr)cout << "Error: memory could not be allocated for w.";
    }
    MPI_Bcast(p, N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, N, MPI_INT, 0, MPI_COMM_WORLD);
}
void KnapSolver::solve() {
    int i, j, pr1, pr2, ps1, ps2, size, rank, m, cnt_r1, cnt_r2, cnt_s1, cnt_s2;
    double start = 0, end = 0, startBT = 0, endBT = 0;
    x = new (nothrow) int [N+2];
    if (x == nullptr)
        cout << "Error: memory could not be allocated for x.";
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
    m = ceil((double) (C + 1) / (double) size);
    int nsize = m * size;
//    std::cout << " m = " << m << ", nsize = " << nsize << "C + 1 =" << C + 1 << "\n";
    a = new (nothrow) int [N * nsize];
    if (a == nullptr)
        cout << "Error: memory could not be allocated for a.";
    start = MPI_Wtime();
    for (i = 0; i < N; i++) {
         //cout<<"\nHello from rank "<< rank <<" of size "<< size <<".\n";
        if (i != 0 && rank != 0) {
            const int pbeg = m * rank - w[i];
            const int pend = m * rank + m - 1 - w[i];
            pr1 = floor((double) (pbeg) / (double) m);
            pr2 = floor((double) (pend) / (double) m);
            if (pr1 >= 0) {
                cnt_r1 = (m * pr1 + (m - 1))- pbeg + 1;
                MPI_Recv(&a[(i - 1) * nsize + pbeg], cnt_r1, MPI_INT, pr1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (pr1 != pr2 && pr2 >= 0 && pr2 < rank) {
                cnt_r2 = pend -(m * pr2) + 1;
                MPI_Recv(&a[(i - 1) * nsize + m * pr2], cnt_r2, MPI_INT, pr2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        for (j = rank * m; j < std::min((rank * m) + m, C + 1); j++) {
            if (j < w[i]) {
                if (j == 0 || i == 0)
                    a[i * nsize + j] = 0;
                else
                    a[i * nsize + j] = a[(i - 1) * nsize + j];
            }
            if (j >= w[i]) {
                if (i == 0)
                    a[i * nsize + j] = p[i];
                else {
                    int k = j - w[i];
                    a[i * nsize + j] = max(a[(i - 1) * nsize + j], a[(i - 1) * nsize + k] + p[i]);
                }
            }
            //if (i == N - 1 && j == C)cout << "\nAns = " << a[i * nsize + j] << ".\n";
        }
        if (i != N - 1 && rank < size - 1) {
            const int pbeg = m * rank + w[i + 1];
            const int pend = m * rank + m - 1 + w[i + 1];
            ps1 = floor((double)pbeg / (double)m);
            ps2 = floor((double)pend / (double)m);
            if (ps1 < size && ps1 > rank) {
                cnt_s1 = (m * ps1 + (m - 1)) - pbeg + 1;
                MPI_Send(&a[i * nsize + (m * rank)], cnt_s1, MPI_INT, ps1, 1, MPI_COMM_WORLD);
                //cout<<"\nSending "<<cnt_s1<< " size, starting from " << a[i * (C+1) + (m*rank)]<<" to p" << ps1 << " from p"<<rank<<" of "<<i+1<<"th object. (ps1)\n";

            }
            if (ps1 != ps2 && ps2 < size && ps2 > rank) {
                cnt_s2 = pend - m * ps2 + 1;
                MPI_Send(&a[i * nsize + m * ps2 - w[i + 1]], cnt_s2, MPI_INT, ps2, 1, MPI_COMM_WORLD);
                //cout<<"\nSending "<<cnt_s2<< " size, starting from "<<a[i * (C+1) + (m*rank+(m-1)-cnt_s2+1)]<<" to p" << ps2 << " from p"<<rank<<" of "<<i+1<<"th object. (ps2)\n";
            }
        }
    }
    end = MPI_Wtime();
    if (rank == 0) {
        //cout << "\nThe process took " << end - start << " seconds to run." << std::endl;
        cout << end - start << "\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    /////Back tracking algorithm./////
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    startBT = MPI_Wtime();
    int j, k = C, EndPointOfP = 0, i = N-1;
    //for(int i = N - 1; i >= 0; i--){
    while(i >= 0){
        if(rank < size-1)
        {
            MPI_Recv(&x, N+2, MPI_INT, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(x[N+1] < m*rank) goto label_1;
        }
        if(i == 0)
        {
            if (a[i * (C+1) + x[N+1]] == 0)
            {
                x[i] = 0;
                i--;
            }
            else
            {
                x[i] = 1;
                i--;
            }
        }
        
        label_1:
        MPI_Send(&x, N+2, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
        else if (a[i * (C+1) + k] != a[(i-1) * (C+1) + k])
        {
            x[i] = 1;
            if(k-w[i]<P1)
            {
                EndPointOfP1 = i;
                MPI_Send(&i, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                break;
            }
            k = k - w[i];
        }
        else
            x[i] = 0;
    }
    endBT = MPI_Wtime();
    /*startBT = MPI_Wtime();
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
                                    MPI_Send(&i, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                                    break;
                            }
                            k = k - w[i];
                    }
                    else
                            x[i] = 0;
            }
            MPI_Send(&x[EndPointOfP1], N - EndPointOfP1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
    if(rank == 0)
    {
            int j;
            MPI_Recv(&j, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            k=k-w[j];
            for (int i = j-1; i >= 0; i--)
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
            MPI_Recv(&x[j], N - j, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //cout<<"The solution vector X = {";
            //for(int i=0; i<=j; i++){
            //	cout<<x[i];
            //	if(i<j)cout<<",";
            //}
            //cout<<"}\n";
    }
    endBT = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);*/
    /*if(rank == 1){
            //cout << "\nThe process took " << end - start << " seconds to run." << std::endl;
            //MPI_Send(end-start, 1, MPT_D);
            double runtime = end - start;
            //cout << runtime << "\n";
            MPI_Send(&runtime, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }*/
    /*if(rank == 0){
            double runtime_all;
            MPI_Recv(&runtime_all, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //cout << endBT - startBT <<"\t";
            runtime_all = runtime_all+ (endBT - startBT);
            cout << runtime_all <<"\n";
    //	cout << "The backtrack algrithm process took " << end - start << " seconds to run." << std::endl;
    }*/
    //MPI_Barrier(MPI_COMM_WORLD);
    //cout<<"x from P1 = ";
    //if(rank==1)for(int i=0; i<N; i++){
    //cout<<x[i]<<",";
    //cout << "The backtrack algrithm process took " << endBT - startBT << " seconds to run." << std::endl;
    //cout << endBT - startBT << "\n";
    //}
    //cout<<endl;
    delete[] p;
    delete[] a;
    delete[] x;
}

int main(int argc, char* argv[]) {
    KnapSolver kp;
    char* str = NULL;
    int rank;
    MPI_Init(&argc, &argv); //initialize MPI operations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        if (argc >= 2) {
            str = argv[1];
        } else {
            fprintf(stderr, "usage: %s <input_file>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    kp.read(rank, str);
    //clock_t begin = clock();
    kp.solve();
    MPI_Finalize(); //finalize MPI operations

    //clock_t end = clock();
    //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //cout << "The process took " << time_spent << " seconds to run.\n";
    //cout << time_spent << "\t";
}
