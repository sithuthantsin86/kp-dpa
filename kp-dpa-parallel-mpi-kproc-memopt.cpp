///Dynamic programming algorithm for Knapsack Problem
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
    int *a, *p, *w, *x, *buff_s1, *buff_s2, *buff_r, C, N;
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
    int i, j, pr1, pr2, ps1, ps2, size, rank, m, cnt_r1 = 0, cnt_r2 = 0, cnt_rAll, cnt_s1, cnt_s2, z, l, v=0;
    double start = 0, end = 0, startBT = 0, endBT = 0;
    x = new (nothrow) int [N+2];
    if (x == nullptr)cout << "Error: memory could not be allocated for x.";
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
    m = ceil((double) (C + 1) / (double) size);
    int nsize = m * size;
    //std::cout << " m = " << m << ", nsize = " << nsize << ", C + 1 =" << C + 1 << "\n";
    a = new (nothrow) int [N * m];
    if (a == nullptr)
        cout << "Error: memory could not be allocated for a.";
    buff_s1 = new (nothrow) int [m];
    if (buff_s1 == nullptr)
        cout << "Error: memory could not be allocated for buff_s1.";
    buff_s2 = new (nothrow) int [m];
    if (buff_s2 == nullptr)
        cout << "Error: memory could not be allocated for buff_s2.";
    buff_r = new (nothrow) int [m];
    if (buff_r == nullptr)
        cout << "Error: memory could not be allocated for buff_r.";
    start = MPI_Wtime();
    for (i = 0; i < N; i++) {
        //cout<<"\nHello from rank "<< rank <<" of size "<< size <<".\n";
        if (i != 0 && rank != 0) {
            const int pbeg = m * rank - w[i]; ///Finding the beginning point of receiving process. 
            const int pend = m * rank + m - 1 - w[i];  ///Finding the ending point of receiving process.
            pr1 = floor((double) (pbeg) / (double) m);
            pr2 = floor((double) (pend) / (double) m);
            if (pr1 >= 0) {
                cnt_r1 = (m * pr1 + (m - 1))- pbeg + 1;
                MPI_Recv(&buff_r[0], cnt_r1, MPI_INT, pr1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //cout<<"\nReceiving "<<cnt_r1<< " size, starting from " << buff_r[0]<<" to p" << rank << " from p"<<pr1<<" for "<<i+1<<"th object. (pr1)\n";
            }
            if (pr1 != pr2 && pr2 >= 0 && pr2 < rank) {
                cnt_r2 = pend -(m * pr2) + 1;
                MPI_Recv(&buff_r[cnt_r1], cnt_r2, MPI_INT, pr2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //cout<<"\nReceiving "<<cnt_r2<< " size, starting from " << buff_r[cnt_r1]<<" to p" << rank << " from p"<<pr2<<" for "<<i+1<<"th object. (pr2)\n";
            }
            cnt_rAll = cnt_r1 + cnt_r2;
            cnt_r1 = cnt_r2 = 0;
            //cout<<"I'm rank "<<rank<<" in i = "<<i<<" and I'm revieving "<<cnt_rAll<<" size.\n";
        }
        for (j = 0; j < m; j++) {
            if ((j + (m*rank)) < w[i] && (j + (m*rank)) < C+1) {
                if ((j + (m*rank)) == 0 || i == 0){
                    a[i * m + j] = 0;
                    //cout<<"\n\n!Work at rank = "<<rank<<", i = "<<i<<", j = "<<j<<", weight = ["<<(j + (m*rank))<<" and that is 0.\n\n";
                }
                else{
                    a[i * m + j] = a[(i - 1) * m + j];
                    //cout<<"\n\n!Work at rank = "<<rank<<", i = "<<i<<", j = "<<j<<", weight = ["<<(j + (m*rank))<<" and that is "<<a[(i - 1) * m + j]<<".\n\n";
                }
            }
            if ((j + (m*rank)) >= w[i] && (j + (m*rank)) < C+1) {
                if (i == 0){
                    a[i * m + j] = p[i];
                    //cout<<"\n\n!Work at rank = "<<rank<<", i = "<<i<<", j = "<<j<<", [ i==0 and fill p[i] = "<<p[i]<<"].\n\n";
                }
                else {
                    int q = j - w[i];
                    if(cnt_rAll > 0 && v < cnt_rAll){
                        a[i * m + j] = max(a[(i - 1) * m + j], buff_r[v] + p[i]);
                        v++;
                        //cout<<"\n\n!Work at rank = "<<rank<<", i = "<<i<<", j = "<<j<<", ["<<a[(i - 1) * m + j]<<" or "<<buff_r[j]<<"+"<<p[i]<<"].\n\n";
                    }
                    else{
                        a[i * m + j] = max(a[(i - 1) * m + j], a[(i - 1) * m + q] + p[i]);
                        //cout<<"\n\n!!Work at rank = "<<rank<<", i = "<<i<<", j = "<<j<<", ["<<a[(i - 1) * m + j]<<" or "<<a[(i - 1) * m + q]<<"+"<<p[i]<<"].\n\n";
                    }

                }
            }
           if (rank == size-1 && i == N - 1 && (j + (m*rank)) == C)cout << "\n\nThe optimal value = " << a[i * m + j] << ".\n\n";
        }
        v=0;
        if (i != N - 1 && rank < size - 1) {
            const int pbeg = m * rank + w[i + 1];
            const int pend = m * rank + m - 1 + w[i + 1];
            ps1 = floor((double)pbeg / (double)m);
            ps2 = floor((double)pend / (double)m);
            l=0;
            if (ps1 < size && ps1 > rank) {
                cnt_s1 = (m * ps1 + (m - 1)) - pbeg + 1;
                for(z=0; z<cnt_s1; z++){
                    buff_s1[z] = a[i * m + l];
                    l++;
                }
                MPI_Send(&buff_s1[0], cnt_s1, MPI_INT, ps1, 1, MPI_COMM_WORLD);
                //cout<<"\nSending "<<cnt_s1<< " size, starting from " << a[i * m + 0]<<" to p" << ps1 << " from p"<<rank<<" of "<<i+1<<"th object. (ps1)\n";
            }
            if (ps1 != ps2 && ps2 < size && ps2 > rank){
                cnt_s2 = pend - m * ps2 + 1;
                l = m * ps2 - w[i + 1] - (m*rank);
                for(z=0; z<cnt_s2; z++){
                    buff_s2[z] = a[i * m + l];
                    l++;
                }
                MPI_Send(&buff_s2[0], cnt_s2, MPI_INT, ps2, 1, MPI_COMM_WORLD);
                //cout<<"\nSending "<<cnt_s2<< " size, starting from "<<a[i * m + (m * ps2 - w[i + 1] - (m*rank))]<<" to p" << ps2 << " from p"<<rank<<" of "<<i+1<<"th object. (ps2)\n";
            }
        }
    }
    end = MPI_Wtime();
    /*if (rank == 0) {
        //cout << "\nThe process took " << end - start << " seconds to run." << std::endl;
        cout <<"\nThe direct calculation runtime = "<< end - start << "\n";
    }*/
    MPI_Barrier(MPI_COMM_WORLD);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    /////Back tracking algorithm./////
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    int k;
    startBT = MPI_Wtime();
    if(rank == size-1)
    {
        k = x[N] = C;
        i = x[N+1] = N-1;
    }
    //cout<<"Process "<<rank<<" of "<<size<<".\n";
    if(rank < size-1)
    {
        MPI_Recv(x, N+2, MPI_INT, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //cout<<"\nProcess "<<rank<<" received i="<<x[N+1]<<" and k="<<x[N]<<" from process "<<rank+1<<".\n";
        k = x[N];
        i = x[N+1];
    }
    //cout<<"\nProcess "<<rank<<" is doing while with i="<<i<<" and k="<<k<<" and StartingPointOfProcess="<<rank*m<<".\n";
    while(i >= 0 && k >= rank*m)
    {
        if(i == 0)
        {
            if (a[i * m + (k - (rank * m))] == 0)
            {
                x[i] = 0;
                i--;
                //cout<<"\nProcess "<<rank<<" is adding "<<x[i+1]<<" to the "<<i+1<<"th place and next k="<<x[N]<<" and i="<<x[N+1]<<".\n";
            }
            else
            {
                x[i] = 1;
                i--;
                //cout<<"\nProcess "<<rank<<" is adding "<<x[i+1]<<" to the "<<i+1<<"th place and next k="<<x[N]<<" and i="<<x[N+1]<<".\n";
            }
        }
        else if (a[i * m + (k - (rank * m))] != a[(i-1) * m + (k - (rank * m))])
        {
            x[i] = 1;
            k = k - w[i];
            i--;
            //cout<<"\nProcess "<<rank<<" is adding "<<x[i+1]<<" to the "<<i+1<<"th place and next k="<<x[N]<<" and i="<<x[N+1]<<".\n";
        }
        else
        {
            x[i] = 0;
            i--;
            //cout<<"\nProcess "<<rank<<" is adding "<<x[i+1]<<" to the "<<i+1<<"th place and next k="<<x[N]<<" and i="<<x[N+1]<<".\n";
        }
    }
    x[N] = k;
    x[N+1] = i;
    if(rank != 0)
    {
        MPI_Send(x, N+2, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
        //cout<<"\nProcess "<<rank<<" sent k="<<x[N]<<" and i="<<x[N+1]<<" to process "<<rank-1<<".\n";
        /*cout<<"\nThe solution vector in "<<i<<"th step = {";
        for(int i=0; i<N+1; i++)
        {
            cout<<x[i];
            if(i!=N-1)cout<<", ";
        }
        cout<<"}.\n";*/
    }
    MPI_Barrier(MPI_COMM_WORLD);
    endBT = MPI_Wtime();
    if(rank == 0)
    {
        //cout <<"\nThe inverse calculation runtime = "<< endBT - startBT << "\n";
        cout << (end - start) + (endBT - startBT) << "\n";
        /*cout<<"\nThe solution vector = {";
        for(int i=0; i<N; i++){
            cout<<x[i];
            if(i!=N-1)cout<<", ";
        }
        cout<<"}.\n";*/
    }
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
    //cout << "\nThe process took " << time_spent << " seconds to run.\n";
}
