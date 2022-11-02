#include <mpi.h>
#include <iostream>
#include <omp.h>

#define MIN(a,b)        ((a)>=(b)?(b):(a))

double foo( const double x )
{
        return 3.*x*x + 2.*x + 1.;
}


int main()
{
using std::cout;
using std::endl;

	int rank, numproc;

	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numproc);


	if( rank == 0 )
	{	cout << "rank : " << rank << endl;
		cout << "proc : " << numproc << endl;
		cout << "rank/proc test end ----------\n" << endl;
	}

	if( rank == 0 ){ cout << "Cyclic Distribution Test" << endl; }

	int TaskSize = 10;
	int idx=0;

	for(int i=rank;i<TaskSize;i=rank + (numproc * idx))
	{
		cout << "rank ( " << rank << " ) / i : " << i << endl;
		idx++;
	}
	if( rank == 0 ){ cout << "Cyclic Distribution Test End ----------\n" << endl; }


	double a = 0, b = 3.;
	int div = 1000000000;
	double dx = (b-a)/div;
	double partial_sum;
	double total_sum = 0.;
	double x1,x2;
	double h1,h2;
	double area;

	double elapsed_t;
	int thc;
	
	elapsed_t = -omp_get_wtime();

	area = 0.;

	////	////	////	BLOCK DISTRIBUTION
	int L = (div-1)/numproc;
	int R = (div-1)%numproc;

	int ista = L*rank + MIN(R,rank);
	int iend = ista + L - 1;
	if( R > rank ){ iend++; }

	for(int i=ista;i<iend+1;i++)
	{
		x1 = a + dx*i;
		x2 = x1 + dx;

		h1 = foo(x1);
		h2 = foo(x2);
	
		total_sum += (h2+h1)*dx/2.;

		idx++;
	}
	// MPI_Allreduce( h_matrix_tmp[n]->data, sp_sys->scf_h_matrix_vs_sp_ion_monopole[n]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
	// Allreduce
	MPI_Allreduce( &total_sum, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	elapsed_t += omp_get_wtime();

	//printf("Answer %lf\n",total_sum);       // should be ~ 39.000 // simple exe with time ... threading is much slower ...
	printf("Answer %lf\n",area);       // should be ~ 39.000 // simple exe with time ... threading is much slower ...
	printf("threads in use : %d / elapsed time : %lf\n",thc,elapsed_t);

	MPI_Finalize();

	return 0;
}

/*
	mpiicpc -O3 -std=c++11 *.cpp -qopenmp
*/


