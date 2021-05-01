#include <iostream>
#include <math.h>
#include <mpi.h> // API MPI

double funcao2(double x) {
    return (2 / sqrt(M_PI)) * exp(-x * x);
}

double funcao(double x) {
    return (cos(0 * M_PI * x) * exp((cos(M_PI * x) - 1) / (2 * M_PI * 1)));
}

bool debug = false ;
double xIni = 0.0;
double xFim = 1.0;
double soma = 0.0;
int n = 1000;
double h;
int m;
int world_size; // número total de processos
int world_rank; // ID (rank) do processo

int main() {

    MPI_Init(NULL, NULL); // Inicializa o MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    h = (xFim - xIni) / (double) n;
    m = (int) n / world_size;
    if (world_rank == 0) {
        soma = funcao(xIni) / 2;
        for (int i = 1; i <= m; i++)
            soma += funcao(i * h);
        soma = soma * h;
        printf("soma(%d) = %20.12f", world_rank,  soma);

        double SomaDosProcessos[world_size];
        SomaDosProcessos[0] = soma;
        if (debug) {
            printf("\nSoma do processo %d eh = %20.12f", 0, SomaDosProcessos[0]);
        }
        for (int i = 1; i < world_size; i++) {
            MPI_Recv(&SomaDosProcessos[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (debug) {
                printf("\nSoma do processo %d eh = %20.12f", 0, SomaDosProcessos[0]);
            }
        }
        double somaTotal = 0.0;
        for (int i = 0; i < world_size; i++)
            somaTotal += SomaDosProcessos[i];
        printf("\nValor da integral eh    = %20.12f", somaTotal);
        printf("\nValor da funcao erro eh = %20.12f", erf(xFim));

    } else if (world_rank == world_size - 1) {
        for (int i = (world_size - 1) * m + 1; i < n; i++)
            soma += funcao(i * h);
        soma += funcao(xFim) / 2;
        soma = soma * h;
        printf("soma(%d) = %20.12f", world_rank, soma);
        MPI_Send(&soma, 1, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD);

    } else {
        for (int i = world_rank * m + 1; i <= (world_rank + 1) * m; i++)
            soma += funcao(i * h);
        soma = soma * h;
        printf("soma(%d) = %20.12f", world_rank, soma);
        MPI_Send(&soma, 1, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
/*
 * Rode no terminal
 * $ mpic++ 334gather.cpp -lgsl -lgslcblas -o 334gather
 * $ mpirun -np 6 -oversubscribe 334gather
 * Processo 0 soma = 15.000000
 *
 * */
