#include <iostream>
#include <math.h>
#include <mpi.h> // API MPI

double funcao2(double x) {
    return (2 / sqrt(M_PI)) * exp(-x * x);
}

double funcao(int n, double v, double x) {
    return (2 * cos(n * M_PI * x) * exp((cos(M_PI * x) - 1) / (2 * M_PI * v)));
}

double soma1(int MM, double v, double x, double t, double C[]) {
    double soma = 0.0;
    for (int i = 1; i <= MM; i++)
        soma += C[i] * exp(-i * i * M_PI * M_PI * v * t) * i * sin(i * M_PI * x);
    return soma;
}

double soma2(int MM, double v, double x, double t, double C[]) {
    double soma = 0.0;
    soma += C[0];
    for (int i = 1; i <= MM; i++)
        soma += C[i] * exp(-i * i * M_PI * M_PI * v * t) * cos(i * M_PI * x);
    return soma;
}


bool debug = false;
double xIni = 0.0;
double xFim = 1.0;
double soma = 0.0;
int MM = 20;
int n;
double v = 1.0;
int NN = 1000;
double h;
int m;
int world_size; // número total de processos
int world_rank; // ID (rank) do processo

double integra(int n, double v, int NN) {

    soma = 0.0;

    h = (xFim - xIni) / (double) NN;
    m = (int) NN / world_size;
    if (world_rank == 0) {

        double somaTotal = 0.0;
        soma = funcao(n, v, xIni) / 2;
        for (int i = 1; i <= m; i++)
            soma += funcao(n, v, i * h);
        soma = soma * h;
        if (debug) {
            printf("soma(%d) = %20.12f", world_rank, soma);
        }

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

        for (int i = 0; i < world_size; i++)
            somaTotal += SomaDosProcessos[i];
        if (debug) {
            printf("\nValor da integral eh  C%d  = %20.12f", n, somaTotal);
        }
        //printf("\nValor da funcao erro eh = %20.12f", erf(xFim));
        return somaTotal;

    } else if (world_rank == world_size - 1) {
        for (int i = (world_size - 1) * m + 1; i < NN; i++)
            soma += funcao(n, v, i * h);
        soma += funcao(n, v, xFim) / 2;
        soma = soma * h;
        if (debug) {
            printf("soma(%d) = %20.12f", world_rank, soma);
        }
        MPI_Send(&soma, 1, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD);

    } else {
        for (int i = world_rank * m + 1; i <= (world_rank + 1) * m; i++)
            soma += funcao(n, v, i * h);
        soma = soma * h;
        if (debug) {
            printf("soma(%d) = %20.12f", world_rank, soma);
        }
        MPI_Send(&soma, 1, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD);
    }


}

int main() {

    MPI_Init(NULL, NULL); // Inicializa o MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double C[MM];
    C[0] = 0.5 * integra(0, v, NN);
    for (int i = 1; i <= MM; i++)
        C[i] = integra(i, v, NN);
    if (world_rank == 0) {
        for (int i = 0; i <= MM; i++)
            printf("\nValor da integral em  C%d  = %20.17es", i, C[i]);
        double u_x_t;
        for (int i=1;i<10;i++) {
            u_x_t = 2 * M_PI * v * soma1(MM, v,(double)(i)*0.1, 0.5, C) / soma2(MM, v, (double)(i)*0.1, 0.5, C);
            printf("\nValor de U(x=%3.1f; t=0.5) = %7.4es", i*0.1,u_x_t);
        }
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
