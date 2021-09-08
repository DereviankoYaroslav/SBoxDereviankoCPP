//
// Created by Ярослав on 08.09.2021.
//

#ifndef SBOXDEREVIANKOCPP_MAIN_H
#define SBOXDEREVIANKOCPP_MAIN_H

#include <stdlib.h>
#include <stdint.h>

void particleSwarmOptimization(void* ddata);

typedef struct dereviankoData {
    int size;
    int count;
    int N;
    int maxIter;
    int mode;
    int* sboxes;
}dereviankoData;

int raiseToPower(int num, int pow);

int *valueToBinary(int i, int rank);

int *elemsForN(int size);

int *binaryElementsApprox(int *arr, int size, int count);

int *SBoxApprox(int *sbox, int size, int count);

int LATMax(int *sbox, int size, int count);

int *SBoxGeneratingDec(int n, int m, int counter);

void FisherYates(int *player, int n);

int myModulusDec(int number, int mod);

#endif

