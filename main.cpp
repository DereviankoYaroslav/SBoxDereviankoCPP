#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#include <pthread.h>
#include <math.h>
#include <conio.h>
#include <thread>
#include "main.h"

CONST int aesSbox[] = {99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125,
                       250,
                       89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52,
                       165,
                       229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117,
                       9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252,
                       177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2,
                       127,
                       80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205,
                       12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144,
                       136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145,
                       149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186,
                       120,
                       37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246,
                       14,
                       97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233,
                       206,
                       85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22};

int main() {
    int mode;
    int N;
    int maxIter;
    int threadsCount;

    FILE *readThreads;
    readThreads = fopen("Threads count.txt", "r");
    fscanf(readThreads, "%d", &threadsCount);
    fclose(readThreads);

    FILE *readMode;
    readMode = fopen("Mode.txt", "r");
    fscanf(readMode, "%d", &mode);
    fclose(readMode);

    FILE *readN;
    readN = fopen("N value (amount of s-boxes in population).txt", "r");
    fscanf(readN, "%d", &N);
    fclose(readN);

    FILE *readMaxIter;
    readMaxIter = fopen("maxIter value (number of iterations).txt", "r");
    fscanf(readMaxIter, "%d", &maxIter);
    fclose(readMaxIter);

    int *sboxes = (int *) calloc(N * 256, sizeof(int));

    dereviankoData data = {
            256,
            8,
            N,
            maxIter,
            mode,
            sboxes
    };
    dereviankoData datas[threadsCount];
    for (int i = 0; i < threadsCount; i++){
        datas[i] = data;
    }

    std::thread threads[threadsCount];
    auto start = std::chrono::system_clock::now();
    printf("\nPSO in progress...\n");
    for (int i = 0; i < threadsCount; i++)
        threads[i] = std::thread(particleSwarmOptimization, &datas[i]);

    for (int i = 0; i < threadsCount; i++)
        threads[i].join();

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    printf("Часу витрачено: %f (секунди) (%f (хвилини))\n", diff.count(), diff.count() / 60.0);

    FILE *file;
    fopen_s(&file, "PSO results.txt", "w");
    if (file == NULL) {
        printf("ERROR: Can't save sbox to file!\n");
        for (;;);
    }
    fprintf(file, "\n\nPSO S-boxes\n\n");
    for (int i = 0; i < threadsCount; i++) {
        for (int q = 1; q < datas[i].N; ++q) {
            for (int w = 0; w < datas[i].size; ++w) {
                fprintf(file, "%d, ", datas[i].sboxes[q * datas[i].size + w]);
            }
            fprintf(file, "\n\n");
        }
        fprintf(file, "\n");

        printf("\n\nPSO S-boxes\n\n");

        for (int q = 1; q < datas[i].N; ++q) {
            for (int w = 0; w < datas[i].size; ++w) {
                printf("%d, ", datas[i].sboxes[q * datas[i].size + w]);
            }
            printf("\n\n");
        }
    }

    system("PAUSE");

    //int *ar = particleSwarmOptimization(256,8,N,maxIter,mode);
}

int raiseToPower(int num, int pow) {
    int res = 1;
    for (int i = 0; i < pow; ++i) {
        res *= num;
    }
    return res;
}

void FisherYates(int *arr, int n) {
    int i, j, tmp;

    for (i = n - 1; i > 0; i--) {
        j = rand() % (i + 1);
        tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
}

int *SBoxGeneratingDec(int n, int m, int counter) {
    int size = raiseToPower(2, n);
    int *dec = (int *) calloc(size, sizeof(int));
    srand((counter * counter) % size);
    for (int i = 0; i < size;) {
        dec[i] = rand() % size;
        int contains = 0;
        for (int j = 0; j < i; ++j) {
            if (dec[i] == dec[j]) {
                contains = 1;
                break;
            }
        }
        if (!contains) {
            i++;
        }
    }
    /*printf("Generated s-box: ");
    for (int i = 0; i < size; ++i) {
        printf("%d, ", dec[i]);
    }
    printf("\n");*/
    FisherYates(dec, size);
    return dec;
}

int *valueToBinary(int i, int rank) {
    int *res = (int *) calloc(rank, sizeof(int));
    for (int j = 0; j < rank; ++j) {
        res[rank - 1 - j] = i >> j & 1;
    }
    return res;
}

int *binaryElementsApprox(int *arr, int size, int count) {
    int *result = (int *) calloc(size * count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(arr[i], count);
        for (int j = count - 1; j >= 0; j--) {
            result[j * size + i] = bin[j];
        }
        free(bin);
    }
    return result;
}

int *SBoxApprox(int *sbox, int size, int count) {
    int *result = binaryElementsApprox(sbox, size, count);
    return result;
}

int *elemsForN(int size) {
    int *result = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        result[i] = i;
    }
    return result;
}

int LATMax(int *sbox, int size, int count) {
    int *ar = SBoxApprox(sbox, size, count);
    int *elems = elemsForN(size);
    int *binelems = binaryElementsApprox(elems, size, count);
    int *temp = (int *) calloc(size, sizeof(int));
    int *temp2 = (int *) calloc(size, sizeof(int));
    int *coefficients = (int *) calloc(size * size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin1 = valueToBinary(i, count);
        for (int k = count - 1; k >= 0; k--) {
            if (bin1[k]) {
                //printf("K===%d ",k);
                //printf("X == \n ");
                for (int l = 0; l < size; ++l) {
                    temp[l] = temp[l] ^ binelems[k * size + l];
                    //printf("%d ",temp[l]);
                }
                //printf("\n ");
            }
        }
        //printf("\n ");
        for (int j = 0; j < size; ++j) {
            int *bin2 = valueToBinary(j, count);
            for (int q = count - 1; q >= 0; q--) {
                if (bin2[q]) {
                    //printf("K===%d ",k);
                    //printf("\nY [%d]== \n ", j);
                    for (int w = 0; w < size; ++w) {
                        temp2[w] = temp2[w] ^ ar[q * size + w];
                        //printf("%d ",temp2[l]);
                    }
                }
            }
            //printf("\n ");
            int calc = 0;
            for (int r = 0; r < size; ++r) {
                temp2[r] = temp2[r] ^ temp[r];
                //printf("%d ", temp2[l]);
                if (temp2[r] == 0) {
                    ++calc;
                }
                temp2[r] = 0;
            }
            int result = 0;
            result = calc - (size / 2);
            //printf("COEFFS = %d ", result);
            coefficients[i * size + j] = result;
            free(bin2);
        }
        for (int t = 0; t < size; ++t) {
            temp[t] = 0;
        }
        //printf("\n ");
        free(bin1);
    }
    for (int n = 0; n < size; ++n) {
        for (int m = 0; m < size; ++m) {
            //printf("%d ", coefficients[n*size+m]);
        }
        //printf("\n");
    }
    int result = 0;
    for (int p = 1; p < size * size; p++) {
        if (abs(coefficients[p]) > result)
            result = abs(coefficients[p]);
    }
    free(ar);
    free(elems);
    free(binelems);
    free(temp);
    free(temp2);
    free(coefficients);
    return result;
}

int myModulusDec(int number, int mod) {
    if (number < 0) {
        while (number < 0) {
            number = number + mod;
        }
    }
    return number % mod;
}

void particleSwarmOptimization(void *ddata) {
    dereviankoData *dd = (dereviankoData *) ddata;

    srand(time(NULL));
    int flag = rand() % dd->size;
    int population[2 * dd->N][dd->size];
    for (int i = 0; i < dd->size; ++i) {
        population[0][i] = aesSbox[i];
    }
    for (int q = 1; q < dd->N; ++q) {
        int *ar1 = SBoxGeneratingDec(dd->count, dd->count, q + flag);
        for (int w = 0; w < dd->size; ++w) {
            population[q][w] = ar1[w];
        }
        free(ar1);
    }
    int arrNL[dd->N];
    for (int q = 0; q < dd->N; ++q) {
        for (int w = 0; w < dd->size; ++w) {
            printf("%d ", population[q][w]);
        }
        int LAT = LATMax(population[q], dd->size, dd->count);
        int NL = raiseToPower(2, dd->count - 1) - LAT;
        printf("\nNon-linearity from LAT = %d \n", NL);
        printf("\n");
        arrNL[q] = NL;
    }
    int g[dd->size];
    for (int i = 0; i < dd->N; ++i) {
        for (int j = dd->N - 1; j > i; --j) {
            if (arrNL[j] > arrNL[j - 1]) {
                int t = arrNL[j - 1];
                arrNL[j - 1] = arrNL[j];
                arrNL[j] = t;
                for (int k = 0; k < dd->size; ++k) {
                    g[k] = population[j - 1][k];
                    population[j - 1][k] = population[j][k];
                    population[j][k] = g[k];
                }
            }
        }
    }
    printf("\n");
    for (int q = 0; q < dd->N; ++q) {
        printf("\n%d ", arrNL[q]);
    }
    printf("\nSORTED BY Non-Linearity\n");
    for (int q = 0; q < dd->N; ++q) {
        for (int w = 0; w < dd->size; ++w) {
            printf("%d, ", population[q][w]);
        }
        printf("\n\n");
    }
    int gBest[dd->size];
    for (int m = 0; m < dd->size; ++m) {
        gBest[m] = population[0][m];
    }
    int pBest[dd->N][dd->size];
    for (int i = 1; i < dd->N; ++i) {
        for (int j = 0; j < dd->size; ++j) {
            pBest[i - 1][j] = population[i][j];
        }
    }
    printf("\n\n");
    printf("\ngBest\n");
    for (int m = 0; m < dd->size; ++m) {
        printf("%d, ", gBest[m]);
    }
    printf("\npBest\n");
    for (int q = 0; q < dd->N - 1; ++q) {
        for (int w = 0; w < dd->size; ++w) {
            printf("%d, ", pBest[q][w]);
        }
        printf("\n\n");
    }
    double weight1 = 0.1;
    double weight2 = 1.6;
    double weightCur;
    int curIter = 0;
    int Vel[dd->N][dd->size];
    int arrNLSorted[dd->size];
    while (dd->maxIter > 0) {
        weightCur = weight1 + (curIter - 1) * ((weight2 - weight1) / dd->maxIter);
        printf("\nWHILE\n");
        int Q = 100;
        int rd1 = rand() % (Q);
        double xr1 = (double) rd1 / Q;
        double c1 = 2 * xr1;
        int rd2 = rand() % (Q);
        double xr2 = (double) rd2 / Q;
        double c2 = 2 * xr2;
        int rd3 = rand() % (Q);
        double xr3 = (double) rd3 / Q;
        double r1 = xr3;
        int rd4 = rand() % (Q);
        double xr4 = (double) rd4 / Q;
        double r2 = xr4;
        /*printf("xr1 = %lf ", xr1);
        printf("xr2 = %lf ", xr2);
        printf("xr3 = %lf ", xr3);
        printf("xr4 = %lf \n", xr4);
        printf("\n\n");*/
        for (int b = 0; b < dd->N; ++b) {
            arrNLSorted[b] = arrNL[b];
        }
        int tempSbox[dd->size];
        int tempSbox2[dd->size];
        for (int i = 0; i < dd->N; ++i) {
            for (int j = 0; j < dd->size;) {
                if (dd->mode == 1) {
                    Vel[i][j] = ceil(weightCur * Vel[i][j] + c1 * r1 * (gBest[j] - population[i][j] +
                                                                        c2 * r2 * (gBest[j] - population[i][j])));
                }
                if (dd->mode == 0) {
                    Vel[i][j] = ceil(weightCur * Vel[i][j] + c1 * r1 * (pBest[i][j] - population[i][j] +
                                                                        c2 * r2 * (gBest[j] - population[i][j])));
                }
                if (Vel[i][j] < 0) {
                    Vel[i][j] = myModulusDec((Vel[i][j] + 256), 256);
                }
                //printf("Vel[%d][%d] = %d ", i,j,Vel[i][j]);
                int X = myModulusDec((population[i][j] + Vel[i][j]), 256);
                int contains;
                if (contains == 0) {
                    tempSbox[j] = X;
                } else {
                    tempSbox[j] = myModulusDec((tempSbox[j] + rand()), 256);
                    /*if (tempSbox[j] == 0){
                        tempSbox[j] = myModulusDec((tempSbox[j]+rand()),256);
                    }*/
                };
                //printf("\n%d temp box b4 cycle [%d]", tempSbox[j],j);
                contains = 0;
                for (int k = 0; k < j; ++k) {
                    if (tempSbox[k] == tempSbox[j]) {
                        //printf("\n%d k = \n", k);
                        //printf("\n%d j = \n", j);
                        //tempSbox[j] = myModulusDec((tempSbox[j]+rand()),256);
                        //tempSbox2[j] = tempSbox[j];
                        //printf("\n%d j", tempSbox[j]);
                        //printf("\n%d k", tempSbox[k]);
                        contains = 1;
                        break;
                    }
                }
                if (!contains) {
                    j++;
                }
            }
            for (int k = 0; k < dd->size; ++k) {
                population[dd->N + i][k] = tempSbox[k];
                //printf("%d ", tempSbox[k]);
            }
            //printf("\n");
        }
        /*printf("\nNEW Arrays\n");
        for (int q = 0; q < 2*N; ++q){
            for(int w = 0; w < size; ++w){
                printf("%d, ",population[q][w]);
            }
            int LAT = LATMax(population[q],size,count);
            int NL = raiseToPower(2, count - 1) - LAT;
            printf( "\nNon-linearity from LAT = %d \n", NL);
            printf("\n");
            printf("\n\n");
        }*/
        int arrNL2[2 * dd->N];
        for (int q = 0; q < 2 * dd->N; ++q) {
            for (int w = 0; w < dd->size; ++w) {
                printf("%d ", population[q][w]);
            }
            int LAT = LATMax(population[q], dd->size, dd->count);
            int NL = raiseToPower(2, dd->count - 1) - LAT;
            printf("\nNon-linearity from LAT = %d \n", NL);
            printf("\n");
            arrNL2[q] = NL;
        }
        int g2[dd->size];
        for (int i = 0; i < 2 * dd->N; ++i) {
            for (int j = (2 * dd->N - 1); j > i; --j) {
                if (arrNL2[j] > arrNL2[j - 1]) {
                    int h = arrNL2[j - 1];
                    arrNL2[j - 1] = arrNL2[j];
                    arrNL2[j] = h;
                    for (int k = 0; k < dd->size; ++k) {
                        g2[k] = population[j - 1][k];
                        population[j - 1][k] = population[j][k];
                        population[j][k] = g2[k];
                    }
                }
            }
        }
        printf("\n");
        for (int q = 0; q < 2 * dd->N; ++q) {
            printf("\n%d ", arrNL2[q]);
        }
        for (int q = dd->N; q < 2 * dd->N; ++q) {
            for (int w = 0; w < dd->size; ++w) {
                population[q][w] = 0;
            }
        }
        printf("\nSORTED BY Non-Linearity\n");
        for (int q = 0; q < dd->N; ++q) {
            for (int w = 0; w < dd->size; ++w) {
                printf("%d, ", population[q][w]);
            }
            printf("\n\n");
        }
        for (int m = 0; m < dd->size; ++m) {
            gBest[m] = population[0][m];
        }
        for (int i = 1; i < dd->N; ++i) {
            for (int j = 0; j < dd->size; ++j) {
                pBest[i - 1][j] = population[i][j];
            }
        }
        printf("\n\n");
        printf("\ngBest\n");
        for (int m = 0; m < dd->size; ++m) {
            printf("%d, ", gBest[m]);
        }
        printf("\npBest\n");
        for (int q = 0; q < dd->N - 1; ++q) {
            for (int w = 0; w < dd->size; ++w) {
                printf("%d, ", pBest[q][w]);
            }
            printf("\n\n");
        }
        dd->maxIter = dd->maxIter - 1;
        dd->mode = 0;
    }
    //printf("\n\nFinal data\n\n");
    for (int q = 0; q < dd->N; ++q) {
        for (int w = 0; w < dd->size; ++w) {
            dd->sboxes[q * dd->size + w] = population[q][w];
            //printf("%d, ", result[q * size + w]);
        }
        //printf("\n\n");
    }
    //return result;
}