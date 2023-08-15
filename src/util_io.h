#ifndef UTILIO_H
#define UTILIO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <complex>
#include <getline.h>

using namespace std;

void dump_1d_real(char* filename, int n, double* arr) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%+6.5f\n", arr[i]);
    }

    fclose(fp);
}

void dump_1d_cplx(char* filename, int n, complex<double>* arr) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%+6.5f %+6.5f\n", real(arr[i]), imag(arr[i]));
    }

    fclose(fp);
}

void dump_2d_real(char* filename, int rows, int cols, double* arr) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            fprintf(fp, "%+6.5f ", arr[i * cols + j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

// using numpy format: matrix of (<real>+<cplx>j)
void dump_2d_cplx(char* filename, int rows, int cols, complex<double>* arr) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    int index;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            index = i * cols + j;

            if (imag(arr[index]) < 0) {
                fprintf(fp, "(%+6.5f%+6.5fj) ", real(arr[index]), imag(arr[index]));
            }
            else {
                fprintf(fp, "(%+6.5f+%6.5fj) ", real(arr[index]), imag(arr[index]));
            }
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

int read_lcfs_header(char* filename, int* mpol, int* ntor, int* mnmax, int* nfp) {
    char* line = NULL;
    char* mpol_line = "# mpol\n";
    char* ntor_line = "# ntor\n";
    char* mnmax_line = "# mnmax\n";
    char* nfp_line = "# nfp\n";

    int found_mpol = 0;
    int found_ntor = 0;
    int found_mnmax = 0;
    int found_nfp = 0;

    size_t len = 0;
    ssize_t read;

    int status = 0;

    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: could not open %st\n", filename);
        status = 1;
        return status;
    }

    int keep_going = 1;

    while (keep_going && (read = getline(&line, &len, fp)) != -1) {

        if (!found_mpol && strcmp(line, mpol_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*mpol) = atoi(line);
                found_mpol = 1;
            }
            else {
                printf("failed to read line after '%s'\n", mpol_line);
                status = 1;
                keep_going = 0;
            }
        }
        else if (!found_ntor && strcmp(line, ntor_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*ntor) = atoi(line);
                found_ntor = 1;
            }
            else {
                printf("failed to read line after '%s'\n", ntor_line);
                status = 1;
                keep_going = 0;
            }
        }
        else if (!found_mnmax && strcmp(line, mnmax_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*mnmax) = atoi(line);
                found_mnmax = 1;
            }
            else {
                printf("failed to read line after '%s'\n", mnmax_line);
                status = 1;
                keep_going = 0;
            }
        }
        else if (!found_nfp && strcmp(line, nfp_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*nfp) = atoi(line);
                found_nfp = 1;
            }
            else {
                printf("failed to read line after '%s'\n", nfp_line);
                status = 1;
                keep_going = 0;
            }
        }

        if (found_mpol && found_ntor && found_mnmax && found_nfp) {
            keep_going = 0;
        }
    }

    fclose(fp);

    if (line) {
        free(line);
    }

    return status;
}

int read_lcfs(char* filename, int mnmax, double* rmnc, double* zmns) {
    char* line = NULL;
    char* rmnc_line = "# rmnc\n";
    char* zmns_line = "# zmns\n";

    int found_rmnc = 0;
    int found_zmns = 0;

    size_t len = 0;
    ssize_t read;

    int status = 0;

    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: could not open %s\n", filename);
        status = 1;
        return status;
    }

    int keep_going = 1;

    while (keep_going && (read = getline(&line, &len, fp)) != -1) {
        if (!found_rmnc && strcmp(line, rmnc_line) == 0) {
            for (int i = 0; keep_going && i < mnmax; ++i) {
                read = getline(&line, &len, fp);
                if (read != -1) {
                    rmnc[i] = atof(line);
                }
                else {
                    printf("failed to read line %d after '%s'\n", i, rmnc_line);
                    status = 1;
                    keep_going = 0;
                }
            }
            found_rmnc = 1;
        }
        else if (!found_zmns && strcmp(line, zmns_line) == 0) {
            for (int i = 0; keep_going && i < mnmax; ++i) {
                read = getline(&line, &len, fp);
                if (read != -1) {
                    zmns[i] = atof(line);
                }
                else {
                    printf("failed to read line %d after '%s'\n", i, zmns_line);
                    status = 1;
                    keep_going = 0;
                }
            }
            found_zmns = 1;
        }

        if (found_rmnc && found_zmns) {
            keep_going = 0;
        }
    }

    fclose(fp);

    if (line) {
        free(line);
    }

    return status;
}

#endif
