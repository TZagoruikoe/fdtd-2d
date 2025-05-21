#include "fdtd-2d.h"
#include <openacc.h>
double bench_t_start, bench_t_end;

#define EPSILON 1e-6

static
double rtclock() {
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
        printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start() {
    bench_t_start = rtclock ();
}

void bench_timer_stop() {
    bench_t_end = rtclock ();
}

void bench_timer_print() {
    printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
}

int compare_array_with_file(const char* filename, int nx, int ny, double arr[nx][ny]) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open file %s for reading.\n", filename);
        return -1;
    }

    double val;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (fscanf(f, "%lf", &val) != 1) {
                fprintf(stderr, "ERROR: Failed to read value at [%d][%d]\n", i, j);
                fclose(f);
                return -1;
            }
            if (fabs(val - arr[i][j]) > EPSILON) {
                fprintf(stderr, "ERROR: Mismatch at [%d][%d]: expected %.10lf, got %.10lf\n", i, j, val, arr[i][j]);
                fclose(f);
                return 1;
            }
        }
    }

    fclose(f);
    return 0;
}

static
void init_array (int tmax,
                 int nx,
                 int ny,
                 double ex[nx][ny],
                 double ey[nx][ny],
                 double hz[nx][ny],
                 double _fict_[tmax]) {
    int i, j;
  

#pragma acc parallel loop copyout(_fict_[0:tmax])
    for (i = 0; i < tmax; i++)
        _fict_[i] = (double) i;

#pragma acc parallel loop copyout(ex[0:nx][0:ny], ey[0:nx][0:ny], hz[0:nx][0:ny]) collapse(2)
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++) {
            ex[i][j] = ((double) i * (j + 1)) / nx;
            ey[i][j] = ((double) i * (j + 2)) / ny;
            hz[i][j] = ((double) i * (j + 3)) / nx;
        }
}

static
void kernel_fdtd_2d(int tmax,
                    int nx,
                    int ny,
                    double ex[nx][ny],
                    double ey[nx][ny],
                    double hz[nx][ny],
                    double _fict_[tmax]) {
    int t, i, j;

#pragma acc data copy(ex[0:nx][0:ny], ey[0:nx][0:ny], hz[0:nx][0:ny], _fict_[0:tmax])
{
    for(t = 0; t < tmax; t++) {
    #pragma acc parallel loop present(ey, _fict_) independent
        for (j = 0; j < ny; j++)
            ey[0][j] = _fict_[t];

    #pragma acc parallel loop present(ex, ey, hz) collapse(2)
        for (i = 1; i < nx; i++)
            for (j = 1; j < ny; j++) {
                ey[i][j] = ey[i][j] - 0.5 * (hz[i][j] - hz[i - 1][j]);
                ex[i][j] = ex[i][j] - 0.5 * (hz[i][j] - hz[i][j - 1]);
            }

    #pragma acc parallel loop present(ey, hz) independent
        for (i = 1; i < nx; i++)
            ey[i][0] = ey[i][0] - 0.5 * (hz[i][0] - hz[i - 1][0]);

    #pragma acc parallel loop present(ex, hz) independent
        for (j = 1; j < ny; j++)
            ex[0][j] = ex[0][j] - 0.5 * (hz[0][j] - hz[0][j - 1]);

    #pragma acc parallel loop present(ex, ey, hz) collapse(2)
        for (i = 0; i < nx - 1; i++)
            for (j = 0; j < ny - 1; j++)
                hz[i][j] = hz[i][j] - 0.7* (ex[i][j + 1] - ex[i][j] + ey[i + 1][j] - ey[i][j]);
    }
}
}

int main(int argc, char** argv) {
    int tmax = TMAX;
    int nx = NX;
    int ny = NY;
    double (*ex)[nx][ny];
    ex = (double(*)[nx][ny])malloc ((nx) * (ny) * sizeof(double));
  
    double (*ey)[nx][ny];
    ey = (double(*)[nx][ny])malloc ((nx) * (ny) * sizeof(double));
  
    double (*hz)[nx][ny];
    hz = (double(*)[nx][ny])malloc ((nx) * (ny) * sizeof(double));
  
    double (*_fict_)[tmax];
    _fict_ = (double(*)[tmax])malloc ((tmax) * sizeof(double));
  
    init_array (tmax, nx, ny, *ex, *ey, *hz, *_fict_);
  
    bench_timer_start();
  
    kernel_fdtd_2d (tmax, nx, ny, *ex, *ey, *hz, *_fict_);
  
    bench_timer_stop();

    int ok_ex = compare_array_with_file("ex_medium.txt", nx, ny, *ex);
    int ok_ey = compare_array_with_file("ey_medium.txt", nx, ny, *ey);
    int ok_hz = compare_array_with_file("hz_medium.txt", nx, ny, *hz);

    if (ok_ex == 0 && ok_ey == 0 && ok_hz == 0) {
        printf("OK\n");
    } else {
        printf("ERROR: Output mismatch detected.\n");
    }

    bench_timer_print();
  
    free((void*)ex);
    free((void*)ey);
    free((void*)hz);
    free((void*)_fict_);
  
    return 0;
}