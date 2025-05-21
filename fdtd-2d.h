#ifndef _FDTD_2D_H
#define _FDTD_2D_H 
#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRA_DATASET) && !defined(MEGA_DATASET)
#define LARGE_DATASET
#endif
#if !defined(TMAX) && !defined(NX) && !defined(NY)
#ifdef MINI_DATASET
#define TMAX 20
#define NX 20
#define NY 30
#endif
#ifdef SMALL_DATASET
#define TMAX 40
#define NX 60
#define NY 80
#endif
#ifdef MEDIUM_DATASET
#define TMAX 100
#define NX 500
#define NY 600
#endif
#ifdef LARGE_DATASET
#define TMAX 100
#define NX 2000
#define NY 2400
#endif
#ifdef EXTRA_DATASET
#define TMAX 100
#define NX 8000
#define NY 9600
#endif
#ifdef MEGA_DATASET
#define TMAX 100
#define NX 12000
#define NY 14400
#endif
#endif
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#endif