/*
cpucycles amd64cpuinfo.h version 20060318
D. J. Bernstein
Public domain.
*/

#ifndef CPUCYCLES_amd64cpuinfo_h
#define CPUCYCLES_amd64cpuinfo_h

long long cpucycles_amd64cpuinfo(void);
long long cpucycles_amd64cpuinfo_persecond(void);

#ifndef cpucycles_implementation
#define cpucycles_implementation "amd64cpuinfo"
#define cpucycles cpucycles_amd64cpuinfo
#define cpucycles_persecond cpucycles_amd64cpuinfo_persecond
#endif

#endif
