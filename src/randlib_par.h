/* randlib_par.h, part of the Global Epidemic Simulation v1.0 BETA
/* Header for parallel random library
/*
/* Copyright 2012, MRC Centre for Outbreak Analysis and Modelling
/* 
/* Licensed under the Apache License, Version 2.0 (the "License");
/* you may not use this file except in compliance with the License.
/* You may obtain a copy of the License at
/*
/*       http://www.apache.org/licenses/LICENSE-2.0
/*
/* Unless required by applicable law or agreed to in writing, software
/* distributed under the License is distributed on an "AS IS" BASIS,
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
/* See the License for the specific language governing permissions and
/* limitations under the License.
*/

#ifndef RNDLIB_PAR_H
#define RNDLIB_PAR_H

long ignbin(long,double);
long ignpoi(double);
long ignbin_mt(long ,double ,int);
long ignpoi_mt(double,int);
double ranf(void);
double ranf_mt(int);
void initSeeds(long ,long);
double sexpo_mt(int);
double sexpo(void);
long mltmod(long ,long ,long );
double snorm(void);
double snorm_mt(int);
double fsign(double, double );
double gengam(double,double);
double gengam_mt(double,double,int);
double sgamma(double);
double sgamma_mt(double,int);

#endif
