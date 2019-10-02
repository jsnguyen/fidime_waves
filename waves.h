// waves.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <fstream>

#define PI 3.14159265359

double WaveEquationExplicitInitial(double f[3], double eta, double timestep, double init_v);
double WaveEquationExplicit(double f[3], double eta, double prev);
double GaussianWavePacket(double x, double start, double c, double t, double sigma);
double GaussianWavePacketDerivative(double x, double start, double c, double t, double sigma);