#pragma once
#define MatrixElement_H
#include <vector>
#include "Eigen/Core"
using namespace Eigen;
class MatrixElement
{

private:
    const double pi = 4.0 * atan(1.0);
    double c12_inv;
    double c23_inv;
    double c13_inv;

public:
    MatrixElement();
    double overlap_single(Vector3d state1, double b11, double b12, double b22);
    double energy_single(Vector3d state1, double b11, double b12, double b22);
    double overlap(Vector3d state1, double b11, double b12, double b22);
    double energy(Vector3d state1, double b11, double b12, double b22);
};