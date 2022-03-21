#include <iostream>
#include "InputData.h"
using namespace luMath;

int main()
{
    setlocale(LC_ALL, "Rus");
    InputData<double> data;
    data.setGaussMethod();
    data.setInverseMatrixByMethod(InputData<double>::GaussMethod);
    std::cout <<"\nПроверка: A * A' = \n" <<std::fixed <<std::setprecision(2) <<std::setw(10) << data.getInverseMatrix() * data.getMainMatrix();
    //data.DecompositionMethod();
    //data.OrtogonalizationMethod();
    //data.SimpleIterationMethod();
    return 0;
}
