#include <iostream>
#include "InputData.h"
using namespace luMath;

int main()
{
    setlocale(LC_ALL, "Rus");
    InputData<double> data;
    data.setMethod(InputData<double>::GaussMethod);
    data.setInverseMatrixByMethod(InputData<double>::GaussMethod);
    std::cout <<"\nПроверка: A * A' = \n" <<std::fixed <<std::setprecision(2) <<std::setw(10) << data.getInverseMatrix() * data.getMainMatrix();
   
    data.setMethod(InputData<double>::DecompositionMethod);
    data.setInverseMatrixByMethod(InputData<double>::DecompositionMethod);
    std::cout << "\nПроверка: A * A' = \n" << std::fixed << std::setprecision(2) << std::setw(10) << data.getInverseMatrix() * data.getMainMatrix();
    
    //data.OrtogonalizationMethod();
    //data.SimpleIterationMethod();
    return 0;
}
