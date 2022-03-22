#include <iostream>
#include "InputData.h"
using namespace luMath;

int main()
{
    setlocale(LC_ALL, "Rus");
    InputData<double> data;

    switch (data.getMethod()) 
    {
    case InputData<double>::METHOD::GAUSS:
        switch (data.getTask())
        {
        case InputData<double>::TASK::ROOT:
            data.getRoot(InputData<double>::GaussMethod);
            break;
        case InputData<double>::TASK::DETERMINANT:
            data.getDeterminant();
            break;
        case InputData<double>::TASK::INVERS:
            data.setInverseMatrixByMethod(InputData<double>::GaussMethod);
            break;
        }
        break;
    case InputData<double>::METHOD::DECOMPOSOTION:
        switch (data.getTask())
        {
        case InputData<double>::TASK::ROOT:
            data.getRoot(InputData<double>::DecompositionMethod);
            break;
        case InputData<double>::TASK::DETERMINANT:
            data.getDeterminant();
            break;
        case InputData<double>::TASK::INVERS:
            data.setInverseMatrixByMethod(InputData<double>::DecompositionMethod);
            break;
        }
        break;
    }
    
   
    
    
    //data.OrtogonalizationMethod();
    //data.SimpleIterationMethod();
    return 0;
}
