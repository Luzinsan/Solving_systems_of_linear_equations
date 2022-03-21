#include <iostream>
#include "InputData.h"
using namespace luMath;

int main()
{
    setlocale(LC_ALL, "Rus");
    InputData data;
    data.setGaussMethod();
    data.setInverseMatrixByMethod(InputData::GaussMethod);
    //data.DecompositionMethod();
    //data.OrtogonalizationMethod();
    //data.SimpleIterationMethod();
    return 0;
}
