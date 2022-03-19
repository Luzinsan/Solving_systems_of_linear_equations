#ifndef INPUTDATA_H
#define INPUTDATA_H
#include <iostream>
#include <fstream>
#include <cmath>
#include "Matrix.h"
#include "Vector.h"

namespace luMath
{
    char getSymbol(std::initializer_list<char> list,
        std::string notification_message = "",
        std::string error_message = "Недопустимое значение, попробуйте ещё раз.\n->");

    double getDouble(double min = -DBL_MAX,
        double max = DBL_MAX,
        std::string notification_message = "",
        std::string error_message = "Недопустимое значение, попробуйте ещё раз.\n->");

    class InputData
    {
    private:
        std::ifstream* _fin;
        Matrix<double>* _matrix;
        Vector<double>* _result;
        int _type;
        int _matrixDimention;
        int _NAfterComma;
    public:

        InputData()
        {
            _fin = new std::ifstream("input.txt");

            *_fin >> _type;
            std::cout << "\n\tТип задачи: " << _type;
            *_fin >> _matrixDimention;
            std::cout << "\n\tПорядок матрицы: " << _matrixDimention;

            double* _array = new double[_matrixDimention * (_matrixDimention + 1)];
            for (int i = 0; i < (_matrixDimention + 1) *_matrixDimention; i++)
                *_fin >> _array[i];

            _matrix = new Matrix<double>(_matrixDimention, _matrixDimention + 1, _array);
            delete[] _array;
            std::cout << "\n\tКоэффициенты считанной матрицы и вектор свободных коэффициентов (последний стоблик):\n\n" << *_matrix;
            delete _fin;
        }
        ~InputData()
        {
            delete[] _matrix;
            delete[] _result;
        }
  

        void GaussMethod() 
        {
            Matrix<double> tempMatrix(*_matrix);
            // Прямой ход метода Гаусса - преобразование матрицы к треугольному виду
            for (int i = 0; i < _matrixDimention; i++) // проходим по всем строкам
            {
                double coeff = tempMatrix[i][i]; // запоминаем коэффициент по диагонали
                for (int j = i; j < _matrixDimention + 1; j++) // проходим по всем элементам текущей строки, включая вектор коэффициентов
                {
                    std::cout << '\n' << tempMatrix[i][j] << " / " << coeff << " = ";
                    tempMatrix[i][j] /= coeff; // преобразуем диагональ матрицы к единичной диагонали (подготовительная стадия)
                    std::cout << tempMatrix[i][j] << '\n';

                    std::cout << '\n' << std::setw(10) << tempMatrix;

                }
                std::cout << '\n' << std::setw(10) << tempMatrix;
                for (int j = i+1; j < _matrixDimention; j++) 
                {
                    coeff = tempMatrix[j][i]; // запоминаем коэффициент преобразования 
                    for (int k = i; k < _matrixDimention+1; k++) // проходим по всем элементам текущей строки, включая вектор коэффициентов
                    {
                        std::cout << '\n' << tempMatrix[j][k] << " - " << coeff << " * " << tempMatrix[i][k] << " = ";
                        tempMatrix[j][k] -= coeff * tempMatrix[i][k]; // преобразуем диагональ матрицы к единичной диагонали (подготовительная стадия)
                        std::cout << tempMatrix[j][k] << '\n';

                        std::cout << '\n' << std::setw(10) << tempMatrix;

                    }
                    std::cout << '\n' << std::setw(10) << tempMatrix;
                }
               
            }
            std::cout << "\nМатрица с единичной диагональю: \n" << std::setw(10) << tempMatrix;
                
            // Обратный ход метода Гаусса
            _result = new Vector<double>(_matrixDimention);
            for (int i = _matrixDimention - 1; i >= 0; i--)
            {
                double sumCoeff = 0;
                for (int j = i + 1; j < _matrixDimention; j++) 
                {
                    std::cout << "\n" << tempMatrix[i][j] << " * " << (*_result)[j] << " + " << sumCoeff << " = ";
                    sumCoeff += tempMatrix[i][j] * (*_result)[j];
                    std::cout << sumCoeff << "\n";
                }
                std::cout << "result: " << tempMatrix[i][_matrixDimention] << " - " << sumCoeff << " = ";
                (*_result)[i] = tempMatrix[i][_matrixDimention] - sumCoeff;
                std::cout << (*_result)[i] << "\n";

            }


            (*_result).transposition();

            std::cout << *_result;


        
        }
    
        void DecompositionMethod() 
        {
            // Раскладываем матрицу _matrix на матрицы B и C так, что A = B*C
            Matrix<double> B(_matrixDimention), C(_matrixDimention);
            for(int j = 0; j < _matrixDimention; j++)
                for (int i = j; j < _matrixDimention; j++) 
                {
                    double sumCoeff = 0;
                    for (int k = 1; k < j - 1; k++)
                        sumCoeff += B[i][k] * C[k][j];

                    B[i][j] = (*_matrix)[i][j] - sumCoeff;
                    
                }
        
        }
    };
}
#endif