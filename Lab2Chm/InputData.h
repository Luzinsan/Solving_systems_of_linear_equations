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
        std::string error_message = "������������ ��������, ���������� ��� ���.\n->");

    double getDouble(double min = -DBL_MAX,
        double max = DBL_MAX,
        std::string notification_message = "",
        std::string error_message = "������������ ��������, ���������� ��� ���.\n->");

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
            std::cout << "\n\t��� ������: " << _type;
            *_fin >> _matrixDimention;
            std::cout << "\n\t������� �������: " << _matrixDimention;

            double* _array = new double[_matrixDimention * (_matrixDimention + 1)];
            for (int i = 0; i < (_matrixDimention + 1) *_matrixDimention; i++)
                *_fin >> _array[i];

            _matrix = new Matrix<double>(_matrixDimention, _matrixDimention + 1, _array);
            delete[] _array;
            std::cout << "\n\t������������ ��������� ������� � ������ ��������� ������������� (��������� �������):\n\n" << *_matrix;
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
            // ������ ��� ������ ������ - �������������� ������� � ������������ ����
            for (int i = 0; i < _matrixDimention; i++) // �������� �� ���� �������
            {
                double coeff = tempMatrix[i][i]; // ���������� ����������� �� ���������
                for (int j = i; j < _matrixDimention + 1; j++) // �������� �� ���� ��������� ������� ������, ������� ������ �������������
                {
                    std::cout << '\n' << tempMatrix[i][j] << " / " << coeff << " = ";
                    tempMatrix[i][j] /= coeff; // ����������� ��������� ������� � ��������� ��������� (���������������� ������)
                    std::cout << tempMatrix[i][j] << '\n';

                    std::cout << '\n' << std::setw(10) << tempMatrix;

                }
                std::cout << '\n' << std::setw(10) << tempMatrix;
                for (int j = i+1; j < _matrixDimention; j++) 
                {
                    coeff = tempMatrix[j][i]; // ���������� ����������� �������������� 
                    for (int k = i; k < _matrixDimention+1; k++) // �������� �� ���� ��������� ������� ������, ������� ������ �������������
                    {
                        std::cout << '\n' << tempMatrix[j][k] << " - " << coeff << " * " << tempMatrix[i][k] << " = ";
                        tempMatrix[j][k] -= coeff * tempMatrix[i][k]; // ����������� ��������� ������� � ��������� ��������� (���������������� ������)
                        std::cout << tempMatrix[j][k] << '\n';

                        std::cout << '\n' << std::setw(10) << tempMatrix;

                    }
                    std::cout << '\n' << std::setw(10) << tempMatrix;
                }
               
            }
            std::cout << "\n������� � ��������� ����������: \n" << std::setw(10) << tempMatrix;
                
            // �������� ��� ������ ������
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
            // ������������ ������� _matrix �� ������� B � C ���, ��� A = B*C
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