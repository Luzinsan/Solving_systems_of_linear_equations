#ifndef INPUTDATA_H
#define INPUTDATA_H
#include <iostream>
#include <fstream>
#include <cmath>
#include "Matrix.h"
#include "Vector.h"

namespace luMath
{
    double invert_unit_matrix_initer(size_t m, size_t n, size_t r, size_t c);

    double unit_matrix_initer(size_t m, size_t n, size_t r, size_t c);

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
        std::ofstream* _fout;
        Matrix<double>* _matrix;
        Vector<double>* _result;
        int _type;
        int _matrixDimention;
        double _determinant = 1;
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

            _fout = new std::ofstream("output.txt", std::ios::app);
            *_fout << "\n����� ������:\n";
            Matrix<double> tempMatrix(*_matrix);
            // ������ ��� ������ ������ - �������������� ������� � ������������ ����
            for (int i = 0; i < _matrixDimention; i++) // �������� �� ���� �������
            {
                *_fout << "i = "<< i << "\n" << std::setw(10) << tempMatrix << "\n";
                double coeff = tempMatrix[i][i]; // ���������� ����������� �� ���������
                _determinant *= coeff;
                //std::cout << "\ncoeff=" << coeff << "\n";
                //std::cout << "\n_determinant=" << _determinant << "\n";
                
                for (int j = i; j < _matrixDimention + 1; j++) // �������� �� ���� ��������� ������� ������, ������� ������ �������������
                {
                    //std::cout << '\n' << tempMatrix[i][j] << " / " << coeff << " = ";
                    tempMatrix[i][j] /= coeff; // ����������� ��������� ������� � ��������� ��������� (���������������� ������)
                    //std::cout << tempMatrix[i][j] << '\n';

                    //std::cout << '\n' << std::setw(10) << tempMatrix;

                }
                //std::cout << '\n' << std::setw(10) << tempMatrix;
                for (int j = i + 1; j < _matrixDimention; j++)
                {
                    coeff = tempMatrix[j][i]; // ���������� ����������� �������������� 
                    for (int k = i; k < _matrixDimention + 1; k++) // �������� �� ���� ��������� ������� ������, ������� ������ �������������
                    {
                        //std::cout << '\n' << tempMatrix[j][k] << " - " << coeff << " * " << tempMatrix[i][k] << " = ";
                        tempMatrix[j][k] -= coeff * tempMatrix[i][k]; // ����������� ��������� ������� � ��������� ��������� (���������������� ������)
                        //std::cout << tempMatrix[j][k] << '\n';

                        //std::cout << '\n' << std::setw(10) << tempMatrix;

                    }
                    //std::cout << '\n' << std::setw(10) << tempMatrix;
                }
            }
            //std::cout << "\n������� � ��������� ����������: \n" << std::setw(10) << tempMatrix;
            //std::cout << "\n������������ �������: " << _determinant << "\n";

            // �������� ��� ������ ������
            _result = new Vector<double>(_matrixDimention);
            for (int i = _matrixDimention - 1; i >= 0; i--)
            {
                double sumCoeff = 0;
                for (int j = i + 1; j < _matrixDimention; j++) 
                {
                    //std::cout << "\n" << tempMatrix[i][j] << " * " << (*_result)[j] << " + " << sumCoeff << " = ";
                    sumCoeff += tempMatrix[i][j] * (*_result)[j];
                    //std::cout << sumCoeff << "\n";
                }
                //std::cout << "result: " << tempMatrix[i][_matrixDimention] << " - " << sumCoeff << " = ";
                (*_result)[i] = tempMatrix[i][_matrixDimention] - sumCoeff;
                //std::cout << (*_result)[i] << "\n";

            }

            (*_result).transposition();

            std::cout << "\n���������:\n" << *_result;
            std::cout << "\n������������: " << _determinant << "\n";
            *_fout << "\n���������:\n" << *_result
                   << "\n������������: " << _determinant << "\n";
            delete _fout;
            

        }
    
        void DecompositionMethod() 
        {
            // ������������ ������� _matrix �� ������� B � C ���, ��� A = B * C
            Matrix<double>* B = new Matrix<double>(_matrixDimention, unit_matrix_initer);
            Matrix<double>* C = new Matrix<double>(_matrixDimention, unit_matrix_initer);
            std::cout << "B: \n" << *B << "\nC:\n" << *C << "\n";
            
            for(int j = 1; j < _matrixDimention; j++)
                for (int i = j; j < _matrixDimention; j++)
                {
                    double sumCoeff = 0;
                    for (int k = 1; k < j - 1; k++) 
                    {
                        std::cout << "\n" << (*B)[i][k] << " * " << (*C)[k][j] << " + " << sumCoeff << " = ";
                        sumCoeff += (*B)[i][k] * (*C)[k][j];
                        std::cout << sumCoeff;
                    }

                    std::cout << "\n" << (*_matrix)[i][j] << " - " << sumCoeff << " = ";
                    (*B)[i][j] = (*_matrix)[i][j] - sumCoeff;
                    std::cout << (*B)[i][j] << "\n";


                    std::cout << "B: \n" << *B << "\nC:\n" << *C << "\n";
                }

            for(int i = 1; i < _matrixDimention - 1; i++)
                for (int j = i + 1; j < _matrixDimention; j++) 
                {
                    double sumCoeff = 0;
                    for (int k = 1; k < i - 1; i++) 
                    {
                        std::cout << "\n" << (*B)[i][k] << " * " << (*C)[k][j] <<" + "<< sumCoeff << " = ";
                        sumCoeff += (*B)[i][k] * (*C)[k][j];
                        std::cout << sumCoeff;
                        
                    }
                    std::cout << "(1 / " << (*B)[i][i] << ") * (" << (*_matrix)[i][j] << " - " << sumCoeff << ")\n";
                    (*C)[i][j] = (1 / (*B)[i][i]) * ((*_matrix)[i][j] - sumCoeff);
                    std::cout << (*C)[i][j];

                    std::cout << "B: \n" << *B << "\nC:\n" << *C << "\n";
                }
               
            std::cout << "\nB * C: \n" << *B << "\n*\n" << *C << "\n" << *B*(*C) << "\n";
        
        }
    };
}
#endif