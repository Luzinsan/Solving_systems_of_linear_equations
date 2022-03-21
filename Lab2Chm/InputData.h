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

    double zero_matrix_initer(size_t m, size_t n, size_t r, size_t c);

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
        // �������� ������
        Matrix<double>* _expandedMatrix;
        Matrix<double>* _inverseMatrix;
        Matrix<double>* A;
        Vector<double>* b;
        Vector<double>* x;
        int m; // ����������� ���������� �������
        // �������� � ��������
        Vector<double>* ResidualVector;
        
        int _type;
        double _determinant = 1;
        double _eps = 1e-12;
        int _NAfterComma;
    public:

        InputData()
        {
            _fin  = new std::ifstream("input.txt");
            _fout = new std::ofstream("output.txt", std::ios::app);
            *_fin >> _type;
            std::cout << "\n\t��� ������: " << _type;
            *_fin >> m;
            std::cout << "\n\t������� �������: " << m;

            double* array = new double[m * (m + 1)];
            for (int i = 0; i < (m + 1) * m; i++)
                *_fin >> array[i];

            _expandedMatrix = new Matrix<double>(m, m + 1, array);
            _inverseMatrix = new Matrix<double>(m);
            delete[] array;
            std::cout << "\n\t������������ ��������� ������� � ������ ��������� ������������� (��������� �������):\n\n" << *_expandedMatrix;
            setA(*_expandedMatrix, m);
            setB(*_expandedMatrix, m);
            x = new Vector<double>(m);
            (*x).transposition();

            ResidualVector = new Vector<double>(m);
            (*ResidualVector).transposition();
            
            delete _fin;
        }

        // ������������� ���������� �������
        void setA(const Matrix<double>& matrix, size_t size) 
        {
            A = new Matrix<double>(size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    (*A)[i][j] = matrix[i][j];
        }

        // ������������� ������� ��������� ������������� - initialization of the vector of free coefficients
        void setB(const Matrix<double>& matrix, size_t size) // size - number of rows 
        {
            b = new Vector<double>(size);
            for (int i = 0; i < size; i++)
                (*b)[i] = matrix[i][size];
            (*b).transposition();
        }

        ~InputData()
        {
            delete _expandedMatrix;
            delete _inverseMatrix;
            delete A;
            delete x;
            delete b;
            delete ResidualVector;
            delete _fout;
        }
  

        Matrix<double>& getExpandedMatrix() { return *_expandedMatrix; }
        //������ �������
        void setResidualVector(const Matrix<double>& _A, const Vector<double>& _x, const Vector<double>& _b)
        {
            *_fout << "\n������ ������� e*:"
                << "\nA:\n" << _A
                << "\n* x:\n" << _x
                << "--------------------------\n" << _A * _x
                << "-\n" << _b
                << "--------------------------\n" << (_A * _x) - _b;
            
            (*ResidualVector) = (_A * _x) - _b;
        }


        void setInverseMatrixByMethod(Vector<double>(*Method)(const Matrix<double>&, const Vector<double>&, double& determinant))
        {
            Vector<Vector<double>> x_temp(m);
            Vector<Vector<double>> E(m); 
            for (int i = 0; i < m; i++)
            {
                E[i] = Vector<double>(m);
                for (int j = 0; j < m; j++)
                    if (i == j)
                        E[i][j] = 1;
                    else
                        E[i][j] = 0;
                E[i].transposition();
            }
            //std::cout << "\nE:\n" << E << "\n";

            for (int i = 0; i < m; i++) 
            {
                x_temp[i] = Method(*A, E[i], _determinant);
                //std::cout << "\nx'[" << i << "]=\n" << x_temp[i];
            }
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    (*_inverseMatrix)[i][j] = x_temp[j][i];
            *_fout << "\n�������� �������:\n" << std::setw(10) << (*_inverseMatrix);
        }

        void setGaussMethod() 
        {
            *_fout << "\n����� ������:\n";
            (*x) = GaussMethod(*_expandedMatrix, _determinant);
            *_fout << "\n���������:\n"   << (*x)
                   << "\n������������: " << _determinant << "\n";
            setResidualVector(*A, *x, *b);
            *_fout << "\n��������� ����� ������� �������: " << (*ResidualVector).getModule() << "\n";
        }

        static Vector<double> GaussMethod(const Matrix<double>& _A, const Vector<double>& _b, double& determinant)
        {
            Matrix<double> expandedMatrix(_A.getRows(), _A.getCols() + 1);
            for (int i = 0; i < expandedMatrix.getRows(); i++)
                for (int j = 0; j < expandedMatrix.getCols(); j++)
                    if (j == expandedMatrix.getCols() - 1)
                        expandedMatrix[i][j] = _b[i];
                    else
                        expandedMatrix[i][j] = _A[i][j];
            //*_fout << "\n����������� ������� = \n" << expandedMatrix;
            return GaussMethod(expandedMatrix, determinant);
        }

        static Vector<double> GaussMethod(const Matrix<double>& expandedMatrix, double& determinant)
        {
            // ������� ������������ ����������� �������, ������� ����������� �
            Matrix<double> tempMatrix(expandedMatrix);
            determinant = 1;
            // ������ ��� ������ ������ - �������������� ������� � ������������ ����
            for (int i = 0; i < tempMatrix.getRows(); i++) // �������� �� ���� �������
            {
                std::cout << "i = "<< i << "\n" << std::setw(10) << tempMatrix << "\n";
                double coeff = tempMatrix[i][i]; // ���������� ����������� �� ���������
                determinant *= coeff;
                //std::cout << "\ncoeff=" << coeff << "\n";
                //std::cout << "\n_determinant=" << _determinant << "\n";
                
                for (int j = i; j < tempMatrix.getRows() + 1; j++) // �������� �� ���� ��������� ������� ������, ������� ������ �������������
                {
                    //std::cout << '\n' << tempMatrix[i][j] << " / " << coeff << " = ";
                    // ���� ��� ������� �� ���������, �� �� ����������� � �������, 
                    // � ���� ����� ������ �� ������� ������, �� ������ ������� �� ���� �����������
                    tempMatrix[i][j] /= coeff; 
                    //std::cout << tempMatrix[i][j] << '\n';
                    //std::cout << '\n' << std::setw(10) << tempMatrix;

                }
                std::cout << '\n' << std::setw(10) << tempMatrix;
                for (int j = i + 1; j < tempMatrix.getRows(); j++)
                {
                    coeff = tempMatrix[j][i]; // ���������� ����������� ��������� 
                    for (int k = i; k < tempMatrix.getCols(); k++) // �������� �� ���� ��������� ������, ��������� �������� ������� ����� ����������
                    {
                        //std::cout << '\n' << tempMatrix[j][k] << " - " << coeff << " * " << tempMatrix[i][k] << " = ";
                        tempMatrix[j][k] -= coeff * tempMatrix[i][k]; // �������� �� ������� ������ ������� i-� ������ ����������� �� coeff
                        //std::cout << tempMatrix[j][k] << '\n'; // � ���������� ������� ��������� ���������� ����� ��� ��������
                        //std::cout << '\n' << std::setw(10) << tempMatrix;
                    }
                    //std::cout << '\n' << std::setw(10) << tempMatrix;
                }
            }
            std::cout << "\n������� � ��������� ����������: \n" << std::setw(10) << tempMatrix;
            std::cout << "\n������������ �������: " << determinant << "\n";


            Vector<double> result(tempMatrix.getRows());
            result.transposition();
            // �������� ��� ������ ������
            for (int i = result.getLength() - 1; i >= 0; i--)
            {
                double sumCoeff = 0;
                for (int j = i + 1; j < result.getLength(); j++)
                {
                    //std::cout << "\n" << tempMatrix[i][j] << " * " << result[j] << " + " << sumCoeff << " = ";
                    sumCoeff += tempMatrix[i][j] * result[j];
                   // std::cout << sumCoeff << "\n";
                }
                //std::cout << "result: " << tempMatrix[i][m] << " - " << sumCoeff << " = ";
                result[i] = tempMatrix[i][result.getLength()] - sumCoeff;
                //std::cout << result[i] << "\n";

            }
            return result;
        }
    
        void DecompositionMethod() 
        {
            // ������������ ������� _matrix �� ������� B � C ���, ��� A = B * C
            Matrix<double>* B = new Matrix<double>(m, unit_matrix_initer);
            Matrix<double>* C = new Matrix<double>(m, unit_matrix_initer);
            std::cout << "B: \n" << *B << "\nC:\n" << *C << "\n";
            
            for (int j = 0; j < m; j++)
            {
                // b_ij = a_ij - sum(b_ik*c_kj)
                for (int i = j; i < m; i++)// �������� �� ��������� ������� ������� B
                {
                    double sumCoeff = 0;
                    for (int k = 0; k < j - 1; k++)
                    {
                        std::cout << "\n" << (*B)[i][k] << " * " << (*C)[k][j] << " + " << sumCoeff << " = ";
                        sumCoeff += (*B)[i][k] * (*C)[k][j];
                        std::cout << sumCoeff;
                    }

                    std::cout << "\n" << (*A)[i][j] << " - " << sumCoeff << " = ";
                    (*B)[i][j] = (*A)[i][j] - sumCoeff;
                    std::cout << (*B)[i][j] << "\n";


                    std::cout << "B: \n" << std::setw(10) << *B << "\nC:\n" << std::setw(10) << *C << "\n";
                }

                std::cout << "\n���������� ������� ������� B: " << j << "\n";

                // c_ij = (1/b_ii)*(a_ij - sum(b_ik*c_kj))
                for (int i = j + 1; i < m; i++)
                {
                    double sumCoeff = 0;
                    for (int k = 0; k < i - 1; k++)
                    {
                        std::cout << "\n" << (*B)[j][k] << " * " << (*C)[k][i] << " + " << sumCoeff << " = ";
                        sumCoeff += (*B)[j][k] * (*C)[k][i];
                        std::cout << sumCoeff;

                    }
                    std::cout << "\n(1 / " << (*B)[j][j] << ") * (" << (*A)[j][i] << " - " << sumCoeff << ") = ";
                    (*C)[j][i] = (1 / (*B)[j][j]) * ((*A)[j][i] - sumCoeff);
                    std::cout << (*C)[j][i] << "\n";

                    std::cout << "B: \n" <<std::setw(10) << *B << "\nC:\n" << std::setw(10) << *C << "\n";
                }
                std::cout << "\n���������� ������ ������� C: " << j << "\n";

            }

               
            std::cout << "\nB * C: \n" << std::setw(10) << *B << "\n*\n" << std::setw(10) << *C << "\n" << std::setw(10) << *B*(*C) << "\n";
        
        }
    
        void OrtogonalizationMethod() 
        {
            Matrix<double> A(m+1, zero_matrix_initer);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m + 1; j++)
                {
                    /*if (j == m)
                        A[i][j] = -(*A)[i][j];
                    else
                        A[i][j] = (*_matrix)[i][j];*/
                }
                std::cout << "\nA = \n" << A;
            }
            A[m][m] = 1;
            std::cout << "\nA = \n" << A;
        
        
        }

        void SimpleIterationMethod() 
        {
            Matrix<double> a(m);
            Vector<double> b(m);
            

            Vector<double> x0(b);
            Vector<double> x1(b);
            Vector<double> mod(m);
            do 
            {
                for (int i = 0; i < m; i++)
                {
                    b[i] = (*A)[i][m - 1] / (*A)[i][i];
                    for (int j = 0; j < m; j++)
                        if (i == j)
                            a[i][j] = 0;
                        else
                            a[i][j] = -(*A)[i][j] / (*A)[i][i];
                    std::cout << "\nb=\n" << std::setw(10) << b << "\na=\n" << std::setw(10) << a;
                }
                x1 = b + a * x0;
                 mod = x1 - x0;
            } while (mod.getModule() );

        }
    };
}
#endif