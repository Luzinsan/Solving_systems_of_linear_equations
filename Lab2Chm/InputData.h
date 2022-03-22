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

    template <class T>
    class InputData
    {
    public:
        enum class METHOD
        {
            GAUSS=1,
            DECOMPOSOTION
        };

        enum class TASK
        {
            ROOT=1,
            DETERMINANT,
            INVERS
        };

    private:
        std::ifstream* _fin;
        std::ofstream* _fout;
        // �������� ������
        Matrix<T>* _expandedMatrix;
        Matrix<T>* _inverseMatrix;
        Matrix<T>* A;
        Vector<T>* b;
        Vector<T>* x;
        int m; // ����������� ���������� �������
        // �������� � ��������
        Vector<T>* ResidualVector;
        Matrix<T>* ResidualMatrix;

        METHOD _method;
        TASK _task;
        
        T _determinant = 0;
        T _eps = 1e-12;
        int _NAfterComma;
    public:
        

        InputData()
        {
            _fin  = new std::ifstream("input.txt");
            _fout = new std::ofstream("output.txt"/*, std::ios::app*/);
            int c;
            *_fin >> c; // ����������� �����
            _method = static_cast<METHOD>(c);
            *_fin >> c; // ����������� ��� ������
            _task = static_cast<TASK>(c);
           
            *_fin >> m;
            std::cout << "\n\t������� �������: " << m;

            T* array = new T[m * (m + 1)];
            for (int i = 0; i < (m + 1) * m; i++)
                *_fin >> array[i];

            _expandedMatrix = new Matrix<T>(m, m + 1, array);
            _inverseMatrix = new Matrix<T>(m);
            delete[] array;
            std::cout << "\n\t������������ ��������� ������� � ������ ��������� ������������� (��������� �������):\n\n" << *_expandedMatrix;
            setA(*_expandedMatrix, m);
            setB(*_expandedMatrix, m);
            x = new Vector<T>(m);
            (*x).transposition();

            ResidualVector = new Vector<T>(m);
            (*ResidualVector).transposition();
            ResidualMatrix = new Matrix<T>(m);
            
            delete _fin;
        }

        METHOD getMethod() { return _method; }
        TASK getTask() { return _task;  }

        // ������������� ���������� �������
        void setA(const Matrix<T>& matrix, size_t size) 
        {
            A = new Matrix<T>(size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    (*A)[i][j] = matrix[i][j];
        }

        // ������������� ������� ��������� ������������� - initialization of the vector of free coefficients
        void setB(const Matrix<T>& matrix, size_t size) // size - number of rows 
        {
            b = new Vector<T>(size);
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
            delete ResidualMatrix;
            delete _fout;
        }
  

        Matrix<T>& getExpandedMatrix() { return *_expandedMatrix; }
        
        //������ �������
        void setResidualVector(const Matrix<T>& _A, const Vector<T>& _x, const Vector<T>& _b)
        {
            *_fout << "\n������ ������� e*:\n"
               /* << "A:\n" << _A
                << "\n* x:\n" << _x
                << "--------------------------\n" << _A * _x
                << "-\n" << _b
                << "--------------------------\n"*/ << (_A * _x) - _b;
            
            (*ResidualVector) = (_A * _x) - _b;
        }

        const Matrix<T>& getInverseMatrix() const {  return *_inverseMatrix;  }
        const Matrix<T>& getMainMatrix() const { return *A; }

        void setInverseMatrixByMethod(Vector<T>(*Method)(const Matrix<T>&, const Vector<T>&, T& determinant, std::ofstream&))
        {
            Vector<Vector<T>> x_temp(m);
            Vector<Vector<T>> E(m); 
            for (int i = 0; i < m; i++)
            {
                E[i] = Vector<T>(m);
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
                x_temp[i] = Method(*A, E[i], _determinant,*_fout);
                *_fout  << "\nx'[" << i << "]=\n" << x_temp[i];
            }
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    (*_inverseMatrix)[i][j] = x_temp[j][i];
            *_fout << "\n�������� �������:\n" << std::setw(10) << (*_inverseMatrix);
            
            *ResidualMatrix = (*_inverseMatrix) * (*A) - Matrix<double>(m, unit_matrix_initer);
            *_fout << "\n������� �������:\n" /*<< std::fixed*/ << std::setprecision(5) << std::setw(15) << *ResidualMatrix;
            *_fout << "\n��������� ����� ������� �������: " << (*ResidualMatrix).getModule() << "\n";
               
        }

        void setRoot(Vector<T>(*Method)(const Matrix<T>&, const Vector<T>&, T& determinant, std::ofstream&))
        {
            (*x) = Method(*A, *b, _determinant, *_fout);
        }
        void getRoot(Vector<T>(*Method)(const Matrix<T>&, const Vector<T>&, T& determinant, std::ofstream&))
        {
            if(Method == InputData::GaussMethod)
                *_fout << "\n����� ������:\n";
            else if(Method == InputData::DecompositionMethod)
                *_fout << "\n����� ������������:\n";
               
            setRoot(Method);
            *_fout << "\n���������:\n" << (*x);
            setResidualVector(*A, *x, *b);
            *_fout << "\n��������� ����� ������� �������: " << (*ResidualVector).getModule() << "\n";
        }

        T getDeterminant() 
        {
            if (_determinant == 0) 
            {
                switch (_method)
                {
                case METHOD::GAUSS:
                    setRoot(InputData::GaussMethod);
                    break;
                case METHOD::DECOMPOSOTION:
                    setRoot(InputData::DecompositionMethod);
                    break;
                }
            }
            *_fout << "\n������������: " << _determinant << "\n";
            return _determinant;
        }

        
        static Vector<T> GaussMethod(const Matrix<T>& A, const Vector<T>& b, T& determinant, std::ofstream& out = std::cout)
        {
            Matrix<T> expandedMatrix(A.getRows(), A.getCols() + 1);
            for (int i = 0; i < expandedMatrix.getRows(); i++)
                for (int j = 0; j < expandedMatrix.getCols(); j++)
                    if (j == expandedMatrix.getCols() - 1)
                        expandedMatrix[i][j] = b[i];
                    else
                        expandedMatrix[i][j] = A[i][j];
            //*_fout << "\n����������� ������� = \n" << expandedMatrix;
            return GaussMethod(expandedMatrix, determinant,out);
        }

        static Vector<T> GaussMethod(const Matrix<T>& expandedMatrix, T& determinant, std::ofstream& out=std::cout)
        {
            // ������� ������������ ����������� �������, ������� ����������� �
            Matrix<T> tempMatrix(expandedMatrix);
            determinant = 1;
            // ������ ��� ������ ������ - �������������� ������� � ������������ ����
            for (int i = 0; i < tempMatrix.getRows(); i++) // �������� �� ���� �������
            {
                out << "i = "<< i << "\n" << std::setw(10) << tempMatrix << "\n";
                T coeff = tempMatrix[i][i]; // ���������� ����������� �� ���������
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
                out << '\n' << std::setw(10) << tempMatrix;
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
            //std::cout << "\n������� � ��������� ����������: \n" << std::setw(10) << tempMatrix;
            //std::cout << "\n������������ �������: " << determinant << "\n";


            Vector<T> result(tempMatrix.getRows());
            result.transposition();
            // �������� ��� ������ ������
            for (int i = result.getLength() - 1; i >= 0; i--)
            {
                T sumCoeff = 0;
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
    
        
        static Vector<T> DecompositionMethod(const Matrix<T>& A, const Vector<T>& b, T& determinant, std::ofstream& out = std::cout)
        {
            // ������������ ������� _matrix �� ������� B � C ���, ��� A = B * C
            Matrix<T> B(A.getRows(), unit_matrix_initer);
            Matrix<T> C(A.getRows(), unit_matrix_initer);
            //std::cout << "B: \n" << B << "\nC:\n" << C << "\n";
            unsigned m = B.getRows();
            for (int j = 0; j < m; j++)
            {
                // b_ij = a_ij - sum(b_ik*c_kj)
                for (int i = j; i < m; i++)// �������� �� ��������� ������� ������� B
                {
                    T sumCoeff = 0;
                    for (int k = 0; k <= j - 1; k++)
                    {
                        //std::cout << "\n" << B[i][k] << " * " << C[k][j] << " + " << sumCoeff << " = ";
                        sumCoeff += B[i][k] * C[k][j];
                        //std::cout << sumCoeff;
                    }
                    //std::cout << "\n" << (*A)[i][j] << " - " << sumCoeff << " = ";
                    B[i][j] = A[i][j] - sumCoeff;
                    //std::cout << B[i][j] << "\n";
                    //std::cout << "B: \n" << std::setw(10) << B << "\nC:\n" << std::setw(10) << C << "\n";
                }
                //std::cout << "\n���������� ������� ������� B: " << j << "\n";
                // c_ij = (1/b_ii)*(a_ij - sum(b_ik*c_kj))
                for (int i = j + 1; i < m; i++)
                {
                    T sumCoeff = 0;
                    for (int k = 0; k <= j - 1; k++)
                    {
                        //std::cout << "\n" << B[j][k] << " * " << C[k][i] << " + " << sumCoeff << " = ";
                        sumCoeff += B[j][k] * C[k][i];
                        //std::cout << sumCoeff;

                    }
                    //std::cout << "\n(1 / " << B[j][j] << ") * (" << (*A)[j][i] << " - " << sumCoeff << ") = ";
                    C[j][i] = (1 / B[j][j]) * (A[j][i] - sumCoeff);
                    //std::cout << C[j][i] << "\n";
                    //std::cout << "B: \n" <<std::setw(10) << B << "\nC:\n" << std::setw(10) << C << "\n";
                }
                //std::cout << "\n���������� ������ ������� C: " << j << "\n";
                out << "\n\ti = "<< j<<"\n\t������� B : \n" << std::setprecision(5) << std::setw(15) << B << "\n\t������� C : \n"  << std::setw(15) << C;
            }

            //std::cout << "\nB * C: \n" << std::setw(10) << B << "\n*\n" << std::setw(10) << C << "\n" << std::setw(10) << B*C << "\n";
            Vector<T> y(m);
            determinant = 1;
            for (int i = 0; i < m; i++)
            {
                T sumCoeff = 0;
                for (int k = 0; k <= i - 1; k++)
                    sumCoeff += B[i][k] * y[k];
                //std::cout << "\ni=" << i << ": y = " << " ( " << (*b)[i] << " - " << sumCoeff << ") / " << B[i][i] << " = ";
                y[i] = (b[i] - sumCoeff) / B[i][i];
                determinant *= B[i][i];
                //std::cout << y[i];
            }
            out << "\ny' = \n" <<y.transposition() << "\n";
            Vector<T> x(m);
            for (int i = m - 1; i >= 0; i--)
            {
                T sumCoeff = 0;
                for (int k = i+1; k < m; k++)
                    sumCoeff += C[i][k] * x[k];
                //std::cout << "\ni=" << i << ": x = "  << y[i] << " - " << sumCoeff << " = ";
                x[i] = y[i] - sumCoeff;
                //std::cout << x[i];
            }
            //std::cout << "\nx=" << x << "\n";
            x.transposition();
            return x;

        }
    
        
        void OrtogonalizationMethod() 
        {
            Matrix<T> A(m+1, zero_matrix_initer);
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