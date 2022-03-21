#ifndef INPUTDATA_H
#define INPUTDATA_H
#include <iostream>
#include <fstream>
#include <cmath>
#include "Matrix.h"
#include "Vector.h"

namespace luMath
{
    template <class T>
    T invert_unit_matrix_initer(size_t m, size_t n, size_t r, size_t c);

    template <class T>
    T unit_matrix_initer(size_t m, size_t n, size_t r, size_t c);

    template <class T>
    T zero_matrix_initer(size_t m, size_t n, size_t r, size_t c);

    char getSymbol(std::initializer_list<char> list,
        std::string notification_message = "",
        std::string error_message = "Недопустимое значение, попробуйте ещё раз.\n->");

    double getDouble(double min = -DBL_MAX,
        double max = DBL_MAX,
        std::string notification_message = "",
        std::string error_message = "Недопустимое значение, попробуйте ещё раз.\n->");

    template <class T>
    class InputData
    {
    private:
        std::ifstream* _fin;
        std::ofstream* _fout;
        // основные данные
        Matrix<T>* _expandedMatrix;
        Matrix<T>* _inverseMatrix;
        Matrix<T>* A;
        Vector<T>* b;
        Vector<T>* x;
        int m; // размерность квадратной матрицы
        // проверка и точность
        Vector<T>* ResidualVector;
        
        int _type;
        T _determinant = 1;
        T _eps = 1e-12;
        int _NAfterComma;
    public:

        InputData()
        {
            _fin  = new std::ifstream("input.txt");
            _fout = new std::ofstream("output.txt", std::ios::app);
            *_fin >> _type;
            std::cout << "\n\tТип задачи: " << _type;
            *_fin >> m;
            std::cout << "\n\tПорядок матрицы: " << m;

            T* array = new T[m * (m + 1)];
            for (int i = 0; i < (m + 1) * m; i++)
                *_fin >> array[i];

            _expandedMatrix = new Matrix<T>(m, m + 1, array);
            _inverseMatrix = new Matrix<T>(m);
            delete[] array;
            std::cout << "\n\tКоэффициенты считанной матрицы и вектор свободных коэффициентов (последний стоблик):\n\n" << *_expandedMatrix;
            setA(*_expandedMatrix, m);
            setB(*_expandedMatrix, m);
            x = new Vector<T>(m);
            (*x).transposition();

            ResidualVector = new Vector<T>(m);
            (*ResidualVector).transposition();
            
            delete _fin;
        }

        // инициализация квадратной матрицы
        void setA(const Matrix<T>& matrix, size_t size) 
        {
            A = new Matrix<T>(size);
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    (*A)[i][j] = matrix[i][j];
        }

        // инициализация вектора свободных коэффициентов - initialization of the vector of free coefficients
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
            delete _fout;
        }
  

        Matrix<T>& getExpandedMatrix() { return *_expandedMatrix; }
        //вектор невязки
        void setResidualVector(const Matrix<T>& _A, const Vector<T>& _x, const Vector<T>& _b)
        {
            *_fout << "\nВектор невязки e*:"
                << "\nA:\n" << _A
                << "\n* x:\n" << _x
                << "--------------------------\n" << _A * _x
                << "-\n" << _b
                << "--------------------------\n" << (_A * _x) - _b;
            
            (*ResidualVector) = (_A * _x) - _b;
        }

        const Matrix<T>& getInverseMatrix() const {  return *_inverseMatrix;  }
        const Matrix<T>& getMainMatrix() const { return *A; }

        void setInverseMatrixByMethod(Vector<T>(*Method)(const Matrix<T>&, const Vector<T>&, T& determinant))
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
                x_temp[i] = Method(*A, E[i], _determinant);
                //std::cout << "\nx'[" << i << "]=\n" << x_temp[i];
            }
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    (*_inverseMatrix)[i][j] = x_temp[j][i];
            *_fout << "\nОбратная матрица:\n" << std::setw(10) << (*_inverseMatrix);
        }

        void setGaussMethod() 
        {
            *_fout << "\nМетод Гаусса:\n";
            (*x) = GaussMethod(*_expandedMatrix, _determinant);
            *_fout << "\nРезультат:\n"   << (*x)
                   << "\nОпределитель: " << _determinant << "\n";
            setResidualVector(*A, *x, *b);
            *_fout << "\nЕвклидова норма вектора невязки: " << (*ResidualVector).getModule() << "\n";
        }

        static Vector<T> GaussMethod(const Matrix<T>& _A, const Vector<T>& _b, T& determinant)
        {
            Matrix<T> expandedMatrix(_A.getRows(), _A.getCols() + 1);
            for (int i = 0; i < expandedMatrix.getRows(); i++)
                for (int j = 0; j < expandedMatrix.getCols(); j++)
                    if (j == expandedMatrix.getCols() - 1)
                        expandedMatrix[i][j] = _b[i];
                    else
                        expandedMatrix[i][j] = _A[i][j];
            //*_fout << "\nРасширенная матрица = \n" << expandedMatrix;
            return GaussMethod(expandedMatrix, determinant);
        }

        static Vector<T> GaussMethod(const Matrix<T>& expandedMatrix, T& determinant)
        {
            // удобнее использовать расширенную матрицу, поэтому зафиксируем её
            Matrix<T> tempMatrix(expandedMatrix);
            determinant = 1;
            // Прямой ход метода Гаусса - преобразование матрицы к треугольному виду
            for (int i = 0; i < tempMatrix.getRows(); i++) // проходим по всем строкам
            {
                std::cout << "i = "<< i << "\n" << std::setw(10) << tempMatrix << "\n";
                T coeff = tempMatrix[i][i]; // запоминаем коэффициент по диагонали
                determinant *= coeff;
                //std::cout << "\ncoeff=" << coeff << "\n";
                //std::cout << "\n_determinant=" << _determinant << "\n";
                
                for (int j = i; j < tempMatrix.getRows() + 1; j++) // проходим по всем элементам текущей строки, включая вектор коэффициентов
                {
                    //std::cout << '\n' << tempMatrix[i][j] << " / " << coeff << " = ";
                    // если это элемент на диагонали, то он вырождается в единицу, 
                    // а если любой другой на текущей строке, то просто делится на этот коэффициент
                    tempMatrix[i][j] /= coeff; 
                    //std::cout << tempMatrix[i][j] << '\n';
                    //std::cout << '\n' << std::setw(10) << tempMatrix;

                }
                std::cout << '\n' << std::setw(10) << tempMatrix;
                for (int j = i + 1; j < tempMatrix.getRows(); j++)
                {
                    coeff = tempMatrix[j][i]; // запоминаем коэффициент умножения 
                    for (int k = i; k < tempMatrix.getCols(); k++) // проходим по всем элементам строки, некоторые элементы которой будут обнуляться
                    {
                        //std::cout << '\n' << tempMatrix[j][k] << " - " << coeff << " * " << tempMatrix[i][k] << " = ";
                        tempMatrix[j][k] -= coeff * tempMatrix[i][k]; // вычитаем из текущей строки верхнюю i-ю строку помноженную на coeff
                        //std::cout << tempMatrix[j][k] << '\n'; // в результате получим некоторое количество нулей под единицей
                        //std::cout << '\n' << std::setw(10) << tempMatrix;
                    }
                    //std::cout << '\n' << std::setw(10) << tempMatrix;
                }
            }
            std::cout << "\nМатрица с единичной диагональю: \n" << std::setw(10) << tempMatrix;
            std::cout << "\nОпределитель матрицы: " << determinant << "\n";


            Vector<T> result(tempMatrix.getRows());
            result.transposition();
            // Обратный ход метода Гаусса
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
    
        
        void DecompositionMethod() 
        {
            // Раскладываем матрицу _matrix на матрицы B и C так, что A = B * C
            Matrix<T>* B = new Matrix<T>(m, unit_matrix_initer<T>);
            Matrix<T>* C = new Matrix<T>(m, unit_matrix_initer<T>);
            std::cout << "B: \n" << *B << "\nC:\n" << *C << "\n";
            
            for (int j = 0; j < m; j++)
            {
                // b_ij = a_ij - sum(b_ik*c_kj)
                for (int i = j; i < m; i++)// проходим по элементам столбца матрицы B
                {
                    T sumCoeff = 0;
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

                std::cout << "\nОбработали столбец матрицы B: " << j << "\n";

                // c_ij = (1/b_ii)*(a_ij - sum(b_ik*c_kj))
                for (int i = j + 1; i < m; i++)
                {
                    T sumCoeff = 0;
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
                std::cout << "\nОбработали строку матрицы C: " << j << "\n";

            }

               
            std::cout << "\nB * C: \n" << std::setw(10) << *B << "\n*\n" << std::setw(10) << *C << "\n" << std::setw(10) << *B*(*C) << "\n";
        
        }
    
        
        void OrtogonalizationMethod() 
        {
            Matrix<T> A(m+1, zero_matrix_initer<T>);
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