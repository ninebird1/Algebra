#pragma once

#include <math.h>
#include <memory>

namespace rs
{
    inline bool nearZero(float a)
    {
        if (fabs(a) <= 1e-6)
        {
            return true;
        }
        return false;
    }

    class Complex
    {
    public:
        float   m_real;
        float   m_imaginary;
        Complex()
        {
            m_real = 0;
            m_imaginary = 0;
        }
        Complex(float real)
        {
            m_real = real;
            m_imaginary = 0;
        }
        Complex(float real, float imagin)
        {
            m_real = real;
            m_imaginary = imagin;
        }
        Complex(const Complex& c)
        {
            m_real = c.m_real;
            m_imaginary = c.m_imaginary;
        }
        Complex& operator = (const Complex& c)
        {
            m_real = c.m_real;
            m_imaginary = c.m_imaginary;
            return *this;
        }
        void add(const Complex& b)
        {
            m_real += b.m_real;
            m_imaginary += b.m_imaginary;
        }
        void sub(const Complex& b)
        {
            m_real -= b.m_real;
            m_imaginary -= b.m_imaginary;
        }
        void mul(const Complex& b)
        {
            float realNum = m_real * b.m_real - m_imaginary * b.m_imaginary;
            m_imaginary = m_real * b.m_imaginary + m_imaginary * b.m_real;
            m_real = realNum;
        }
        void mul(float b)
        {
            m_imaginary *= b;
            m_real *= b;
        }
        void div(const Complex& b)
        {
            Complex c = b.conjugate();
            this->mul(c);
            double n = b.norm();
            n = n * n;
            if (!nearZero(n))
            {
                m_real /= n;
                m_imaginary /= n;
            }
        }

        Complex add(const Complex& b) const
        {
            Complex a = *this;
            a.add(b);
            return a;
        }
        Complex sub(const Complex& b) const
        {
            Complex a = *this;
            a.sub(b);
            return a;
        }
        Complex mul(const Complex& b) const
        {
            Complex a = *this;
            a.mul(b);
            return a;
        }
        Complex div(const Complex& b) const
        {
            Complex a = *this;
            a.div(b);
            return a;
        }

        void conjugate()
        {
            m_imaginary = -m_imaginary;
        }
        float norm() const
        {
            return sqrt(m_real * m_real + m_imaginary * m_imaginary);
        }
        Complex conjugate() const
        {
            Complex a = *this;
            a.conjugate();
            return a;
        }
        Complex operator -() const
        {
            Complex a(-m_real, -m_imaginary);
            return a;
        }
        operator float()
        {
            return norm();
        }
    };

    Complex operator + (const Complex& a, const Complex& b)
    {
        return a.add(b);
    }
    Complex& operator += (Complex& a, const Complex& b)
    {
        a.add(b);
        return a;
    }
    Complex operator - (const Complex& a, const Complex& b)
    {
        return a.sub(b);
    }
    Complex& operator -= (Complex& a, const Complex& b)
    {
        a.sub(b);
        return a;
    }
    Complex operator * (const Complex& a, const Complex& b)
    {
        return a.mul(b);
    }
    Complex& operator *= (Complex& a, const Complex& b)
    {
        a.mul(b);
        return a;
    }
    Complex operator / (const Complex& a, const Complex& b)
    {
        return a.div(b);
    }
    Complex& operator /= (Complex& a, const Complex& b)
    {
        a.div(b);
        return a;
    }
   
    template <typename Number>
    Complex operator * (const Complex& a, Number b)
    {
        Complex c = a;
        c.m_real *= b;
        c.m_imaginary *= b;
        return c;
    }
    template <typename Number>
    Complex operator * (Number a, const Complex& b)
    {
        Complex c = b;
        c.m_real *= a;
        c.m_imaginary *= a;
        return c;
    }
    template <typename Number>
    Complex operator / (const Complex& a, Number b)
    {
        Complex c = a;
        c.m_real /= b;
        c.m_imaginary /= b;
        return c;
    }
    template <typename Number>
    Complex operator / (Number a, const Complex& b)
    {
        Complex c = b.conjugate();
        float n = b.norm();
        n = n * n;
        float d = a / n;
        c.mul(d);
        return c;
    }

    template <class T = float>
    class Matrix;

    template <class T>
    class Vector
    {
    public:
        Vector(T* t, int N)
        {
            m_N = N;
            m_data = new T[N];
            memcpy(m_data, t, sizeof(T) * N);
        }
        explicit Vector(int N)
        {
            m_N = N;
            m_data = new T[N];
            memset(m_data, 0, sizeof(T) * N);
        }
        Vector(T t, int N)
        {
            m_N = N;
            m_data = new T[N];
            for (int i = 0; i < N; ++i)
            {
                m_data[i] = t;
            }
        }
        Vector& operator = (const Vector& b)
        {
            for (int i = 0; i < m_N; ++i)
            {
                m_data[i] = b.m_data[i];
            }
            return *this;
        }
        ~Vector()
        {
            delete m_data;
            m_data = nullptr;
        }
    public:
        int len() const
        {
            return m_N;
        }
        void normalize()
        {
            float norm = magnitude();
            if (norm > 0)
            {
                for (int i = 0; i < m_N; ++i)
                {
                    m_data[i] /= norm;
                }
            }
        }
        T& operator [](int i)
        {
            return m_data[i];
        }
        void mul(Matrix<T>& m)
        {
            if (m_N != m.m_row)
            {
                return;
            }
            for (int i = 0; i < m.m_column; ++i)
            {
                auto v = m.column(i);
                m_data[i] = this->dot(v);
            }
        }
        void mul(T t)
        {
            for (int i = 0; i < m_N; ++i)
            {
                m_data[i] *= t;
            }
        }
    public:
        std::shared_ptr<Matrix<T>> mulTranspose() const
        {
            Matrix<T> a(m_data, m_N, 1);
            Matrix<T> b(m_data, 1, m_N);
            return a.mul(b);
        }
        float magnitude() const
        {
            float norm = 0;
            for (int i = 0; i < m_N; ++i)
            {
                norm += m_data[i] * m_data[i];
            }
            norm = sqrt(norm);
            return norm;
        }
        T operator [](int i) const
        {
            return m_data[i];
        }
        T* minItem(bool onlyMag = true) const
        {
            T m = m_data[0];
            int k = 0;
            for (int i = 1; i < m_N; ++i)
            {
                if (m > m_data[i])
                {
                    m = m_data[i];
                    k = i;
                }
            }
            return m_data + k;
        }
        float magnitudeMinItem(int* pos) const
        {
            float m = abs(m_data[0]);
            int k = 0;
            for (int i = 1; i < m_N; ++i)
            {
                if (m > abs(m_data[i]))
                {
                    m = abs(m_data[i]);
                    k = i;
                }
            }
            if (pos)
            {
                *pos = k;
            }
            return m;
        }
        T* maxItem() const
        {
            T m = m_data[0];
            int k = 0;
            for (int i = 1; i < m_N; ++i)
            {
                if (m < m_data[i])
                {
                    m = m_data[i];
                    k = i;
                }
            }
            return m_data + k;
        }
        float magnitudeMaxItem(int* pos) const
        {
            float m = abs(m_data[0]);
            int k = 0;
            for (int i = 1; i < m_N; ++i)
            {
                if (m < abs(m_data[i]))
                {
                    m = abs(m_data[i]);
                    k = i;
                }
            }
            if (pos)
            {
                *pos = k;
            }
            return m;
        }
        T dot(const Vector& b) const
        {
            T d = 0;
            for (int i = 0; i < m_N; ++i)
            {
                d += m_data[i] * b.m_data[i];
            }
            return d;
        }
        bool equal(const Vector& b) const
        {
            if (m_N != b.m_N)
            {
                return false;
            }
            for (int i = 0; i < m_N; ++i)
            {
                if (!nearZero(m_data[i] - b.m_data[i]))
                {
                    return false;
                }
            }
            return true;
        }
        bool isNature() const
        {
            float norm = magnitude();
            if (!nearZero(norm - 1))
            {
                return false;
            }
            int k = 0;
            magnitudeMaxItem(&k);
            if (!nearZero(abs(m_data[k]) - 1))
            {
                return false;
            }
            for (int i = 0; i < m_N; ++i)
            {
                if (i != k)
                {
                    if (!nearZero(m_data[i]))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
    private:
        T*      m_data;
        int     m_N;
        friend class Matrix<T>;
    };

    template <class T>
    class Matrix
    {
        typedef typename std::conditional<std::is_integral<T>::value, float, T>::type E;
    public:
        Matrix(int m, int n)
        {
            m_row = m;
            m_column = n;
            m_data = new T*[m_row];
            for (int i = 0; i < m_row; i++)
            {
                m_data[i] = new T[m_column];
                memset(m_data[i], 0, sizeof(T) * m_column);
            }
        }
        Matrix(T* t, int m, int n)
        {
            m_row = m;
            m_column = n;
            m_data = new T*[m_row];
            for (int i = 0; i < m_row; i++)
            {
                m_data[i] = new T[m_column];
                for (int j = 0; j < m_column; j++)
                {
                    m_data[i][j] = t[i * m_column + j];
                }
            }
        }
        Matrix(const Matrix& b)
        {
            m_row = b.m_row;
            m_column = b.m_column;
            m_data = new T*[m_row];
            for (int i = 0; i < m_row; i++)
            {
                m_data[i] = new T[m_column];
                for (int j = 0; j < m_column; j++)
                {
                    m_data[i][j] = b[i][j];
                }
            }
        }
        Matrix& operator = (const Matrix& b)
        {
            for (int i = 0; i < m_row; i++)
            {
                for (int j = 0; j < m_column; j++)
                {
                    m_data[i][j] = b.m_data[i][j];
                }
            }
            return *this;
        }

        template <typename SrcType>
        Matrix& operator = (const Matrix<SrcType>& b)
        {
            for (int i = 0; i < m_row; i++)
            {
                for (int j = 0; j < m_column; j++)
                {
                    m_data[i][j] = b[i][j];
                }
            }
            return *this;
        }

        ~Matrix()
        {
            for (int i = 0; i < m_row; i++)
            {
                delete m_data[i];
            }
            delete m_data;
            m_data = nullptr;
            m_row = 0;
            m_column = 0;
        }
        Matrix& identity()
        {
            int n = fmin(m_row, m_column);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        m_data[i][j] = 1;
                    }
                    else
                    {
                        m_data[i][j] = 0;
                    }
                }
            }
            return *this;
        }
        Matrix& mul(E a)
        {
            for (int i = 0; i < m_row; i++)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    m_data[i][j] *= a;
                }
            }
            return *this;
        }
        Matrix& mulRow(int i, E a)
        {
            for (int j = 0; j < m_column; ++j)
            {
                m_data[i][j] *= a;
            }
            return *this;
        }
        Matrix& add(const Matrix& b)
        {
            for (int i = 0; i < m_row; i++)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    m_data[i][j] += b.m_data[i][j];
                }
            }
            return *this;
        }
        Matrix& add(E a)
        {
            for (int i = 0; i < m_row; ++i)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    m_data[i][j] += a;
                }
            }
            return *this;
        }
        Matrix& replace(const Matrix& mat, int m0 = 0, int n0 = 0)
        {
            for (int i = 0; i < mat.m_row; i++)
            {
                for (int j = 0; j < mat.m_column; j++)
                {
                    m_data[i + m0][j + n0] = mat.m_data[i][j];
                }
            }
            return *this;
        }
        template <typename U>
        Matrix& replaceRow(const Vector<U>& v, int i)
        {
            if (v.m_N != m_column)
            {
                return *this;
            }
            for (int j = 0; j < m_column; ++j)
            {
                m_data[i][j] = v[j];
            }
            return *this;
        }
        template <typename U>
        Matrix& replaceColumn(const Vector<U>& v, int i)
        {
            if (v.len() != m_row)
            {
                return *this;
            }
            for (int j = 0; j < m_row; ++j)
            {
                m_data[j][i] = v[j];
            }
            return *this;
        }
        Matrix& addRow(int dst, int src, E k, int n0 = 0)
        {
            for (int i = n0; i < m_column; ++i)
            {
                m_data[dst][i] += m_data[src][i] * k;
            }
            return *this;
        }
        Matrix& swapRow(int a, int b)
        {
            if (a == b)
            {
                return *this;
            }
            for (int i = 0; i < m_column; ++i)
            {
                std::swap(m_data[a][i], m_data[b][i]);
            }
            return *this;
        }
        Matrix& simplest()
        {
            for (int i = 0; i < m_row; ++i)
            {
                //find first none zero item in row i column i
                T tmp = m_data[i][i];
                if (nearZero(tmp))
                {
                    for (int j = i + 1; j < m_row; ++j)
                    {
                        tmp = m_data[j][i];
                        if (!nearZero(tmp))
                        {
                            //swap with the row that ij is not zero
                            swapRow(i, j);
                            break;
                        }
                    }
                }
                if (!nearZero(tmp))
                {
                    this->mulRow(i, 1 / tmp);
                    for (int j = 0; j < m_row; j++)
                    {
                        if (j != i)
                        {
                            this->addRow(j, i, -m_data[j][i]);
                        }
                    }
                }
            }
            return *this;
        }
        Matrix& inverse()
        {
            T a;
            inverse(a);
            return *this;
        }
        bool isOrthogonal() const
        {
            for (int i = 0; i < m_column; ++i)
            {
                auto v = this->column(i);
                if (!nearZero(v->magnitude() - 1))
                {
                    return false;
                }
            }
            for (int i = 0; i < m_column - 1; ++i)
            {
                auto u = this->column(i);
                for (int j = i + 1; j < m_column; ++j)
                {
                    auto v = this->column(j);
                    T d = u->dot(*v);
                    if (!nearZero(d))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        operator T** () const
        {
            return m_data;
        }
        T det() const
        {
            if (m_row != m_column)
            {
                return 0;
            }
            std::shared_ptr<Matrix> L, U;
            int swapCnt = 0;
            LUDecomposite(L, U, &swapCnt);
            T d = 1;
            for (int i = 0; i < m_row; ++i)
            {
                d *= (*U)[i][i];
            }
            if (nearZero(d))
            {
                return d;
            }
            //行交换
            if (swapCnt % 2)
            {
                d = -d;
            }
            return d;
        }
        bool isSysmetric() const
        {
            if (m_row != m_column)
            {
                return false;
            }
            T a;
            for (int i = 0; i < m_row; ++i)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    if (!isSysmetric(i, j, a))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
    public:
        std::shared_ptr<Vector<T>> row(int i, int from = 0, int to = 0) const
        {
            if (to == 0)
            {
                to = m_column;
            }
            int columns = to - from;
            if (columns > m_column)
            {
                return nullptr;
            }
            auto v = std::make_shared<Vector<T>>(m_data[i] + from, columns);
            return v;
        }
        std::shared_ptr<Vector<T>> column(int j, int from = 0, int to = 0) const
        {
            if (to == 0)
            {
                to = m_row;
            }
            int row = to - from;
            if (row > m_row)
            {
                return nullptr;
            }
            auto v = std::make_shared<Vector<T>>(row);
            for (int i = 0; i < row; i++)
            {
                v->m_data[i] = m_data[i + from][j];
            }
            return v;
        }
        std::shared_ptr<Vector<T>> mul(const Vector<T>& v) const
        {
            if (m_column != v.m_N)
            {
                return nullptr;
            }
            auto P = std::make_shared<Vector<T>>(v.m_N);
            for (int i = 0; i < m_row; ++i)
            {
                auto r = row(i);
                P->m_data[i] = r->dot(v);
            }
            return P;
        }
        std::shared_ptr<Vector<T>> diagonal() const
        {
            int m = fmin(m_row, m_column);
            auto v = std::make_shared<Vector<T>>(m);
            for (int i = 0; i < m; ++i)
            {
                v->m_data[i] = m_data[i][i];
            }
            return v;
        }
    public:
        std::shared_ptr<Matrix> mul(const Matrix& b) const
        {
            if (m_column != b.m_row)
            {
                return nullptr;
            }
            auto c = std::make_shared<Matrix<T>>(m_row, b.m_column);
            for (int i = 0; i < m_row; i++)
            {
                for (int j = 0; j < b.m_column; ++j)
                {
                    for (int k = 0; k < m_column; ++k)
                    {
                        c->m_data[i][j] += m_data[i][k] * b.m_data[k][j];
                    }
                }
            }
            return c;
        }
        std::shared_ptr<Matrix> subMatrix(int m, int n, int m0 = 0, int n0 = 0) const
        {
            auto mat = std::make_shared<Matrix<T>>(m, n);
            mat->identity();
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; ++j)
                {
                    mat->m_data[i][j] = m_data[i + m0][j + n0];
                }
            }
            return mat;
        }
        std::shared_ptr<Matrix> expand(int m, int n, int m0 = 0, int n0 = 0) const
        {
            auto mat = std::make_shared<Matrix<T>>(m, n);
            mat->identity();
            for (int i = 0; i < m_row; i++)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    mat->m_data[i + m0][j + n0] = m_data[i][j];
                }
            }
            return mat;
        }
        std::shared_ptr<Matrix> transpose() const
        {
            auto c = std::make_shared<Matrix<T>>(m_column, m_row);
            for (int i = 0; i < m_column; i++)
            {
                for (int j = 0; j < m_row; ++j)
                {
                    c->m_data[i][j] = transpose<T>(i, j);
                }
            }
            return c;
        }
        std::shared_ptr<Matrix> combine(const Matrix& b) const
        {
            if (m_row != b.m_row)
            {
                return nullptr;
            }
            std::shared_ptr<Matrix> P = std::make_shared<Matrix>(m_row, m_column + b.m_column);
            for (int i = 0; i < m_row; ++i)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    P->m_data[i][j] = m_data[i][j];
                }
                for (int j = 0; j < b.m_column; ++j)
                {
                    P->m_data[i][m_column + j] = b.m_data[i][j];
                }
            }
            return P;
        }
        std::shared_ptr<Matrix> simplest() const
        {
            auto M = std::make_shared<Matrix>(m_row, m_column);
            *M = *this;
            M->simplest();
            return M;
        }
        std::shared_ptr<Matrix> inverse() const
        {
            auto M = std::make_shared<Matrix>(m_row, m_column);
            *M = *this;
            M->inverse();
            return M;
        }
        std::shared_ptr<Matrix> add(const Matrix& b) const
        {
            auto c = std::make_shared<Matrix>(m_row, m_column);
            *c = *this;
            c->add(b);
            return c;
        }
    public:
        void QRDecomposite(std::shared_ptr<Matrix>& Q, std::shared_ptr<Matrix>& R) const
        {
            householder(Q, R);
        }
        void LUDecomposite(std::shared_ptr<Matrix>& L, std::shared_ptr<Matrix>& U, int* rowSwapCnt = nullptr) const
        {
            L = std::make_shared<Matrix>(m_row, m_column);
            U = std::make_shared<Matrix>(m_row, m_column);
            auto P = std::make_shared<Matrix>(m_row, m_column);
            L->identity();
            P->identity();
            *U = *this;
            int swapCnt = 0;
            //预处理，获取P矩阵
            for (int i = 0; i < m_row - 1; ++i)
            {
                //主非0元
                if (nearZero(U->m_data[i][i]))
                {
                    for (int j = i + 1; j < m_row; ++j)
                    {
                        //swap two row
                        if (!nearZero(U->m_data[j][i]))
                        {
                            P->swapRow(i, j);
                            ++swapCnt;
                            break;
                        }
                    }
                }
                for (int j = i + 1; j < m_row; ++j)
                {
                    E k = -(U->m_data[j][i] / U->m_data[i][i]);
                    U->addRow(j, i, k, i);
                }
            }
            if (rowSwapCnt)
            {
                *rowSwapCnt = swapCnt;
            }
            //permutate
            U = P->mul(*this);
            for (int i = 0; i < m_row - 1; ++i)
            {
                //a[i][i] must not equal to zero after swapped row by left multiply P matrix
                for (int j = i + 1; j < m_row; ++j)
                {
                    E k = -(U->m_data[j][i] / U->m_data[i][i]);
                    U->addRow(j, i, k, i);
                    L->m_data[j][i] = -k;
                }
            }
            //PA=LU -> A=P(-1)LU
            L = P->transpose()->mul(*L);
        }
        void SingularValueDecomposite(std::shared_ptr<Matrix<Complex>>& U, std::shared_ptr<Matrix<float>>& E, std::shared_ptr<Matrix<Complex>>& V) const
        {
            const Matrix A = *this;
            const Matrix B = *(A.transpose());
            auto BA = B.mul(A);

            //eigen vectors of left singular matrix
            std::vector<Complex> BAEigenValues, ABEigenValues;
            V = BA->eigen(BAEigenValues);

            //eigen vectors of right singular matrix
            auto AB = A.mul(B);
            U = AB->eigen(ABEigenValues);

            //matrix E
            E = std::make_shared<Matrix<float>>(m_row, m_column);
            int n = fmin(m_row, m_column);
            for (int i = 0; i < n; ++i)
            {
                (*E)[i][i] = sqrt(ABEigenValues[i]);
            }
        }

        template <typename EigenValueType>
        std::shared_ptr<Vector<EigenValueType>> inversePowerMethod(EigenValueType eigenValue) const
        {
            Matrix<EigenValueType> A(m_row, m_column);
            A = *this;
            Matrix<EigenValueType> I(m_row, m_column);
            I.identity();
            I.mul(-eigenValue);
            A.add(I);
            std::shared_ptr<Matrix<EigenValueType>> L, U;
            A.LUDecomposite(L, U);
            L->inverse();
            U->inverse();

            //Uu = L(-1)v
            Vector<EigenValueType> v(1, m_row);
            std::shared_ptr<Vector<EigenValueType>> u = U->mul(v);
            int loopCnt = 0;
            while (1)
            {
                ++loopCnt;
                int m = 0;
                u->magnitudeMaxItem(&m);
                EigenValueType um = (*u)[m];
                u->mul(1.0 / um);
                //迭代终止条件
                if (u->equal(v))
                {
                    break;
                }
                v = *u;
                u = L->mul(*u);
                u = U->mul(*u);
            }
            u->normalize();
            return u;
        }
        std::shared_ptr<Matrix<Complex>> eigen(std::vector<Complex>& eigenValues) const
        {
            eigenValues.resize(m_row);
            std::shared_ptr<Matrix> Q = std::make_shared<Matrix>(m_row, m_column);
            Q->identity();
            Matrix A = *this;
            std::shared_ptr<Vector<T>>      preDia;
            Vector<int>                     rootMask(1, m_row);
            std::vector<std::pair<int, int>> complexMask;
            int loopCnt = 0; 
            bool sysmetric = this->isSysmetric();
            while (1)
            {
                ++loopCnt;
                std::shared_ptr<Matrix> r;
                std::shared_ptr<Matrix> q;
                A.QRDecomposite(q, r);
                A = *r->mul(*q);
                /*A2 = Q1(-1)A1Q1
                A3 = Q2(-1)A2Q2 = Q2(-1)Q1(-1)A1Q1Q2=(Q1Q2)(-1)A1(Q1Q2)
                */
                *Q = *Q->mul(*q);

                //q -> I
                auto dia = q->diagonal();
                if (!preDia)
                {
                    preDia = dia;
                    continue;
                }
                //equal to last diagonal
                if (dia->equal(*preDia))
                {
                    break;
                }
                else
                {
                    complexMask.clear();
                    Vector<int> mask(m_row);
                    int i = 0;
                    for (; i < m_row; ++i)
                    {
                        if (mask[i] == 0)
                        {
                            //real root
                            if (nearZero(preDia->m_data[i] - dia->m_data[i]) && nearZero(abs(dia->m_data[i]) - 1))
                            {
                                mask[i] = 1;
                            }
                            //sysmetric matrix do not have complex eigen value
                            else if (!sysmetric)
                            {
                                int j = i + 1;
                                for (; j < m_row; ++j)
                                {
                                    //complex root
                                    if (nearZero(abs(dia->m_data[i]) - abs(dia->m_data[j])))
                                    {
                                        mask[i] = 2;
                                        mask[j] = 2;
                                        complexMask.push_back(std::make_pair(i, j));
                                        break;
                                    }
                                }
                                if (j == m_row)
                                {
                                    break;
                                }
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    if (i == m_row)
                    {
                        rootMask = mask;
                        break;
                    }
                }
                preDia = dia;
            }
            //compute complex eigen value
            for (int i = 0; i < complexMask.size(); ++i)
            {
                int m = complexMask[i].first;
                int n = complexMask[i].second;
                double sum = A.m_data[m][m] + A.m_data[n][n];
                double det = A.m_data[m][m] * A.m_data[n][n] - A.m_data[m][n] * A.m_data[n][m];
                //a b is complex,then a*b = det
                double realNum = sum / 2;
                double imagineNum = sqrt(det - realNum * realNum);
                eigenValues[m] = Complex(realNum, imagineNum);
                eigenValues[n] = Complex(realNum, -imagineNum);
            }
            auto eigenVectors = std::make_shared<Matrix<Complex>>(m_row, m_column);
            for (int i = 0; i < m_row; ++i)
            {
                if (!sysmetric)
                {
                    //eigen vectors of real eigen value
                    if (rootMask[i] == 1)
                    {
                        eigenValues[i] = Complex(A.m_data[i][i]);
                        auto realEigenVec = inversePowerMethod(A.m_data[i][i]);
                        eigenVectors->replaceColumn(*realEigenVec, i);
                    }
                    //eigen vectors of complex eigen value
                    else
                    {
                        auto complexEigenVec = inversePowerMethod(eigenValues[i]);
                        eigenVectors->replaceColumn(*complexEigenVec, i);
                    }
                }
                //sysmetric matrix
                else
                {
                    eigenValues[i] = Complex(A.m_data[i][i]);
                    eigenVectors->replaceColumn(*Q->column(i), i);
                }
            }
            return eigenVectors;
        }
    private:
        void householder(std::shared_ptr<Matrix>& Q, std::shared_ptr<Matrix>& R) const
        {
            //A=QR
            Q = std::make_shared<Matrix>(m_row, m_column);
            R = std::make_shared<Matrix>(m_row, m_column);
            Q->identity();
            *R = *this;
            for (int i = 0; i < m_row - 1; ++i)
            {
                auto v = R->column(i, i);
                float mag = v->magnitude();
                /*|y| = |x|
                v = (x-y)/|x-y|
                */
                //sign(u0) * |v|
                float norm = (*v)[0];
                T y = (*v)[0] / norm * mag;
                v->m_data[0] -= y;
                v->normalize();
                //h = I - 2VVt
                Matrix a(v->m_data, v->m_N, 1);
                auto b = a.transpose();
                Matrix I(m_row - i, m_row - i);
                I.identity().add(*(a.mul(-2).mul(*b)));

                auto Ai = R->subMatrix(m_row - i, m_column - i, i, i);
                Ai = I.mul(*Ai);

                auto h = I.expand(m_row, m_row, i, i);
                //h = h(i) * h(i-1)
                Q = h->mul(*Q);
                R->replace(*Ai, i, i);
            }
            Q = Q->transpose();
        }
    private:
        template <typename Element>
        Element transpose(int i, int j) const
        {
            return m_data[j][i];
        }
        template <typename Element>
        void inverse(Element)
        {
            Matrix I(m_row, m_column);
            I.identity();
            auto M = this->combine(I);
            M->simplest();
            *this = *(M->subMatrix(m_row, m_column, 0, m_column));
        }
        template <>
        Complex transpose(int i, int j) const
        {
            auto c = m_data[j][i];
            c.conjugate();
            return c;
        }
        template <>
        void inverse(Complex)
        {
            /*
            C = A + Bi
            C^ = (A + BA^B)^ - (A^B)(A + BA^B)^i
            */
            Matrix<float> A(m_row, m_column);
            Matrix<float> B(m_row, m_column);
            for (int i = 0; i < m_row; ++i)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    A[i][j] = m_data[i][j].m_real;
                    B[i][j] = m_data[i][j].m_imaginary;
                }
            }
            Matrix<float> invA = A;
            invA.inverse();
            auto A_plus_BinvAB = B.mul(invA)->mul(B)->add(A).inverse();
            auto invA_B_A_plus_BinvAB = invA.mul(B)->mul(A_plus_BinvAB);
            for (int i = 0; i < m_row; ++i)
            {
                for (int j = 0; j < m_column; ++j)
                {
                    m_data[i][j].m_real = A_plus_BinvAB[i][j];
                    m_data[i][j].m_imaginary = -(*invA_B_A_plus_BinvAB)[i][j];
                }
            }
        }
        template <typename Element>
        bool isSysmetric(int i, int j, Element a) const
        {
            if (nearZero(m_data[i][j] - m_data[j][i]))
            {
                return true;
            }
            return false;
        }
        template <>
        bool isSysmetric(int i, int j, Complex a) const
        {
            a = m_data[i][j] + m_data[j][i];
            auto b = m_data[i][j] - m_data[j][i];
            if (nearZero(a.m_imaginary) && nearZero(b.m_real))
            {
                return true;
            }
            return false;
        }
    private:
        T**     m_data;
        int     m_row;
        int     m_column;
    };
}