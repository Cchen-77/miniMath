#ifndef MINIMATH_H
#define MINIMATH_H
#include<vector>
#include<iostream>
#include<cmath>
namespace miniMath{
    typedef double Real;
}
namespace miniMath{
    class MatrixInitializer;
    class Matrix{
        friend class MatrixInitializer;
    public:
        uint32_t dim_0,dim_1;
        Matrix();
        Matrix(uint32_t dim_0,uint32_t dim_1);
        Matrix(const Matrix& rhs);
        ~Matrix();
    public:
        Real& at(uint32_t i,uint32_t j);
        Real at(uint32_t i,uint32_t j) const;
        Matrix operator+(const Matrix& rhs) const;
        Matrix operator-() const;
        Matrix operator-(const Matrix& rhs) const;
        Matrix operator*(const Matrix& rhs) const;
        Matrix& operator=(const Matrix& rhs);
        MatrixInitializer operator<<(Real item);
    public:
        void swapRow(uint32_t r1,uint32_t r2);
    public:
        static Matrix I(uint32_t dim);
    private:
        Real* items;
    };
    class Vector:public Matrix{
    public:
        Vector();
        Vector(uint32_t dim);
        Vector(const Matrix& matrix);
    public:
        using Matrix::at;
        Real& at(uint32_t i);
        Real at(uint32_t i) const;
    };

    class MatrixInitializer{
        friend class Matrix;
    private:
        MatrixInitializer(Matrix& mat,uint32_t idx);
        Matrix& matrix;
        uint32_t index;
    public:
        MatrixInitializer operator,(Real item);
    };

    std::ostream& operator<<(std::ostream& os,Matrix& mat);
    Matrix operator+(const Real& lhs,const Matrix& rhs);
    Matrix operator+(const Matrix& lhs,const Real& rhs);
    Matrix operator-(const Real& lhs,const Matrix& rhs);
    Matrix operator-(const Matrix& lhs,const Real& rhs);
    Matrix operator*(const Real& lhs,const Matrix& rhs);
    Matrix operator*(const Matrix& lhs,const Real& rhs);
    Matrix operator/(const Matrix& lhs,const Real& rhs);
}
namespace miniMath{
    Matrix transpose(const Matrix& mat);
    Real length(const Vector& vector);
    Real dot(const Vector& lhs,const Vector& rhs);
}
//Solvers:
namespace miniMath{
    void LUSplit(const Matrix& mat,Matrix& L,Matrix& U,Matrix& P,Matrix& invP);
    void QRSplit_MSchmidt(const Matrix& mat,Matrix& Q,Matrix& R);
    void QRSplit_HouseHolder(const Matrix& mat,Matrix& Q,Matrix& R);
    Vector JacobiSolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector GSSolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector SORSolve(const Matrix& A,const Vector& b,Real omega,bool printInfo=false);
    Vector GESolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector LUSolve(const Matrix& A,const Vector& b,bool printInfo=false);

}
namespace miniMath{
    class Polynomial{
    public:
        Polynomial();
        Polynomial(const Polynomial& other);
        Polynomial(std::initializer_list<Real> list);
    public:
        Polynomial operator+(const Polynomial& rhs) const;
        Polynomial operator*(const Polynomial& rhs) const;
        Polynomial operator*(Real rhs) const;
        Polynomial& operator=(const Polynomial& rhs);
        std::vector<Real> coeffs;
    public:
        void printCoeffs();
    };
    std::ostream& operator<<(std::ostream& os,const Polynomial& p);
    Polynomial operator*(Real lambda,const Polynomial& p);
}
//interpolation
namespace miniMath{
    Polynomial LagrangeInterpolation(std::vector<Real>& xi,std::vector<Real>& yi);
    std::vector<Polynomial> CubeSplineInterpolation(std::vector<Real>& xi,std::vector<Real>& yi);
}
#endif


#ifdef MINIMATH_IMPLEMENTATION
namespace miniMath{
    Matrix::Matrix():dim_0(0),dim_1(0),items(nullptr)
    {

    }

    Matrix::Matrix(uint32_t dim_0,uint32_t dim_1){
        this->dim_0 = dim_0;
        this->dim_1 = dim_1;
        items = new Real[dim_0*dim_1];
        memset(items,0,sizeof(Real)*dim_0*dim_1);
    }
    Matrix::Matrix(const Matrix& rhs){
        dim_0 = rhs.dim_0;
        dim_1 = rhs.dim_1;
        items = new Real[dim_0*dim_1];
        memset(items,0,sizeof(Real)*dim_0*dim_1);
        for(uint32_t i=0;i<dim_0;++i){
            for(uint32_t j=0;j<dim_1;++j){
                at(i,j) = rhs.at(i,j);
            }
        }
    }
    Matrix::~Matrix(){
        delete items;
    }
    Real& Matrix::at(uint32_t i,uint32_t j){
        if(i>=dim_0 || j>=dim_1){
            throw std::runtime_error("bad index!");
        }
        return items[i*dim_1+j];

    }
    Real Matrix::at(uint32_t i, uint32_t j) const
    {
        if(i>=dim_0 || j>=dim_1){
            throw std::runtime_error("bad index!");
        }
        return items[i*dim_1+j];
    }
    Matrix Matrix::operator+(const Matrix &rhs) const
    {
        if(rhs.dim_0!=dim_0 || rhs.dim_1!=dim_1){
            throw std::runtime_error("bad operation:add!");
        }
        Matrix mat(dim_0,dim_1);
        for(uint32_t i=0;i<dim_0;++i){
            for(uint32_t j=0;j<dim_1;++j){
                mat.at(i,j) = at(i,j) + rhs.at(i,j);
            }
        }
        return mat;

    }
    Matrix Matrix::operator-() const
    {
        Matrix mat(dim_0,dim_1);
        for(uint32_t i=0;i<dim_0;++i){
            for(uint32_t j=0;j<dim_1;++j){
                mat.at(i,j) = -at(i,j);
            }
        }
        return mat;

    }
    Matrix Matrix::operator-(const Matrix &rhs) const
    {
        return operator+(-rhs);
    }
    Matrix Matrix::operator*(const Matrix &rhs) const
    {
        if(dim_1 != rhs.dim_0){
            throw std::runtime_error("bad operator:mutiply!");
        }
        Matrix mat(dim_0,rhs.dim_1);
        for(uint32_t i=0;i<dim_0;++i){
            for(uint32_t j=0;j<rhs.dim_1;++j){
                for(uint32_t k=0;k<dim_1;++k){
                    mat.at(i,j) += at(i,k)*rhs.at(k,j);
                }

            }
        }
        return mat;
    }
    Matrix &Matrix::operator=(const Matrix &rhs)
    {
        delete items;
        dim_0 = rhs.dim_0;
        dim_1 = rhs.dim_1;
        items = new Real[dim_0*dim_1];
        memset(items,0,sizeof(Real)*dim_0*dim_1);
        for(uint32_t i=0;i<dim_0;++i){
            for(uint32_t j=0;j<dim_1;++j){
                at(i,j) = rhs.at(i,j);
            }
        }
        return *this;
    }

    MatrixInitializer Matrix::operator<<(Real item)
    {
        items[0] = item;
        return MatrixInitializer(*this,1);
    }

    void Matrix::swapRow(uint32_t r1, uint32_t r2)
    {
        for(uint32_t c=0;c<dim_1;++c){
            std::swap(at(r1,c),at(r2,c));
        }
    }

    Matrix Matrix::I(uint32_t dim)
    {
        Matrix matI = Matrix(dim,dim);
        for(uint32_t i=0;i<dim;++i)
            matI.at(i,i) = 1;
        return matI;
    }

    Vector::Vector():Matrix()
    {
        
    }

    Vector::Vector(uint32_t dim) : Matrix(dim, 1)
    {
    }

    Vector::Vector(const Matrix &matrix):Matrix(matrix)
    {
        if(matrix.dim_1 != 1){
            throw std::runtime_error("bad vector init with matrix!");
        }
        
    }
    Real &Vector::at(uint32_t i)
    {
        return Matrix::at(i,0);
    }
    Real Vector::at(uint32_t i) const
    {
        return Matrix::at(i,0);
    }
    std::ostream &operator<<(std::ostream &os, Matrix &mat)
    {
        for(uint32_t i=0;i<mat.dim_0;++i){
            for(uint32_t j=0;j<mat.dim_1;++j){
                os<<mat.at(i,j)<<",";
            }
            os<<std::endl;
        }
        return os;
    }

    Matrix operator+(const Real &lhs, const Matrix &rhs)
    {
        Matrix result = rhs;
        for(int i=0;i<rhs.dim_0;++i){
            for(int j=0;j<rhs.dim_1;++j){
                result.at(i,j)+=lhs;
            }
        }
        return result;
    }
    Matrix operator+(const Matrix &lhs, const Real &rhs)
    {
        Matrix result = lhs;
        for(int i=0;i<lhs.dim_0;++i){
            for(int j=0;j<lhs.dim_1;++j){
                result.at(i,j)+=rhs;
            }
        }
        return result;
    }
    Matrix operator-(const Real &lhs, const Matrix &rhs)
    {
        Matrix result = rhs;
        for(int i=0;i<rhs.dim_0;++i){
            for(int j=0;j<rhs.dim_1;++j){
                result.at(i,j)-=lhs;
            }
        }
        return result;
    }
    Matrix operator-(const Matrix &lhs, const Real &rhs)
    {
        Matrix result = lhs;
        for(int i=0;i<lhs.dim_0;++i){
            for(int j=0;j<lhs.dim_1;++j){
                result.at(i,j)-=rhs;
            }
        }
        return result;
    }
    Matrix operator*(const Real &lhs, const Matrix &rhs)
    {
        Matrix result = rhs;
        for(int i=0;i<rhs.dim_0;++i){
            for(int j=0;j<rhs.dim_1;++j){
                result.at(i,j)*=lhs;
            }
        }
        return result;
    }
    Matrix operator*(const Matrix &lhs, const Real &rhs)
    {
        Matrix result = lhs;
        for(int i=0;i<lhs.dim_0;++i){
            for(int j=0;j<lhs.dim_1;++j){
                result.at(i,j)*=rhs;
            }
        }
        return result;
    }
    Matrix operator/(const Matrix &lhs, const Real &rhs)
    {
        if(rhs==0){
            throw std::runtime_error("matrix divide by zero!");
        }
        Matrix result = lhs;
        for(int i=0;i<lhs.dim_0;++i){
            for(int j=0;j<lhs.dim_1;++j){
                result.at(i,j)/=rhs;
            }
        }
        return result;
    }

    MatrixInitializer::MatrixInitializer(Matrix& mat,uint32_t idx):matrix(mat),index(idx){

    }
    MatrixInitializer MatrixInitializer::operator,(Real item)
    {
        if(index >= matrix.dim_0*matrix.dim_1){
            throw std::runtime_error("fail to do matrix(vector) comma init!");
        }
        matrix.items[index] = item;       
        return MatrixInitializer(matrix,index+1);
    }
}
namespace miniMath{
    Matrix transpose(const Matrix &mat)
    {
        Matrix result(mat.dim_1,mat.dim_0);
        for(int i=0;i<mat.dim_0;++i){
            for(int j=0;j<mat.dim_1;++j){
                result.at(j,i) = mat.at(i,j);
            }
        }
        return result;
    }

    Real length(const Vector& vector){
        if(vector.dim_1!=1){
            throw std::runtime_error("bad vector:maybe you should not give a matrix!");
        }
        Real len = 0;
        for(uint32_t i=0;i<vector.dim_0;++i){
            len += vector.at(i)*vector.at(i);
        }
        return sqrt(len);
    }
    Real dot(const Vector &lhs, const Vector &rhs)
    {
        if(lhs.dim_0 != rhs.dim_0){
            throw std::runtime_error("lhs vector should have same length with rhs vector!");
        }
        Real result = 0;
        for(uint32_t i=0;i<lhs.dim_0;++i){
            result += lhs.at(i)*rhs.at(i);
        }
        return result;
    }
}
namespace miniMath{
    void QRSplit_MSchmidt(const Matrix& mat,Matrix& Q,Matrix& R){
        int N = mat.dim_1;
        int M = mat.dim_0;
        std::vector<Vector> a(N);
        std::vector<Vector> e(N);
        for(int i=0;i<N;++i){
            a[i] = Vector(M);
            for(int j=0;j<M;++j){
                a[i].at(j) = mat.at(j,i);
            }
        }
        R = Matrix(N,N);
        for(int i=0;i<N;++i){
            Real rii = length(a[i]);
            R.at(i,i) = rii;
            if(rii!=0){
                e[i] = a[i]/rii;
            }
            else{
                e[i] = a[i];
            }
            for(int j=1;j<N;++j){
                Real rij = dot(e[i],a[j]);
                a[j] = a[j] - rij*e[i];
                R.at(i,j) = rij;
            }
        }
        Q = Matrix(M,N);
        for(int i=0;i<N;++i){
            for(int j=0;j<M;++j){
                Q.at(j,i) = e[i].at(j);
            }
        }
    }
    void QRSplit_HouseHolder(const Matrix& mat,Matrix& Q,Matrix& R){

    }
    void LUSplit(const Matrix& mat,Matrix& L,Matrix& U,Matrix& P,Matrix& invP){
        if(mat.dim_0!=mat.dim_1){
            throw std::runtime_error("bad matrix A or bad vector b!");
        }
        if(mat.dim_0>100){
            throw std::runtime_error("matrix too big!");
        }
        uint32_t numRows = mat.dim_0;
        uint32_t numColumns = mat.dim_1;

        Matrix tempMat = mat;
        L = Matrix::I(numRows);
        P = Matrix::I(numRows);
        invP = Matrix::I(numRows);
        uint32_t finishrow_count = 0;

        Real pivot = std::abs(tempMat.at(0,0));
        uint32_t pivot_row = 0;
        Matrix invPk = Matrix::I(numRows);
        for(uint32_t i=1;i<numRows;++i){
            if(std::abs(tempMat.at(i,0)) > pivot){
                pivot = std::abs(tempMat.at(i,0));
                pivot_row = i;
            }
        }
        tempMat.swapRow(pivot_row,finishrow_count);
        invPk.swapRow(pivot_row,finishrow_count);
        invP = invPk*invP;
        P = P*invPk;
        ++finishrow_count;
    
        while(finishrow_count < numRows){
            if(pivot == 0){
                throw std::runtime_error("Bad matrix with more than one solution!");
            }
            Matrix Lk = Matrix::I(numRows);
            uint32_t pivot_column = finishrow_count-1;
            for(uint32_t row=finishrow_count;row<numRows;++row){
                Real scale = tempMat.at(row,pivot_column)/pivot;
                Lk.at(row,pivot_column) += scale;
                tempMat.at(row,pivot_column) = 0;
                for(uint32_t column=pivot_column+1;column<numColumns;++column){
                    tempMat.at(row,column) -= scale*tempMat.at(pivot_row,column);
                }
            }
            L = L*Lk;
            Matrix invPk = Matrix::I(numRows);
            pivot = 0;
            for(uint32_t row=finishrow_count;row<numRows;++row){
                if(std::abs(tempMat.at(row,pivot_column+1)) > pivot){
                    pivot = std::abs(tempMat.at(row,pivot_column+1));
                    pivot_row = row;
                }
            }
            tempMat.swapRow(pivot_row,finishrow_count);
            invPk.swapRow(pivot_row,finishrow_count);
            invP = invPk*invP;
            P = P*invPk;
            ++finishrow_count;
        }
        if(pivot == 0){
            throw std::runtime_error("Bad matrix could not be LU splited!");
        }

        U = tempMat;
       
    }

    Vector JacobiSolve(const Matrix& A,const Vector& b,bool printInfo){
        //you should init A for aii !=0
        if(A.dim_0!=A.dim_1 || b.dim_0 != A.dim_0){
            throw std::runtime_error("bad matrix A or bad vector b!");
        }
        uint32_t numRows = A.dim_0;
        uint32_t numColumns = A.dim_1;
        Vector x[2];
        uint32_t curframe = 0;
        for(uint32_t i=0;i<2;++i){
            x[i] = Vector(numRows);
        }
        uint32_t lastframe;
        uint32_t iter=0;
        do{
            if(printInfo){
                printf("iteration %d:\n",iter++);
                std::cout<<x[curframe];
            }
            lastframe = curframe;
            curframe^=1;
            for(uint32_t row=0;row<numRows;++row){
                Real t = b.at(row);
                for(uint32_t column=0;column<numColumns;++column){
                    if(column != row){
                        t -= x[lastframe].at(column)*A.at(row,column);
                    }
                }
                x[curframe].at(row) = t/A.at(row,row);
            }
        }while(length(x[curframe]-x[lastframe])>5e-5);
        if(printInfo){
            printf("total iter:%d\n",iter);
        }
        return x[curframe];
    }
    Vector GSSolve(const Matrix& A,const Vector& b,bool printInfo){
        //you should init A for aii !=0
        if(A.dim_0!=A.dim_1 || b.dim_0 != A.dim_0){
            throw std::runtime_error("bad matrix A or bad vector b!");
        }
        uint32_t numRows = A.dim_0;
        uint32_t numColumns = A.dim_1;
        Vector x[2];
        uint32_t curframe = 0;
        for(uint32_t i=0;i<2;++i){
            x[i] = Vector(numRows);
        }
        uint32_t lastframe;
        uint32_t iter=0;
        do{
            if(printInfo){
                printf("iteration %d:\n",iter++);
                std::cout<<x[curframe];
            }
            lastframe = curframe;
            curframe^=1;
            for(uint32_t row=0;row<numRows;++row){
                Real t = b.at(row);
                for(uint32_t column=0;column<row;++column){
                    t -= x[curframe].at(column)*A.at(row,column);
                }
                for(uint32_t column=row+1;column<numRows;++column){
                    t -= x[lastframe].at(column)*A.at(row,column);
                }
                x[curframe].at(row) = t/A.at(row,row);
            }
        }while(length(x[curframe]-x[lastframe])>5e-5);
        if(printInfo){
            printf("total iter:%d\n",iter);
        }
        return x[curframe];
    }
    Vector SORSolve(const Matrix& A,const Vector& b,Real omega,bool printInfo){
        //you should init A for aii !=0
        if(A.dim_0!=A.dim_1 || b.dim_0 != A.dim_0){
            throw std::runtime_error("bad matrix A or bad vector b!");
        }
        uint32_t numRows = A.dim_0;
        uint32_t numColumns = A.dim_1;
        Vector x[2];
        uint32_t curframe = 0;
        for(uint32_t i=0;i<2;++i){
            x[i] = Vector(numRows);
        }
        uint32_t lastframe;
        uint32_t iter=0;
        do{
            if(printInfo){
                printf("iteration %d:\n",iter++);
                std::cout<<x[curframe];
            }
            lastframe = curframe;
            curframe^=1;
            for(uint32_t row=0;row<numRows;++row){
                Real t = b.at(row);
                for(uint32_t column=0;column<row;++column){
                    t -= x[curframe].at(column)*A.at(row,column);
                }
                for(uint32_t column=row+1;column<numRows;++column){
                    t -= x[lastframe].at(column)*A.at(row,column);
                }
                x[curframe].at(row) = omega*t/A.at(row,row) + (1-omega)*x[lastframe].at(row);
            }
        }while(length(x[curframe]-x[lastframe])>5e-5);
        if(printInfo){
            printf("total iter:%d\n",iter);
        }
        return x[curframe];
    }
    Vector GESolve(const Matrix& A,const Vector& b,bool printInfo){
        if(A.dim_0!=A.dim_1 || A.dim_0 != b.dim_0){
            throw std::runtime_error("bad matrix A or bad vector b!");
        }
        if(A.dim_0>100){
            throw std::runtime_error("matrix too big!");
        }
        uint32_t numRows = A.dim_0;
        uint32_t numColumns = A.dim_1;
        Matrix tempA = A;
        Vector tempb = b;

        uint32_t finishrow_count = 0;

        Real pivot = std::abs(tempA.at(0,0));
        uint32_t pivot_row = 0;

        for(uint32_t i=1;i<numRows;++i){
            if(std::abs(tempA.at(i,0)) > pivot){
                pivot = std::abs(tempA.at(i,0));
                pivot_row = i;
            }
        }
        tempA.swapRow(pivot_row,finishrow_count);
        tempb.swapRow(pivot_row,finishrow_count);
        ++finishrow_count;
   
        while(finishrow_count < numRows){
            if(pivot == 0){
                throw std::runtime_error("Bad matrix with more than one solution!");
            }
            uint32_t pivot_column = finishrow_count-1;
            for(uint32_t row=finishrow_count;row<numRows;++row){
               
                Real scale = tempA.at(row,pivot_column)/pivot;
                tempA.at(row,pivot_column) = 0;
                for(uint32_t column=pivot_column+1;column<numColumns;++column){
                    tempA.at(row,column) -= scale*tempA.at(pivot_row,column);
                }
                tempb.at(row) -= scale*tempb.at(pivot_row);

            }
            pivot = 0;
            for(uint32_t row=finishrow_count;row<numRows;++row){
                if(std::abs(tempA.at(row,pivot_column+1)) > pivot){
                    pivot = std::abs(tempA.at(row,pivot_column+1));
                    pivot_row = row;
                }
            }
            tempA.swapRow(pivot_row,finishrow_count);
            tempb.swapRow(pivot_row,finishrow_count);
            ++finishrow_count;
        }
        if(pivot == 0){
            throw std::runtime_error("Bad matrix with more than one solution!");
        }
        Vector solution(tempb.dim_0);
        for(int row=numRows-1;row>=0;--row){
            Real solution_temp = tempb.at(row);
            for(int column = numColumns -1;column>row;--column){
                solution_temp -= tempA.at(row,column)*solution.at(column);
            }
            solution_temp /= tempA.at(row,row);
            solution.at(row) = solution_temp;
        }
        return solution;
    }
    Vector LUSolve(const Matrix& A,const Vector& b,bool printInfo){
        if(A.dim_0!=A.dim_1 || A.dim_0 != b.dim_0){
            throw std::runtime_error("bad matrix A or bad vector b!");
        }
        if(A.dim_0>100){
            throw std::runtime_error("matrix too big!");
        }
        Matrix L,U,P,invP;
        LUSplit(A,L,U,P,invP);
        if(printInfo){
            printf("===========infos===============\n");
            printf("L Matrix(%dx%d):\n",L.dim_0,L.dim_1);
            std::cout<<L;
            printf("U Matrix(%dx%d):\n",U.dim_0,U.dim_1);
            std::cout<<U;
            printf("invP Matrix(%dx%d):\n",invP.dim_0,invP.dim_1);
            std::cout<<invP;
            printf("===============================\n");
        }
        Vector tempb = invP*b;
        uint32_t numRows = A.dim_0;
        uint32_t numColumns = A.dim_1;
        Vector solution_u(numRows);
        for(uint32_t row=0;row<numRows;++row){
            Real t = tempb.at(row);
            for(uint32_t column=0;column<row;++column){
                t-=L.at(row,column)*solution_u.at(column);
            }
            t/=L.at(row,row);
            solution_u.at(row) = t;
        }
        Vector solution(numRows);
        for(int row=numRows-1;row>=0;--row){
            Real solution_temp = solution_u.at(row);
            for(int column = numColumns -1;column>row;--column){
                solution_temp -= U.at(row,column)*solution.at(column);
            }
            solution_temp /= U.at(row,row);
            solution.at(row) = solution_temp;
        }
        return P*solution;
    }
}
namespace miniMath{
    Polynomial::Polynomial(){
        
    }
    Polynomial::Polynomial(const Polynomial &other)
    {
        coeffs.resize(other.coeffs.size());
        for(int i=0;i<coeffs.size();++i){
            coeffs[i] = other.coeffs[i];
        }
    }
    Polynomial::Polynomial(std::initializer_list<Real> list)
    {
        if(list.begin()==list.end()){
            throw std::runtime_error("empty polynomial is not allowed!plz use polynomial()");
        }
        for(auto it = list.begin();it!=list.end();++it){
            coeffs.push_back(*it);
        }
    }
    Polynomial Polynomial::operator+(const Polynomial& rhs) const{
        int maxSize = std::max(coeffs.size(),rhs.coeffs.size());
        Polynomial result;
        for(int i=0;i<maxSize;++i){
            Real sum = 0.0;
            if(i<coeffs.size()){
                sum += coeffs[i];
            }
            if(i<rhs.coeffs.size()){
                sum+=rhs.coeffs[i];
            }
            result.coeffs.push_back(sum);
        }
        return result;
    }
    Polynomial Polynomial::operator*(const Polynomial& rhs) const{
        int resultSize = coeffs.size()+rhs.coeffs.size()-1;
        Polynomial result;
        if(coeffs.size()==0||rhs.coeffs.size()==0){
            return result;
        }
        result.coeffs.resize(resultSize,0);
        for(int i=0;i<coeffs.size();++i){
            for(int j=0;j<rhs.coeffs.size();++j){
                int dst = i+j;
                result.coeffs[dst] += coeffs[i]*rhs.coeffs[j];
            }
        }
        return result;
    }
    Polynomial Polynomial::operator*(Real rhs) const{
        Polynomial result;
        result = *this;
        for(int i=0;i<result.coeffs.size();++i){
            result.coeffs[i]*=rhs;
        }
        return result;
    }
    Polynomial operator*(Real lambda, const Polynomial &p)
    {
        return p*lambda;
    }
    Polynomial &Polynomial::operator=(const Polynomial &rhs)
    {
        coeffs.resize(rhs.coeffs.size());
        for(int i=0;i<coeffs.size();++i){
            coeffs[i] = rhs.coeffs[i];
        }
        return *this;
    }
    void Polynomial::printCoeffs()
    {
        int count = 0;
        for(int i=0;i<coeffs.size();++i){
            std::cout<<coeffs[i]<<",";
             count=(count+1)%5;
            if(count==0){
                std::cout<<'\n';
            }
        }
        std::cout<<'\n';
    }
    std::ostream &operator<<(std::ostream &os, const Polynomial &p)
    {
        int count = 1;
        os<<p.coeffs[0];
        for(int i=1;i<p.coeffs.size();++i){
            os<<"+"<<p.coeffs[i]<<"x^"<<i;
            count=(count+1)%5;
            if(count==0){
                os<<'\n';
            }
        }
        os<<'\n';
        return os;
    }

}
namespace miniMath{
    Polynomial LagrangeInterpolation(std::vector<Real>& xi,std::vector<Real>& yi){
        if(xi.size()!=yi.size()){
            throw std::runtime_error("xi and yi should be in same size!");
        }
        if(xi.size()<2){
            std::runtime_error("xi/yi size should be greater than 1!");
        }
        int degree = xi.size()-1;
        Polynomial result;
        for(int i=0;i<=degree;++i){
            Polynomial t{yi[i]};
            for(int j=0;j<=degree;++j){
                if(i!=j){
                    t = t*Polynomial{-xi[j]/(xi[i]-xi[j]),1.0/(xi[i]-xi[j])};
                }
            }
            result = result + t;
        }
        return result;
    }
    std::vector<Polynomial> CubeSplineInterpolation(std::vector<Real>& xi,std::vector<Real>& yi){
        
        if(xi.size()!=yi.size()){
            throw std::runtime_error("xi and yi should be in same size!");
        }
        if(xi.size()<2){
            std::runtime_error("xi/yi size should be greater than 1!");
        }
        int N = xi.size()-1;
        Real df0 = 0.0;
        Real dfn = 0.0;
        std::vector<Real> hi;
        std::vector<Real> di;
        hi.resize(N+1,0.0);
        di.resize(N+1,0.0);
        for(int i=1;i<=N;++i){
            hi[i] = xi[i] - xi[i-1];
        }
        //we use natural spline
        di[0] = 0.0;
        di[N] = 0.0;
        for(int i=1;i<N;++i){
            di[i] = 6.0*((yi[i+1]-yi[i])/hi[i+1] - (yi[i]-yi[i-1])/hi[i])/(hi[i]+hi[i+1]);
        }
        Vector d(N+1);
        for(int i=0;i<=N;++i){
            d.at(i) = di[i];
        }
        std::vector<Real> mii;
        std::vector<Real> lambdai;
        mii.resize(N+1,0.0);
        lambdai.resize(N+1,0.0);
        for(int i=0;i<N;++i){
            lambdai[i] = hi[i+1]/(hi[i]+hi[i+1]);
            mii[i] = 1.0 - lambdai[i];
        }
        Matrix A(N+1,N+1);
        A.at(0,0) = 1;
        for(int i=1;i<N;++i){
            A.at(i,i-1) = mii[i];
            A.at(i,i) = 2;
            A.at(i,i+1) = lambdai[i];
        }
        A.at(N,N) = 1;
        
        //Gauss-Seidel convergence with Dominant main diagonal
        Vector m = GSSolve(A,d);
        std::vector<Polynomial> Polynomials; 
        for(int i=1;i<=N;++i){
            Polynomial Si;
            Polynomial t1{xi[i],-1.0};
            Polynomial t2{-xi[i-1],1.0};
            Si = (m.at(i-1)/(6.0*hi[i]))*(t1*t1*t1) + (m.at(i)/(6.0*hi[i]))*(t2*t2*t2)
                + (yi[i-1]-m.at(i-1)/6.0*hi[i]*hi[i])/hi[i]*t1 + (yi[i]-m.at(i)/6.0*hi[i]*hi[i])/hi[i]*t2;
            Polynomials.push_back(Si);
        }
        return Polynomials;
    }
}
#endif