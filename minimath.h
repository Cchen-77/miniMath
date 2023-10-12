#ifndef MINIMATH_H
#define MINIMATH_H
#include<vector>
#include<iostream>
#include<cmath>
namespace miniMath{
    typedef double Real;
}
namespace miniMath{
    Real calcRoot(Real x,int t);
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
}
namespace miniMath{
    void LUSplit(const Matrix& mat,Matrix& L,Matrix& U,Matrix& P,Matrix& invP);
    Real det(const Matrix& mat);
    Matrix inverse(const Matrix& mat);

    Vector JacobiSolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector GSSolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector SORSolve(const Matrix& A,const Vector& b,Real omega,bool printInfo=false);
    Vector GESolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector LUSolve(const Matrix& A,const Vector& b,bool printInfo=false);

}
#endif
#define MINIMATH_IMPLEMENTATION
#ifdef MINIMATH_IMPLEMENTATION
namespace miniMath{

    //HelperFuncs Implementation
    Real calcRoot(Real x,int t){
        if(x<0&&t%2==0){
            throw std::runtime_error("bad input!");
        }
        Real low,high=1e9;
        if(t%2==0){
            low = 0;
        }
        else{
            low = -1e9;
        }
        while(high-low >= 1e-10){ 
            Real mid = (low+high)/2;
            if(std::pow(mid,t)>x){
                high = mid;
            }
            else{
                low = mid;
            }
        }
        return low;
    }
}

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

    std::ostream& operator<<(std::ostream& os,Matrix& mat){
        for(uint32_t i=0;i<mat.dim_0;++i){
            for(uint32_t j=0;j<mat.dim_1;++j){
                os<<mat.at(i,j)<<",";
            }
            os<<std::endl;
        }
        return os;
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
    Real det(const Matrix& mat){
        return 0.0f;
    }
    Matrix inverse(const Matrix& mat){
        return mat;
    }

    Vector JacobiSolve(const Matrix& A,const Vector& b,bool printInfo){
        return b;
    }
    Vector GSSolve(const Matrix& A,const Vector& b,bool printInfo){
        return b;
    }
    Vector SORSolve(const Matrix& A,const Vector& b,Real omega,bool printInfo){
        return b;
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
#endif