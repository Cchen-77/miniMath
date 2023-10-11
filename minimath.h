#ifndef MINIMATH_H
#define MINIMATH_H
#include<vector>
#include<iostream>
#include<cmath>
namespace miniMath{
    double calcRoot(double x,int t);
}
namespace miniMath{
    class MatrixInitializer;
    class Matrix{
        friend class MatrixInitializer;
    public:
        uint32_t dim_0,dim_1;
        Matrix(uint32_t dim_0,uint32_t dim_1);
        Matrix(const Matrix& rhs);
        ~Matrix();
    public:
        float& at(uint32_t i,uint32_t j);
        float at(uint32_t i,uint32_t j) const;
        Matrix operator+(const Matrix& rhs) const;
        Matrix operator-() const;
        Matrix operator-(const Matrix& rhs) const;
        Matrix operator*(const Matrix& rhs) const;
        Matrix& operator=(const Matrix& rhs);
        MatrixInitializer operator<<(float item);

    private:
        float* items;
    };
    class Vector:public Matrix{
    public:
        Vector(uint32_t dim);
        float& at(uint32_t i);
    };

    class MatrixInitializer{
        friend class Matrix;
    private:
        MatrixInitializer(Matrix& mat,uint32_t idx);
        Matrix& matrix;
        uint32_t index;
    public:
        MatrixInitializer operator,(float item);
    };

    std::ostream& operator<<(std::ostream& os,Matrix& mat);
}
namespace miniMath{
    void LUSplit(const Matrix& mat,Matrix& L,Matrix& U);
    float det(const Matrix& mat);
    Matrix inverse(const Matrix& mat);

    Vector JacobiSolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector GSSolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector SORSolve(const Matrix& A,const Vector& b,float omega,bool printInfo=false);
    Vector GESolve(const Matrix& A,const Vector& b,bool printInfo=false);
    Vector LUSolve(const Matrix& A,const Vector& b,bool printInfo=false);
    
    Vector LSolve(const Matrix& A,const Vector& b);
    Vector USovle(const Matrix& A,const Vector& b);
}
#endif
#define MINIMATH_IMPLEMENTATION
#ifdef MINIMATH_IMPLEMENTATION
namespace miniMath{

    //HelperFuncs Implementation
    double calcRoot(double x,int t){
        if(x<0&&t%2==0){
            throw std::runtime_error("bad input!");
        }
        double low,high=1e9;
        if(t%2==0){
            low = 0;
        }
        else{
            low = -1e9;
        }
        while(high-low >= 1e-10){ 
            double mid = (low+high)/2;
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
    Matrix::Matrix(uint32_t dim_0,uint32_t dim_1){
        this->dim_0 = dim_0;
        this->dim_1 = dim_1;
        items = new float[dim_0*dim_1];
        memset(items,0,sizeof(float)*dim_0*dim_1);
    }
    Matrix::Matrix(const Matrix& rhs){
        dim_0 = rhs.dim_0;
        dim_1 = rhs.dim_1;
        items = new float[dim_0*dim_1];
        memset(items,0,sizeof(float)*dim_0*dim_1);
        for(uint32_t i=0;i<dim_0;++i){
            for(uint32_t j=0;j<dim_1;++j){
                at(i,j) = rhs.at(i,j);
            }
        }
    }
    Matrix::~Matrix(){
        delete items;
    }
    float& Matrix::at(uint32_t i,uint32_t j){
        if(i>=dim_0 || j>=dim_1){
            throw std::runtime_error("bad index!");
        }
        return items[i*dim_1+j];

    }
    float Matrix::at(uint32_t i, uint32_t j) const
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
        items = new float[dim_0*dim_1];
        memset(items,0,sizeof(float)*dim_0*dim_1);
        for(uint32_t i=0;i<dim_0;++i){
            for(uint32_t j=0;j<dim_1;++j){
                at(i,j) = rhs.at(i,j);
            }
        }
        return *this;
    }

    MatrixInitializer Matrix::operator<<(float item)
    {
        items[0] = item;
        return MatrixInitializer(*this,1);
    }

    Vector::Vector(uint32_t dim):Matrix(dim,1){

    }

    float &Vector::at(uint32_t i)
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
    MatrixInitializer MatrixInitializer::operator,(float item)
    {
        if(index >= matrix.dim_0*matrix.dim_1){
            throw std::runtime_error("fail to do matrix(vector) comma init!");
        }
        matrix.items[index] = item;       
        return MatrixInitializer(matrix,index+1);
    }
}

namespace miniMath{
    void LUSplit(const Matrix& mat,Matrix& L,Matrix& U){
        return;
    }
    float det(const Matrix& mat){
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
    Vector SORSolve(const Matrix& A,const Vector& b,float omega,bool printInfo){
        return b;
    }
    Vector GESolve(const Matrix& A,const Vector& b,bool printInfo){
        if(A.dim_0!=A.dim_1 || A.dim_0 != b.dim_0){
            throw std::runtime_error("bad matrix A and bad vector b!");
        }
        if(A.dim_0>100){
            throw std::runtime_error("matrix too big!");
        }
        uint32_t numRows = A.dim_0;
        uint32_t numColumns = A.dim_1;
        Matrix tempA = A;
        Vector tempb = b;

        int row_picked[101];
        memset(row_picked,0,sizeof(row_picked));
        std::vector<uint32_t> picked_rows;

        float pivot = std::abs(tempA.at(0,0));
        uint32_t pivot_row = 0;

        for(uint32_t i=1;i<numRows;++i){
            if(std::abs(tempA.at(i,0)) > pivot){
                pivot = std::abs(tempA.at(i,0));
                pivot_row = i;
            }
        }

        picked_rows.push_back(pivot_row);
        row_picked[pivot_row] = 1;
        while(picked_rows.size() < numRows){
            if(pivot == 0){
                throw std::runtime_error("Bad matrix with more than one solution!");
            }
            uint32_t pivot_column = picked_rows.size()-1;
            for(uint32_t row=0;row<numRows;++row){
                if(!row_picked[row]){
                    float scale = tempA.at(row,pivot_column)/pivot;
                    tempA.at(row,pivot_column) = 0;
                    
                    for(uint32_t column=pivot_column+1;column<numColumns;++column){
                        tempA.at(row,column) -= scale*tempA.at(pivot_row,column);
                    }
                    tempb.at(row) -= scale*tempb.at(pivot_row);
                }
            }
            pivot = 0;
            for(uint32_t row=0;row<numRows;++row){
                if(!row_picked[row]&&std::abs(tempA.at(row,pivot_column+1)) > pivot){
                    pivot = std::abs(tempA.at(row,pivot_column+1));
                    pivot_row = row;
                }
            }
            picked_rows.push_back(pivot_row);
            row_picked[pivot_row] = 1;
        }
        if(pivot == 0){
            throw std::runtime_error("Bad matrix with more than one solution!");
        }
        Vector solution(tempb.dim_0);
        for(int row=numRows-1;row>=0;--row){
            float solution_temp = tempb.at(picked_rows[row]);
            for(int column = numColumns -1;column>row;--column){
                solution_temp -= tempA.at(picked_rows[row],column)*solution.at(picked_rows[column]);
            }
            solution_temp /= tempA.at(picked_rows[row],row);
            solution.at(picked_rows[row]) = solution_temp;
        }
        return solution;
    }
    Vector LUSolve(const Matrix& A,const Vector& b,bool printInfo){
        return b;
    }
    
    Vector LSolve(const Matrix& A,const Vector& b){
        return b;
    }
    Vector USovle(const Matrix& A,const Vector& b){
        return b;
    }
}
#endif