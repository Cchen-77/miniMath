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
    bool LUSplit(const Matrix& mat,Matrix& L,Matrix& U);
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
#endif