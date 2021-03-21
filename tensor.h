#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <math.h>
#include <iomanip>  


/**
 * @brief 
 * 
 * @tparam T 
 */
template <class T>
class tensor2{
    public:


        // ------------------------------------------
        //              CONSTRUCTOR
        // ------------------------------------------
        tensor2() = default;
        /**
         * @brief Construct a new tensor2 object
         * 
         * @param px value1
         * @param py value2
         */
        tensor2(const T& px, const T& py): __px(px), __py(py){}
        tensor2(const tensor2<T>& other);
        tensor2(tensor2<T>&& other);
        ~tensor2(){}


        // -------------------------------------------
        //        ARITHMETIC OPERATOR OVERLOAD
        // -------------------------------------------
        tensor2<T>& operator=(const tensor2<T>& other);
        tensor2<T> operator+(const tensor2<T>& other);
        tensor2<T>& operator+=(const tensor2<T>& other);
        tensor2<T> operator-(const tensor2<T>& other);
        tensor2<T>& operator-=(const tensor2<T>& other);
        tensor2<T> operator/(const double& divisor);
        tensor2<T> operator*(const double& multiplier);
        tensor2<T>& operator/=(const double& divisor);
        tensor2<T>& operator*=(const double& multiplier);

        // -------------------------------------------
        //               TENSOR OPERATIONS
        // -------------------------------------------
        T norm();
        T manhattan();

         // -------------------------------------------
        //             GETTERS & SETTERS
        // -------------------------------------------
        T x1(){return __px;}
        T x2(){return __py;}

        void setx1(const T& x1){__px = x1;}
        void setx2(const T& x2){__py = x2;}



        // -------------------------------------------
        //               REPR
        // -------------------------------------------
        friend std::ostream& operator<<(std::ostream& os, tensor2<T>& t){
            os<<"("<<std::setprecision(6)<<t.__px<<", "<<std::setprecision(6)<<t.__py<<")";
            return os;}



    private:

        T __px = 0.0;//value1
        T __py = 0.0;//value2






};



// -------------------------------------------
//        ARITHMETIC OPERATOR OVERLOAD
// -------------------------------------------
template <class T>
tensor2<T>& tensor2<T>::operator=(const tensor2<T>& other){
    __px = other.__px;
    __py = other.__py; 
    return *this;
}

template <class T>
tensor2<T>& tensor2<T>::operator+=(const tensor2<T>& other){
    __px += other.__px;
    __py += other.__py; 
    return *this;
}

template <class T>
tensor2<T>& tensor2<T>::operator-=(const tensor2<T>& other){
    __px -= other.__px;
    __py -= other.__py; 
    return *this;
}

template <class T>
tensor2<T> tensor2<T>::operator+(const tensor2<T>& other){
    tensor2<T> newpt(__px+other.__px,__py+other.__py);
    return std::move(newpt);
}

template <class T>
tensor2<T> tensor2<T>::operator-(const tensor2<T>& other){
    tensor2<T> newpt(__px-other.__px,__py-other.__py);
    return std::move(newpt);
}
template <class T>
tensor2<T> tensor2<T>::operator/(const double& divisor){
    tensor2<T> newpt(__px/divisor,__py/divisor);
    return std::move(newpt);
}

template <class T>
tensor2<T> tensor2<T>::operator*(const double& multiplier){
    tensor2<T> newpt(__px*multiplier,__py*multiplier);
    return std::move(newpt);
}


template <class T>
tensor2<T>& tensor2<T>::operator/=(const double& divisor){
    __px/=divisor;
    __py/=divisor;
    return  *this;
}


template <class T>
tensor2<T>& tensor2<T>::operator*=(const double& multiplier){
    __px*=multiplier;
    __py*=multiplier;
    return  *this;

}

// -------------------------------------------
//               TENSOR OPERATIONS
// -------------------------------------------

/**
 * @brief Calculates the Euclidean (L2) Norm of a tensor
 * 
 * @tparam T tensor type (int double etc.)
 * @return T returns a scalar of type T
 */
template<class T>
T tensor2<T>::norm(){
    return T(sqrt(pow(__px,2)+pow(__py,2)));
}

// ------------------------------------------
//              CONSTRUCTORS
// ------------------------------------------

/**
 * @brief Calculates the Manhattan (L1) Norm of a tensor
 * 
 * @tparam T tensor type (int double etc.)
 * @return T returns a scalar of type T
 */
template<class T>
T tensor2<T>::manhattan(){
    return __px+__py;
}


/**
 * @brief Construct a new tensor2<T>::tensor2 object
 * 
 * @tparam T 
 * @param other 
 */
template <class T>
tensor2<T>::tensor2(const tensor2<T>& other){
    __px = other.__px;
    __py = other.__py; 
}


/**
 * @brief Construct a new tensor2<T>::tensor2 object
 * 
 * @tparam T 
 * @param other 
 */
template <class T>
tensor2<T>::tensor2(tensor2<T>&& other){
    __px = std::move(other.__px);
    __py = std::move(other.__py); 
}
#endif