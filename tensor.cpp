// #include "tensor.h"




// // -------------------------------------------
// //        ARITHMETIC OPERATOR OVERLOAD
// // -------------------------------------------
// template <class T>
// tensor2<T>& tensor2<T>::operator=(const tensor2<T>& other){
//     __px = other.__px;
//     __py = other.__py; 
//     return *this;
// }

// template <class T>
// tensor2<T>& tensor2<T>::operator+=(const tensor2<T>& other){
//     __px += other.__px;
//     __py += other.__py; 
//     return *this;
// }

// template <class T>
// tensor2<T>& tensor2<T>::operator-=(const tensor2<T>& other){
//     __px -= other.__px;
//     __py -= other.__py; 
//     return *this;
// }

// template <class T>
// tensor2<T> tensor2<T>::operator+(const tensor2<T>& other){
//     tensor2<T> newpt(__px+other.__px,__py+other.__py);
//     return std::move(newpt);
// }

// template <class T>
// tensor2<T> tensor2<T>::operator-(const tensor2<T>& other){
//     tensor2<T> newpt(__px-other.__px,__py-other.__py);
//     return std::move(newpt);
// }
// template <class T>
// tensor2<T> tensor2<T>::operator/(const double& divisor){
//     tensor2<T> newpt(__px/divisor,__py/divisor);
//     return std::move(newpt);
// }

// template <class T>
// tensor2<T> tensor2<T>::operator*(const double& multiplier){
//     tensor2<T> newpt(__px*multiplier,__py*multiplier);
//     return std::move(newpt);
// }


// template <class T>
// tensor2<T>& tensor2<T>::operator/=(const double& divisor){
//     __px/=divisor;
//     __py/=divisor;
//     return  *this;
// }


// template <class T>
// tensor2<T>& tensor2<T>::operator*=(const double& multiplier){
//     __px*=multiplier;
//     __py*=multiplier;
//     return  *this;

// }

// // -------------------------------------------
// //               TENSOR OPERATIONS
// // -------------------------------------------

// /**
//  * @brief Calculates the Euclidean (L2) Norm of a tensor
//  * 
//  * @tparam T tensor type (int double etc.)
//  * @return T returns a scalar of type T
//  */
// template<class T>
// T tensor2<T>::norm(){
//     return T(sqrt(pow(__px,2)+pow(__py,2)));
// }

// // ------------------------------------------
// //              CONSTRUCTORS
// // ------------------------------------------

// /**
//  * @brief Calculates the Manhattan (L1) Norm of a tensor
//  * 
//  * @tparam T tensor type (int double etc.)
//  * @return T returns a scalar of type T
//  */
// template<class T>
// T tensor2<T>::manhattan(){
//     return __px+__py;
// }


// /**
//  * @brief Construct a new tensor2<T>::tensor2 object
//  * 
//  * @tparam T 
//  * @param other 
//  */
// template <class T>
// tensor2<T>::tensor2(const tensor2<T>& other){
//     __px = other.__px;
//     __py = other.__py; 
// }


// /**
//  * @brief Construct a new tensor2<T>::tensor2 object
//  * 
//  * @tparam T 
//  * @param other 
//  */
// template <class T>
// tensor2<T>::tensor2(tensor2<T>&& other){
//     __px = std::move(other.__px);
//     __py = std::move(other.__py); 
// }






// int main(){

//     tensor2<double> p1(10.0,20.0);
//     tensor2<double> p2 = p1*3.0;
//     tensor2<double> p3 = p2/100;
//     tensor2<double> p4(p3);
//     p4*=1000.0;

//     // tensor2<int> p2(p1);

//     // tensor2<int> p3(p1);


//     p3 = (p2+p1)/1000;
//     std::cout<<p2.norm()<<std::endl;
//     std::cout<<p1<<std::endl;
//     std::cout<<p2<<std::endl;
//     std::cout<<p3<<std::endl;
//     std::cout<<p4<<std::endl;

//     std::cout<<(p4-p2).norm()<<std::endl;
  






// }

