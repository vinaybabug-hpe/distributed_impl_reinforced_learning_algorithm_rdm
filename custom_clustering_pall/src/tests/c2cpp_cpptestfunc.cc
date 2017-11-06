/*
 * callc2cpp_wrapper.cc
 *
 *  Created on: Aug 27, 2015
 *      Author: vinaya
 */




#include <iostream>


extern "C" void printHello_cpp();



using namespace std;
void printHello_cpp()
{


    std::cout << "Hello World!: FROM CPP.." << std::endl;


    return;
}
