//
//  main.cpp
//  colorDescriptor
//
//  Copyright 2022 Grigorios Piperagkas. All rights reserved.
//  3-clause BSD License
//

#include <stdio.h>
#include<iostream>
#include <string>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include <math.h>
#include <vector>
#include <stdlib.h>  
#include <exception>      

#include "classifier.cpp"


using namespace std;
using namespace cv;

const int nn=10;
//std::string path[nn]= {"1.jpg", "2.jpg", "3.jpg", "4.jpg", "5.jpg", "6.jpg", "7.jpg", "8.jpg", "9.jpg", "10.jpg"};
std::string path[nn]= {"1.png", "2.png", "3.png", "4.png", "5.png", "6.png", "7.png", "8.png", "9.png", "10.png"};


int main(int argc, char **argv)
{
    
    classifier cl(path, nn); //initialize classifier
    
    
    cl.calcfVec(); //build vector of base images and feature vectors
 
    
    cl.calcDistances(); //bring new image "test.jpg", build feature vector and calculate distances
 
    
    cl.printdists(); //print distances for RGB channels and Ex Elx Ellx for each base image to test.jpg

    return 0;
}



