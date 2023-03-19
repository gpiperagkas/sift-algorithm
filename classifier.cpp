//
//  classifier.cpp
//  colorDescriptor
//
//  Copyright 2022 Grigorios Piperagkas.  All rights reserved.
//  3-clause BSD License
//
//
#include<iostream>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <exception>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>

#include "imagergb.cpp"
#include "featurevector.cpp"



class classifier{
    
private:
    int s,l,m,ll;
    int sizen;
    double **dists;
    double **distse;
    std::vector<imagergb*> base;
    std::vector<featureVector*> fV;
    std::string pathnew;
    
public:

    classifier(const std::string path[], int nn )
    {
        sizen=nn;
        fV.reserve(nn);
        base.reserve(nn);
        for (s=0;s<nn;s++)
        {

            imagergb* img= new imagergb(path, s );
            base.push_back(img);

        }
        
        pathnew= "test.png";  //test image path

        
        dists = new double*[3];
        for (m=0;m<3;m++)
            dists[m]= new double[nn];
        
        distse = new double*[3];
        for (m=0;m<3;m++)
            distse[m]= new double[nn];
        
        
    }
  
    
    
    void calcfVec(){                  //build a vector of feature vectors
        try {

            for(unsigned int ll=0; ll<base.size(); ll++)
            {

                base.at(ll)->produceDoG();
                
                base.at(ll)->produceSDQG();
                
                base.at(ll)->eliminateEdgeResp();

                base.at(ll)->extremaLocalization();

                featureVector* fVec= new featureVector(base.at(ll)->mins, base.at(ll)->maxs, base.at(ll)->mines, base.at(ll)->maxes, base.at(ll)->sizex, base.at(ll)->sizey);

                fVec->buildfVector(base.at(ll)->mins, base.at(ll)->maxs, base.at(ll)->mines, base.at(ll)->maxes);

                fV.push_back(fVec);
                
                
            }
        } catch (const std::exception& e) {
                std::cerr << "exception calcfVec caught: " << e.what() << std::endl;
        }
    }
  
    
    
    void calcDistances(){  //produce feature vector of new image and calculate euclidean distances
       
        int c,ci,co;
        int ssize;
        double ddmin, ddmax, dmin, dmax, ff;
        double tmp;
        double tmp1, tmp2;
        
        try {
            imagergb* newimg = new imagergb(pathnew);
        
            newimg->produceDoG();
            
            newimg->produceSDQG();
            
            newimg->eliminateEdgeResp();
        
            newimg->extremaLocalization();
            
            featureVector* newfV = new featureVector(newimg->mins, newimg->maxs, newimg->mines, newimg->maxes, newimg->sizex, newimg->sizey);
        
            newfV->buildfVector(newimg->mins, newimg->maxs, newimg->mines, newimg->maxes);
        
            newimg->printall();
    
            co=0;
            for(unsigned int ll=0; ll<fV.size(); ll++) //distances for rgb sift
            {
                ssize = fV.at(ll)->fsize;
                if (newfV->fsize < ssize)
                    ssize = newfV->fsize;
                for (c=0;c<3;c++)
                {
                    dists[c][co]=0;
                    
                    for (ci=0;ci<ssize;ci++)
                    {
                        ddmin = fV.at(ll)->fmins[c][ci];
                        ddmax = fV.at(ll)->fmaxs[c][ci];
                        dmin = newfV->fmins[c][ci];
                        dmax = newfV->fmaxs[c][ci];

                        if ((dmax+ddmax)==0)
                        {
                            tmp2=0;
                            if ((dmin+ddmin)!=0)
                            {
                                tmp1= std::abs(ddmin-dmin);
                                tmp= (0.5)*((tmp1*tmp1)/(ddmin+dmin));
                            }else{
                                tmp=0;
                            }
                        }
                        else if ((dmin+ddmin)==0)
                        {
                            tmp1=0;
                            if ((ddmax+dmax)!=0){
                                tmp2 = std::abs(ddmax-dmax);
                                tmp = (0.5)*((tmp2*tmp2)/(ddmax+dmax));
                            }else{
                                tmp=0;
                            }
                        }
                        else{
                            tmp1= std::abs(ddmin-dmin);
                            tmp2= std::abs(ddmax-dmax);
                            tmp = (0.5)*((tmp1*tmp1)/(ddmin+dmin))+(0.5)*((tmp2*tmp2)/(ddmax+dmax));
                        }

                        dists[c][co] = dists[c][co]+ tmp;
                    }

                }
             co++;
            }
            
            tmp=0;
            tmp1=0;
            tmp2=0;
            co=0;
            for(unsigned int ll=0; ll<fV.size(); ll++) //distances for csift
            {
                ssize = fV.at(ll)->fsizee;
                if (newfV->fsizee < ssize)
                    ssize = newfV->fsizee;
                for (c=0;c<3;c++)
                {
                    distse[c][co]=0;
                    
                    for (ci=0;ci<ssize;ci++)
                    {
                        ddmin = fV.at(ll)->fmines[c][ci];
                        ddmax = fV.at(ll)->fmaxes[c][ci];
                        dmin = newfV->fmines[c][ci];
                        dmax = newfV->fmaxes[c][ci];
                        
                        if ((dmax+ddmax)==0)
                        {
                            tmp2=0;
                            if ((dmin+ddmin)!=0)
                            {
                                tmp1= std::abs(ddmin-dmin);
                                tmp= (0.5)*((tmp1*tmp1)/(ddmin+dmin));
                            }else{
                                tmp=0;
                            }
                        }
                        else if ((dmin+ddmin)==0)
                        {
                            tmp1=0;
                            if ((ddmax+dmax)!=0){
                                tmp2 = std::abs(ddmax-dmax);
                                tmp = (0.5)*((tmp2*tmp2)/(ddmax+dmax));
                            }else{
                                tmp=0;
                            }
                        }
                        else{
                            tmp1= std::abs(ddmin-dmin);
                            tmp2= std::abs(ddmax-dmax);
                            tmp = (0.5)*((tmp1*tmp1)/(ddmin+dmin))+(0.5)*((tmp2*tmp2)/(ddmax+dmax));
                        }
                        distse[c][co] = distse[c][co]+ tmp;
                    }
                }
                co++;
            }

            
        } catch (const std::exception& e) {
            std::cerr << "exception calcDistances caught: " << e.what() << std::endl;
        }

        

    
    }
    
    

    
    void printdists() // print distances of base images to new image
    {
        int c,co;
  
        std::cout << "distances per color and image (rgb SIFT, 10 base images): " << std::endl;
        std::cout << " "<< std::endl;

        for (c=0;c<3;c++)
        {
            for (co=0;co<sizen;co++)
            {
                std::cout <<" " << dists[c][co];
            }
            std::cout <<" " <<std::endl;
            
        }
        
        std::cout << " "<< std::endl;
        std::cout << "distances per Ex, Elx, Ellx and image (CSIFT - 10 base images):" <<std::endl;
        std::cout << " "<< std::endl;
        
        for (c=0;c<3;c++)
        {
            for (co=0;co<sizen;co++)
            {
                std::cout <<" " << distse[c][co];
            }
            std::cout <<" " <<std::endl;
            
        }

    }
    
    ~classifier()
    {
        fV.clear();
        base.clear();
        
        for (l=0;l<3;l++)
            delete dists[l];
        delete[] dists;
        for (l=0;l<3;l++)
            delete distse[l];
        delete[] distse;
        
    }
    
    
};

