//
//  imagergb.cpp
//  colorDescriptor
//
//  Copyright 2022 Grigorios Piperagkas.  All rights reserved.
//  3-clause BSD License
//
//
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
//#include <stdlib.h>

class imagergb
{
private:
    int i,j,z;
    int sx,sy;
    int sizek;
    double *imgcx;
    int **imgrgb; //initial RGB channels
    double **sdq; //E, El, Ell : spectral differential quotients
    double **sdqg; //spatial differential quotients
    double **sdqg2; // second scale
    double **sdqdog;
    double **imgsc1; //first scale gaussian
    double **imgsc2; //second scale gaussian
    double **imgsc3; //third scale gaussian
    double **imgdog; //first d-o-g
    double **imgdog2; //second d-o-g
    double *k;
    double **Hmin; //hessian mins
    double **Hmax; //hessian maxs
    double sigma;

public:
    double **mins; //mins rgb
    double **maxs; //maxs rgb
    double **mines; //mins Ex Elx Ellx
    double **maxes; //maxs Ex Elx Ellx
    int sizex;
    int sizey;
    
    imagergb(std::string pathnew) //test img constructor
    {
        
        sigma=1.6;
        cv::Mat img;
        img = cv::imread(pathnew, cv::IMREAD_COLOR);
        
        sizex = img.rows;
        sizey = img.cols;
        sx=sizex;
        sy=sizey;
        
        sizek=3;
        
        
        imgcx = new double[sizex*sizey];
        
        
        sdq = new double*[3]; //E, El, Ell
        for (i=0;i<3;i++)
            sdq[i]= new double[sizex*sizey];
        
        sdqg = new double*[3]; //Ex, Elx, Ellx
        for (i=0;i<3;i++)
            sdqg[i]= new double[sizex*sizey];
        
        sdqg2 = new double*[3]; //Ex, Elx, Ellx
        for (i=0;i<3;i++)
            sdqg2[i]= new double[sizex*sizey];
        
        sdqdog = new double*[3];
        for (i=0;i<3;i++)
            sdqdog[i]= new double[sizex*sizey];
        
        imgrgb = new int*[3];
        for (i=0;i<3;i++)
            imgrgb[i] = new int[sizex*sizey];
        
        imgsc1 = new double*[3];
        for (i=0;i<3;i++)
            imgsc1[i]= new double[sizex*sizey];
        
        imgsc2 = new double*[3];
        for (i=0;i<3;i++)
            imgsc2[i]=new double[sizex*sizey];
        
        imgsc3 = new double*[3];
        for (i=0;i<3;i++)
            imgsc3[i]=new double[sizex*sizey];
        
        imgdog = new double*[3];
        for (i=0; i < 3; i++)
            imgdog[i] = new double[sizex*sizey];
        
        imgdog2 = new double*[3];
        for (i=0; i < 3; i++)
            imgdog2[i] = new double[sizex*sizey];
        
        mins = new double*[3];
        for (i=0;i<3;i++)
            mins[i] = new double[sizex*sizey];
        
        maxs = new double*[3];
        for (i=0;i<3;i++)
            maxs[i] = new double[sizex*sizey];
        
        mines = new double*[3];
        for (i=0;i<3;i++)
            mines[i] = new double[sizex*sizey];
        
        maxes = new double*[3];
        for (i=0;i<3;i++)
            maxes[i] = new double[sizex*sizey];

        
        Hmin = new double*[2];
        for (i=0;i<2;i++)
            Hmin[i] =new double[2];
        
        Hmax = new double*[2];
        for (i=0;i<2;i++)
            Hmax[i] =new double[2];
        
        
        k = new double[3*3];
        
        for (i=0;i<sizex;i++)
        {
            for (j=0;j<sizey;j++)
            {
                
                cv::Vec3b intensity = img.at<cv::Vec3b>(i,j);
                
                for (z=2;z>=0;z--)
                {
                    imgrgb[z][i*sizey+j] = int(intensity.val[z]);
                    
                }
                // fill the matrix of spectral differential quotients
                sdq[0][i*sizey+j]=0.06*imgrgb[0][i*sizey+j] + 0.63*imgrgb[1][i*sizey+j] + 0.27*imgrgb[2][i*sizey+j];
                sdq[1][i*sizey+j]=0.3*imgrgb[0][i*sizey+j] + 0.04*imgrgb[1][i*sizey+j] - 0.35*imgrgb[2][i*sizey+j];
                sdq[2][i*sizey+j]=0.34*imgrgb[0][i*sizey+j] -0.6*imgrgb[1][i*sizey+j] + 0.17*imgrgb[2][i*sizey+j];
            }
            
        }
        
    }

    
    imagergb(const std::string path[], int pp ) //base constructor
    {

        sigma=1.6;
        cv::Mat img;
        std::string pathh= path[pp];
        img = cv::imread(pathh, cv::IMREAD_COLOR);
        
        sizex = img.rows;
        sizey = img.cols;
        sx=sizex;
        sy=sizey;
        sizek=3;

        imgcx = new double[sizex*sizey];
        
        sdq = new double*[3]; //E, El, Ell
        for (i=0;i<3;i++)
            sdq[i]= new double[sizex*sizey];
        
        sdqg = new double*[3]; //Ex, Elx, Ellx
        for (i=0;i<3;i++)
            sdqg[i]= new double[sizex*sizey];
        
        sdqg2 = new double*[3]; //Ex, Elx, Ellx
        for (i=0;i<3;i++)
            sdqg2[i]= new double[sizex*sizey];
        
        sdqdog = new double*[3];
        for (i=0;i<3;i++)
            sdqdog[i]= new double[sizex*sizey];
        
        imgrgb = new int*[3];
        for (i=0;i<3;i++)
            imgrgb[i] = new int[sizex*sizey];
        
        imgsc1 = new double*[3];
        for (i=0;i<3;i++)
            imgsc1[i]= new double[sizex*sizey];
 
        imgsc2 = new double*[3];
        for (i=0;i<3;i++)
            imgsc2[i]=new double[sizex*sizey];
        
        imgsc3 = new double*[3];
        for (i=0;i<3;i++)
            imgsc3[i]=new double[sizex*sizey];
        
        imgdog = new double*[3];
        for (i=0; i < 3; i++)
            imgdog[i] = new double[sizex*sizey];

        imgdog2 = new double*[3];
        for (i=0; i < 3; i++)
            imgdog2[i] = new double[sizex*sizey];
        
        mins = new double*[3];
        for (i=0;i<3;i++)
            mins[i] = new double[sizex*sizey];

        maxs = new double*[3];
        for (i=0;i<3;i++)
            maxs[i] = new double[sizex*sizey];
        
        mines = new double*[3];
        for (i=0;i<3;i++)
            mines[i] = new double[sizex*sizey];

        maxes = new double*[3];
        for (i=0;i<3;i++)
            maxes[i] = new double[sizex*sizey];
        
        
        
        Hmin = new double*[2];
        for (i=0;i<2;i++)
            Hmin[i] =new double[2];
        
        Hmax = new double*[2];
        for (i=0;i<2;i++)
            Hmax[i] =new double[2];

        k = new double[3*3];

            for (i=0;i<sizex;i++)
            {
                for (j=0;j<sizey;j++)
                {

                    cv::Vec3b intensity = img.at<cv::Vec3b>(i,j);

                    for (z=2;z>=0;z--)
                    {
                        imgrgb[z][i*sizey+j] = int(intensity.val[z]);
                    }
                    // fill the matrix of spectral differential quotients
                    sdq[0][i*sizey+j]=0.06*imgrgb[0][i*sizey+j] + 0.63*imgrgb[1][i*sizey+j] + 0.27*imgrgb[2][i*sizey+j];
                    sdq[1][i*sizey+j]=0.3*imgrgb[0][i*sizey+j] + 0.04*imgrgb[1][i*sizey+j] - 0.35*imgrgb[2][i*sizey+j];
                    sdq[2][i*sizey+j]=0.34*imgrgb[0][i*sizey+j] - 0.6*imgrgb[1][i*sizey+j] + 0.17*imgrgb[2][i*sizey+j];

                }
    
            }

    }
 
    
    void convolutionG(int kk, int cl, int algflag) //produce scale images for each color
    {
        int qq, ww, ii, jj, mm, nn, mmm, nnn, iii, jjj;
        int kcenterx,kcentery;
        double gaus;
        const double pi=3.141592653589793238463;

        try {
            kcenterx=sizek/2;
            kcentery=sizek/2;
            

            //initialize
          for (qq=0;qq<sizex;qq++){
                for (ww=0;ww<sizey;ww++){
                    imgcx[qq*sizey+ww]=0;
                }
            }
      

            //convolution with gaussians with 3x3 kernel

            for (ii=1;ii<(sizex-1);ii++)
            {
                for (jj=1;jj<(sizey-1);jj++)
                {
                    for (mm=0;mm<sizek;mm++)
                    {
                        mmm=sizek-1-mm;
                        
                        for (nn=0;nn<sizek;nn++)
                        {
                            nnn=sizek-1-nn;

                            iii= ii + (mm - kcentery);
                            jjj= jj + (nn - kcenterx);
                            
                            gaus = (0.5/double(pi*kk*sigma))*exp(-double(mm*mm + nn*nn)/double(2*kk*kk*sigma*sigma));
                            
                            if (iii>=0 && iii<sizex && jjj>=0 && jjj<sizey)
                            {
                                if (algflag==0)
                                    imgcx[ii*sizey+jj] = imgcx[ii*sizey+jj]+ imgrgb[cl][iii*sizey+jjj]*gaus;
                                else
                                    imgcx[ii*sizey+jj] = imgcx[ii*sizey+jj] + sdq[cl][iii*sizey+jjj]*gaus;
                                
                            }
                            
                            
                        }
                    }
                }
            }

            
            
        } catch (const std::exception& t) {
            std::cerr << "exception convolutionG caught: " << t.what() << std::endl;
            
        }

        
    }
    
    void produceDoG() //produce difference of gaussians
    {
        int cc,ii,jj,kk;

        try {
            for (cc=0;cc<3;cc++)//colors
            {
                kk=1;
                this->convolutionG(kk,cc,0);
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        imgsc1[cc][ii*sizey+jj] = imgcx[ii*sizey+jj]; //first scale first
                    }
                }
                
            }
            for (cc=0;cc<3;cc++)//colors
            {
                kk=2;
                this->convolutionG(kk,cc,0);
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        imgsc2[cc][ii*sizey+jj] = imgcx[ii*sizey+jj]; // first scale second
                        
                    }
                }
                
            }
            for (cc=0;cc<3;cc++)//colors
            {
                kk=3;
                this->convolutionG(kk,cc,0);
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        imgsc3[cc][ii*sizey+jj] = imgcx[ii*sizey+jj]; //first scale third
                    }
                }
                
            }
            
            //subtract them to produce difference-of-gaussians.
            for (cc=0;cc<3;cc++)
            {
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        imgdog[cc][ii*sizey+jj] = std::abs(imgsc2[cc][ii*sizey+jj] - imgsc1[cc][ii*sizey+jj]);
                        imgdog2[cc][ii*sizey+jj] = std::abs(imgsc3[cc][ii*sizey+jj] - imgsc2[cc][ii*sizey+jj]);
                    }
                }
            }

            
            
        } catch (const std::exception& r) {
            std::cerr << "exception produceDoG caught: " << r.what() << std::endl;
            
        }
        
        
    }
    
    void produceSDQG() //produce Ex, Elx, Ellx and D-o-G
    {
        int cc,ii,jj,kk;
        try {
            for (cc=0;cc<3;cc++)
            {
                kk=1;
                this->convolutionG(kk,cc,1);
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        sdqg[cc][ii*sizey+jj] = imgcx[ii*sizey+jj];
                    }
                }
            }
            for (cc=0;cc<3;cc++)
            {
                kk=2;
                this->convolutionG(kk,cc,1);
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        sdqg2[cc][ii*sizey+jj] = imgcx[ii*sizey+jj];
                    }
                }
            }
                
            for (cc=0;cc<3;cc++)
            {
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        sdqdog[cc][ii*sizey+jj] = std::abs(sdqg2[cc][ii*sizey+jj] - sdqg[cc][ii*sizey+jj]);
                    }
                }
            }
   
                    

                } catch (const std::exception& b) {
                    std::cerr << "exception for SDQG caught:" << b.what() <<std::endl;
        }
        
        
    }
    
    void eliminateEdgeResp()
    {
        int sizeker=10;
        int ii,jj,c;
        double a,b;
        double r=10;
        double tr,det,check;
        
        check = (r+1)*(r+1)/r;
        
        for (c=0;c<3;c++)
        {
            for (i=5;i<(sx-5);i++)
            {
                for (j=5;j<(sy-5);j++)
                {
                    for (ii=0;ii<2;ii++)
                        for (jj=0;jj<2;jj++)
                        {
                            Hmin[ii][jj] =0;
                            Hmax[ii][jj] =0;
                        }
                    
                    if (mins[c][i*sy+j]<255){
                        for (ii=-5;ii<5;ii++)
                        {
                            if (mins[c][(i+ii)*sy+j]<255)
                                Hmin[0][0] = std::abs(mins[c][(i+ii)*sy+j] - mins[c][i*sy+j]);
                            
                        }
                        for (jj=-5;jj<5;jj++)
                        {
                            if (mins[c][i*sy+(j+jj)]<255)
                                Hmin[1][1] = std::abs(mins[c][i*sy+(j+jj)] - mins[c][i*sy+j]);
                            
                        }
                        for (ii=-5;ii<5;ii++)
                        {
                            for (jj=-5;jj<5;jj++)
                            {
                                if (ii!=0&&jj!=0)
                                {
                                    if (mins[c][(i+ii)*sy+(j+jj)]<255)
                                    {
                                        Hmin[0][1] = std::abs(mins[c][(i+ii)*sy+(j+jj)] - mins[c][i*sy+j]);
                                        Hmin[1][0] = Hmin[0][1];
                                    }
                                }
                            }
                        }
                        tr = Hmin[0][0] + Hmin[1][1];
                        det = Hmin[0][0]*Hmin[1][1] - Hmin[0][1]*Hmin[0][1];
                        if (det!=0)
                        {
                            if ((tr*tr/det)>check)
                                mins[c][i*sy+j]=255;
                        }
                        
                    }
                    
                    if (maxs[c][i*sy+j]>0){
                        for (ii=-5;ii<5;ii++)
                        {
                            if (maxs[c][(i+ii)*sy+j]>0)
                                Hmax[0][0] = std::abs(maxs[c][(i+ii)*sy+j] - maxs[c][i*sy+j]);
                            
                        }
                        for (jj=-5;jj<5;jj++)
                        {
                            if (maxs[c][i*sy+(j+jj)]>0)
                                Hmax[1][1] = std::abs(maxs[c][i*sy+(j+jj)] - maxs[c][i*sy+j]);
                            
                        }
                        for (ii=-5;ii<5;ii++)
                        {
                            for (jj=-5;jj<5;jj++)
                            {
                                if (ii!=0&&jj!=0)
                                {
                                    if (maxs[c][(i+ii)*sy+(j+jj)]>0)
                                    {
                                        Hmax[0][1] = std::abs(maxs[c][(i+ii)*sy+(j+jj)] - maxs[c][i*sy+j]);
                                        Hmax[1][0] = Hmax[0][1];
                                    }
                                }
                            }
                        }
                        tr = Hmax[0][0] + Hmax[1][1];
                        det = Hmax[0][0]*Hmax[1][1] - Hmax[0][1]*Hmax[0][1];
                        if (det!=0)
                        {
                            if ((tr*tr/det)>check)
                                maxs[c][i*sy+j]=0;
                        }
                        
                    }
                    
                }
            }
            
        }
    }

    
            
    void extremaLocalization()
    {
        int ii,jj,cc,iii,jjj,l,m,ll,mm,lll,mmm,skx,sky;
        double **cnvs;
        double **cnvs2;
        double **cnvssdq;
        double min;
        double max;
        double mine;
        double maxe;
        
        
        
        skx=10;
        sky=10;
        cnvs = new double*[skx];
        for (i=0;i<skx;i++)
            cnvs[i] = new double[sky];
        cnvs2 = new double*[skx];
        for (i=0;i<skx;i++)
            cnvs2[i] = new double[sky];
        
        cnvssdq = new double*[skx];
        for (i=0;i<skx;i++)
            cnvssdq[i]= new double[sky];
        
        
        try {
            //initialization
            for (cc=0;cc<3;cc++)
            {
                for (ii=0;ii<sizex;ii++)
                {
                    for (jj=0;jj<sizey;jj++)
                    {
                        mins[cc][ii*sizey+jj]= 255;
                        maxs[cc][ii*sizey+jj]= 0;
                    }
                }
            }
            
            
            for (cc=0;cc<3;cc++)
            {
                for (iii=5;iii<(sizex-5);iii++)
                {
                    for (jjj=5;jjj<(sizey-5);jjj++)
                    {
                        for (l=-5;l<5;l++)
                        {
                            for (m=-5;m<5;m++)
                            {
                                cnvs[l+5][m+5]=imgdog[cc][(iii+l)*sizey+(jjj+m)];
                                cnvs2[l+5][m+5]=imgdog2[cc][(iii+l)*sizey+(jjj+m)];
                                cnvssdq[l+5][m+5]=sdqdog[cc][(iii+l)*sizey+(jjj+m)];
                                
                            }
                        }
                        min=255;
                        max=0;
                        
                        mine=255;
                        maxe=-255;
                        
                        for (l=0;l<skx;l++) // find mins and maxs of skx x sky kernel
                        {
                            for (m=0;m<sky;m++)
                            {
                                if (cnvs[l][m]<min)
                                {
                                    min=cnvs[l][m];
                                    ll=l + iii -5;
                                    mm=m + jjj -5;
                                }
                                if (cnvs2[l][m]<min)
                                {
                                    min=cnvs2[l][m];
                                    ll=l + iii -5;
                                    mm=m + jjj -5;
                                }
                                if (cnvssdq[l][m]<mine)
                                {
                                    mine=cnvssdq[l][m];
                                    ll=l + iii -5;
                                    mm=m + jjj -5;
                                }
                                if (cnvs[l][m]>max)
                                {
                                    max=cnvs[l][m];
                                    lll=l + iii -5;
                                    mmm=m + jjj -5;
                                }
                                if (cnvs2[l][m]>max)
                                {
                                    max=cnvs2[l][m];
                                    lll=l + iii -5;
                                    mmm=m + jjj -5;
                                }
                                if (cnvssdq[l][m]>maxe)
                                {
                                    maxe=cnvssdq[l][m];
                                    ll=l + iii -5;
                                    mm=m + jjj -5;
                                }

                            }
                        }
                        
                        mins[cc][ll*sizey+mm]=min;
                        maxs[cc][lll*sizey+mmm]=max;
                        mines[cc][ll*sizey+mm]=mine;
                        maxes[cc][ll*sizey+mm]=maxe;
                        
                        
                    }
                    
                }
                
            }

        
        
        } catch (const std::exception& w) {
            std::cerr << "exception extremaLocalization caught: " << w.what() << std::endl;
            
        }
        

        
    }
    
    
    
    void printall() //print all stages of image to pngs
    {
       
       
        cv::Mat out1 = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat out2 = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat out3 = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat dog1 = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat dog2 = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat minss = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat maxss = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat miness = cv::Mat::zeros(sizex,sizey, CV_8UC3);
        cv::Mat maxess = cv::Mat::zeros(sizex,sizey, CV_8UC3);

        for (i=0;i<sizex;i++)
        {
            for (j=0;j<sizey;j++)
            {
                out1.at<cv::Vec3b>(i,j) = cv::Vec3b(int(imgsc1[0][i*sizey+j]),int(imgsc1[1][i*sizey+j]),int(imgsc1[2][i*sizey+j]));
                out2.at<cv::Vec3b>(i,j) = cv::Vec3b(int(imgsc2[0][i*sizey+j]),int(imgsc2[1][i*sizey+j]),int(imgsc2[2][i*sizey+j]));
                out3.at<cv::Vec3b>(i,j) = cv::Vec3b(int(imgsc3[0][i*sizey+j]),int(imgsc3[1][i*sizey+j]),int(imgsc3[2][i*sizey+j]));
                dog1.at<cv::Vec3b>(i,j) = cv::Vec3b(int(imgdog[0][i*sizey+j]),int(imgdog[1][i*sizey+j]),int(imgdog[2][i*sizey+j]));
                dog2.at<cv::Vec3b>(i,j) = cv::Vec3b(int(imgdog2[0][i*sizey+j]),int(imgdog2[1][i*sizey+j]),int(imgdog2[2][i*sizey+j]));
                minss.at<cv::Vec3b>(i,j) = cv::Vec3b(int(mins[0][i*sizey+j]),int(mins[1][i*sizey+j]),int(mins[2][i*sizey+j]));
                maxss.at<cv::Vec3b>(i,j) = cv::Vec3b(int(maxs[0][i*sizey+j]),int(maxs[1][i*sizey+j]),int(maxs[2][i*sizey+j]));
                miness.at<cv::Vec3b>(i,j) = cv::Vec3b(int(mines[0][i*sizey+j]),int(mines[1][i*sizey+j]),int(mines[2][i*sizey+j]));
                maxess.at<cv::Vec3b>(i,j) = cv::Vec3b(int(maxes[0][i*sizey+j]),int(maxes[1][i*sizey+j]),int(maxes[2][i*sizey+j]));

            }
        }

        cv::imwrite("imgsc1.png",out1);
        cv::imwrite("imgsc2.png",out2);
        cv::imwrite("imgsc3.png",out3);
        cv::imwrite("imgdog.png",dog1);
        cv::imwrite("imgdog2.png",dog2);
        cv::imwrite("mins.png",minss);
        cv::imwrite("maxs.png",maxss);
        cv::imwrite("minsE.png", miness);
        cv::imwrite("maxsE.png", maxess);



    }
    
    /*
    */
    ~imagergb() //destructor
    {
        
        delete[] k;
        delete[] imgcx;
        
        for (i=0;i<3;i++)
            delete imgrgb[i];
        delete[] imgrgb;
        for (i=0;i<3;i++)
            delete imgsc1[i];
        delete[] imgsc1;
        for (i=0;i<3;i++)
            delete imgsc2[i];
        delete[] imgsc2;
        for (i=0;i<3;i++)
            delete imgsc3[i];
        delete[] imgsc3;
        for (i=0;i<3;i++)
            delete imgdog[i];
        delete[] imgdog;
        for (i=0;i<3;i++)
            delete imgdog2[i];
        delete[] imgdog2;
        for (i=0;i<3;i++)
            delete mins[i];
        delete[] mins;
        for (i=0;i<3;i++)
            delete maxs[i];
        delete[] maxs;
        for (i=0;i<3;i++)
            delete sdq[i];
        delete[] sdq;
        for (i=0;i<3;i++)
            delete sdqg[i];
        delete[] sdqg;
        for (i=0;i<3;i++)
            delete sdqg2[i];
        delete[] sdqg2;
        for (i=0;i<3;i++)
            delete mines[i];
        delete[] mines;
        for (i=0;i<3;i++)
            delete maxes[i];
        delete[] maxes;
        
        for (i=0;i<2;i++)
            delete Hmin[i];
        for (i=0;i<2;i++)
            delete Hmax[i];
        
        delete[] Hmin;
        delete[] Hmax;


        
    }
    
};
