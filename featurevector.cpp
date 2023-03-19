//
//  featurevector.cpp
//  colorDescriptor
//
//  Copyright 2022 Grigorios Piperagkas.  All rights reserved.
//  3-clause BSD License
//
//
//#include <stdlib.h>

class featureVector{
    
private:
    int **x;
    int **y;
    int **xe;
    int **ye;
    double **mins;
    double **maxs;
    int i,ii,j,jj,c,sx,sy;
    
public:
    
    int **fmins;
    int **fmaxs;
    int fsize;
    int **fmines;
    int **fmaxes;
    int fsizee;
    
    featureVector(double **mins, double **maxs, double **mines, double **maxes, int sizex, int sizey) //constructor
    {
        sx=sizex;
        sy=sizey;
        c=0;
        fsize=0;
        ii=0;
        for (i=0;i<sx;i++)
        {
            for (j=0;j<sy;j++)
            {
                if (mins[c][i*sy+j]<255 || maxs[c][i*sy+j]>0)
                    fsize++;
                if (mines[c][i*sy+j]<255 || maxes[c][i*sy+j]>-255)
                    fsizee++;
            }
        }
 
        fmins = new int*[3];
        for (i=0;i<3;i++)
            fmins[i] = new int[fsize];
        fmaxs = new int*[3];
        for (i=0;i<3;i++)
            fmaxs[i]=new int[fsize];
        
        fmines = new int*[3];
        for (i=0;i<3;i++)
            fmines[i] = new int[fsizee];
        fmaxes = new int*[3];
        for (i=0;i<3;i++)
            fmaxes[i]=new int[fsizee];

        
        x = new int*[3];
        for (i=0;i<3;i++)
            x[i]=new int[fsize];
        y = new int*[3];
        for (i=0;i<3;i++)
            y[i]=new int[fsize];
        
        xe = new int*[3];
        for (i=0;i<3;i++)
            xe[i]=new int[fsizee];
        ye = new int*[3];
        for (i=0;i<3;i++)
            ye[i]=new int[fsizee];

    }
 
    
    void buildfVector(double **mins, double **maxs, double **mines, double **maxes) //build extrema vectors of size fsize with coordinates x y for each color RGB and Ex, Elx, Ellx
    {
        //try {
            
            for (c=0;c<3;c++)
            {
                ii=0;
                jj=0;
                for (i=0;i<sx;i++)
                {
                    for (j=0;j<sy;j++)
                    {
                        if ((mins[c][i*sy+j]<255 || maxs[c][i*sy+j]>0)&&(ii<fsize))
                        {
                            x[c][ii]=i;
                            y[c][ii]=j;
                            fmins[c][ii]=int(mins[c][i*sy+j]);
                            fmaxs[c][ii]=int(maxs[c][i*sy+j]);
                            ii++;
                        }
                        if ((mines[c][i*sy+j]<255 || maxes[c][i*sy+j]>-255)&&(jj<fsizee))
                        {
                            xe[c][ii]=i;
                            ye[c][ii]=j;
                            fmines[c][ii]=int(mines[c][i*sy+j]);
                            fmaxes[c][ii]=int(maxes[c][i*sy+j]);
                            jj++;
                            
                        }

                    }
                }
            }

        
        
        //} catch (const std::exception& e) {
        //    std::cerr << "exception buildfVector caught: " << e.what() << std::endl;

        //}
        
    
        
    }
    

    
 
    

    ~featureVector() //destructor
    {
        
        for (i=0;i<3;i++){
            delete fmins[i];
            delete fmaxs[i];
            delete x[i];
            delete y[i];
            delete fmines[i];
            delete fmaxes[i];
            delete xe[i];
            delete ye[i];
            
        }

        delete[] fmines;
        delete[] fmaxes;
        delete[] xe;
        delete[] ye;
        delete[] fmins;
        delete[] fmaxs;
        delete[] x;
        delete[] y;

    }
};
