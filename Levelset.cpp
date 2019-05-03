#include <iostream>
#include <math.h>
#include <iomanip>
#include <array>
#include <algorithm> 
#include <fstream>

struct Grad{
    double** gradx;
    double** grady;
};


double** make_img(int nx, int ny)
{
    double** img = new double *[ny];
    for (int i=0;i<ny;i++) 
    {
        img[i] = new double[nx];
    }
    
    
    return img;
}


double **gaussian_mask(int shape=31, double sigma=6)
{
    double** res=make_img(shape,shape);
    double sum=0;
    int mid = (shape-1)/2;

    for (int i=0;i<shape;i++)
    {
        for (int j=0;j<shape;j++)
        {
            res[i][j] = exp(-((i-mid)*(i-mid)+(j-mid)*(j-mid))/(2*sigma*sigma));
            sum += res[i][j];
        }
    }

    for (int i=0;i<shape;i++)
    {
        for (int j=0;j<shape;j++)
        {
            res[i][j] = res[i][j]/sum;
        }
    }
    
    return res;
}

Grad gradient(double **img, int nx, int ny)
{
    
    double** gradx = make_img(nx,ny);
    double** grady = make_img(nx,ny);

    for (int i=1;i<ny-1;i++)
    {
        for (int j=1;j<nx-1;j++)
        {
            gradx[i][j]=(img[i][j+1]-img[i][j-1])/2;
            grady[i][j]=(img[i+1][j]-img[i-1][j])/2;   
        }
    }
    for (int i=0;i<ny;i++)
    {
        gradx[i][0]=0;
        gradx[i][nx-1]=0;
    }
    for (int i=0;i<nx;i++)
    {
        grady[0][i]=0;
        grady[ny-1][i]=0;
    }
    
    Grad gradimg;
    gradimg.gradx=gradx;
    gradimg.grady=grady;
    
    return gradimg;
}

double Part(double** img,double** u,int nx,int ny)
{
    double sum1=0;
    double sum2=0;
    int count1=0;
    int count2=0;
    double res;
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
        {
            if(u[i][j]>0)
            {
                sum1+=img[i][j];
                count1++;
            }
            else 
                sum2+=img[i][j];
            
        } 
    }
    
    count2=nx*ny-count1;
    res=(sum1/count1+sum2/count2)/2;
    
    return res;
}

void conv(double **img,int nx,int ny,int shape=5)
{
    int mid=(shape-1)/2;
    double** gaussian;
    gaussian=gaussian_mask(shape);
    
    double** u = make_img(nx+2*mid,ny+2*mid);
    
    for (int i=0;i<ny+2*mid;i++) 
        for (int j=0;j<nx+2*mid;j++)
            u[i][j]=0;
    
    for (int i;i<ny;i++) 
    {
        for (int j;j<nx;j++)
            u[i+mid][j+mid]=img[i][j];
    }
    
    for (int i;i<ny;i++) 
    {
        for (int j;j<nx;j++)
        {
            int sum=0;
            for(int m=0;m<shape;m++)
            {
                for(int n=0;n<shape;n++)
                {
                    sum+=gaussian[m][n]*u[i+m][j+n];
                }
            }
            u[i][j]=sum;
        }
           
    }
    
}

double** countour(double **img,int nx,int ny,int Iter=40)
{
    double** phi = make_img(nx,ny);
    double** u = make_img(nx,ny);
    double** spf= make_img(nx,ny);
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
        {
            phi[i][j]=-1;
            u[i][j]=1;
        }
    }
    for (int i=0;i<ny;i++)
    {
        u[i][0]=-1;
        u[i][nx-1]=-1;
        phi[i][0]=1;
        phi[i][nx-1]=1;
    }
    for (int i=0;i<nx;i++)
    {
        u[0][i]=-1;
        u[ny-1][i]=-1;
        phi[0][i]=1;
        phi[ny-1][i]=1;
    }
 
    
    
    double delt=0.4;
    double mu=15;
    double valueinter=delt*mu;
    
    double c;
    
    for (int i=0;i<Iter;i++)
    {

        Grad gradimg;
        gradimg= gradient(img,nx,ny);
        double M=-1;
        c=Part(img,u,nx,ny);
        
        for (int i=0;i<ny;i++) 
        {   
            for (int j=0;j<nx;j++)
            {
                spf[i][j]=img[i][j]-c;
                M=std::max(std::abs(spf[i][j]),M);
            }
        }
        
        for (int i=0;i<ny;i++) 
        {   
            for (int j=0;j<nx;j++)
            {
                spf[i][j]=spf[i][j]/M;    
            }
        }
       
        for (int i=0;i<ny;i++) 
        {   
            for (int j=0;j<nx;j++)
            {
                u[i][j]+=valueinter*spf[i][j]*std::sqrt(gradimg.gradx[i][j]*gradimg.gradx[i][j]+gradimg.grady[i][j]*gradimg.grady[i][j]); 
            }
        }

        
        for (int i=0;i<ny;i++) 
        {
            for (int j=0;j<nx;j++)
            {
                if (u[i][j]>0) u[i][j]=1;
                if (u[i][j]<0) u[i][j]=-1;
            } 
        }
        
        conv(u,nx,ny,5);

        //convolution
    }
    
    return u;
}

void save(double **img,int nx,int ny)
{
    std::ofstream f ("toto.pgm");
    if (!f.is_open())
        std::cout << "Impossible d'ouvrir le fichier en écriture !" << std::endl;
    else
    {
        f<<"P2"<<std::endl;
        f<<nx<<' '<<ny<<std::endl;
        f<<3<<std::endl;
        for (int i=0;i<ny;i++) 
        {
             for (int j=0;j<nx;j++)
            {
                f<<img[i][j]+1<<' '; 
            } 
            f<<std::endl;
        }
    }
    f.close();
}

    
int main()
{
    int shape=5;
    double** gaussian;
    gaussian=gaussian_mask(shape);
    double** u;
    int nx=40;
    int ny=40;
    
    double** img =make_img(nx,ny);
    
    for (int i=0;i<ny;i++) 
    {
         for (int j=0;j<nx;j++)
        {
            if (i==20)
                img[i][j]=100;
            if (j==15)
                img[i][j]=100;
            if ((i-20)*(i-20)+(j-25)*(j-25)<50)
                img[i][j]=200;
             
        } 
    }
    
    u=countour(img,nx,ny,10);
    save(u,nx,ny);
    
    
    return 0;
}

