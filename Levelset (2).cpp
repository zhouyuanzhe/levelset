#include <iostream>
#include <math.h>
#include <array>
#include <algorithm> 
#include <fstream>
#include <vector>


struct Grad
{
    double** gradx;
    double** grady;
};

struct Image
{
    double **image_;
    int nx;
    int ny;
    
};


double** Make_img(int nx, int ny)
{
    double** img = new double *[ny];
    for (int i=0;i<ny;i++) 
    {
        img[i] = new double[nx];
    }
    
    
    return img;
}


double** Gaussian_mask(int shapeGaus=5, double sigma=6)
{
    double** res=Make_img(shapeGaus,shapeGaus);
    double sum = 0;
    int mid = (shapeGaus-1)/2;

    double sig=2*sigma*sigma;
    for (int i=0;i<shapeGaus;i++)
    {
        for (int j=0;j<shapeGaus;j++)
        {
            res[i][j] = exp(-((i-mid)*(i-mid)+(j-mid)*(j-mid))/sig);
            sum += res[i][j];
        }
    }

    for (int i=0;i<shapeGaus;i++)
    {
        for (int j=0;j<shapeGaus;j++)
        {
            res[i][j] = res[i][j]/sum;
        }
    }
    
    return res;
}

Grad Gradient(Image image)
{
    int nx = image.nx;
    int ny = image.ny;
    double **img = image.image_;
    
    
    double** gradx = Make_img(nx,ny);
    double** grady = Make_img(nx,ny);

    for (int i=1;i<ny-1;i++)
    {
        for (int j=1;j<nx-1;j++)
        {
            gradx[i][j] = (img[i][j+1]-img[i][j-1])/2;
            grady[i][j] = (img[i+1][j]-img[i-1][j])/2;   
        }
    }
    for (int i=0;i<ny;i++)
    {
        gradx[i][0] = 0;
        gradx[i][nx-1] = 0;
    }
    for (int i=0;i<nx;i++)
    {
        grady[0][i] = 0;
        grady[ny-1][i] = 0;
    }
    
    Grad gradimg;
    gradimg.gradx=gradx;
    gradimg.grady=grady;
    
    return gradimg;
}

double Partition(Image image, double** u)
{
    int nx = image.nx;
    int ny = image.ny;
    double** img = image.image_;
    
    double res;
    double sum1 = 0;
    double sum2 = 0;
    int count1 = 0;
    int count2 = 0;
    
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
        {
            if(u[i][j]>0)
            {
                sum1 += img[i][j];
                count1++;
            }
            else 
                sum2 += img[i][j];
            
        } 
    }
    
    count2 = nx*ny-count1;
    res = (sum1/count1+sum2/count2)/2;
    
    return res;
}

double Medianval(std::vector<double> v)
{
    std::sort (v.begin(),v.end()); 
    
    return v[v.size()/2+1];
}
//Median valeur d'une liste de data

void Median(double** img, int nx, int ny , int shapeMedian=5)
{
    int mid = (shapeMedian-1)/2;
    
    std::vector<double> v;
    
    int Ny = ny+2*mid;
    int Nx = nx+2*mid;
    
    double** u = Make_img(Nx, Ny);
    
    for (int i=0;i<Ny;i++) 
        for (int j=0;j<Nx;j++)
            u[i][j] = 0;
    
    
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
            u[i+mid][j+mid] = img[i][j];
    }
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
        {
            v.resize(0);
            for(int m=0;m<shapeMedian;m++)
            {
                for(int n=0;n<shapeMedian;n++)
                {
                    v.push_back(u[i+m][j+n]);
                }
            }
            u[i][j] = Medianval(v);
        }
           
    }
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
            img[i][j] = u[i+mid][j+mid];
    }
    
}

void Conv(double** img, int nx, int ny, int shapeGaus=5)
{    
    int mid = (shapeGaus-1)/2;
    double** gaussian;
    gaussian = Gaussian_mask(shapeGaus);
    
    int sum;
    int Ny = ny+2*mid;
    int Nx = nx+2*mid;
    
    double** u = Make_img(Nx,Ny);
    
    for (int i=0;i<mid;i++)
    {
        for (int j=0;j<Nx;j++)
            u[i][j] = 0;
    }
    for (int i=ny+mid;i<Ny;i++)
    {
        for (int j=0;j<Ny;j++)
            u[i][j] = 0;
    }
    for (int i=mid;i<ny+mid;i++)
    {
        for (int j=0;j<mid;j++)
            u[i][j] = 0;
        for (int j=nx+mid;j<Nx;j++)
            u[i][j] = 0;
    }

    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
            u[i+mid][j+mid] = img[i][j];
    }
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
        {
            sum = 0;
            for(int m=0;m<shapeGaus;m++)
            {
                for(int n=0;n<shapeGaus;n++)
                {
                    sum += gaussian[m][n]*u[i+m][j+n];
                }
            }
            u[i][j] = sum;
        }
           
    }
        
    for (int i;i<ny;i++) 
    {
        for (int j;j<nx;j++)
            img[i][j] = u[i+mid][j+mid];
    }
    
}

double** Countour(Image image,int Iter=40,int gaus=5, double delt=0.4, double mu=15)
{
    int nx=image.nx;
    int ny=image.ny;
    double** img=image.image_;
    
    double** phi = Make_img(nx,ny);
    double** u = Make_img(nx,ny);
    double** spf = Make_img(nx,ny);
    
    for (int i=0;i<ny;i++) 
    {
        for (int j=0;j<nx;j++)
        {
            phi[i][j] = -1.;
            u[i][j] = 1.;
        }
    }
    for (int i=0;i<ny;i++)
    {
        u[i][0] = -1.;
        u[i][nx-1] = -1.;
        phi[i][0] = 1.;
        phi[i][nx-1] = 1.;
    }
    for (int i=0;i<nx;i++)
    {
        u[0][i] = -1.;
        u[ny-1][i] = -1.;
        phi[0][i] = 1.;
        phi[ny-1][i] = 1.;
    }
    //set default values
    
    double valueinter=delt*mu;
        
    double c;
    
    for (int i=0;i<Iter;i++)
    {

        Grad gradimg;
        gradimg= Gradient(image);
        double M = -1;
        c = Partition(image,u);
        
        for (int i=0;i<ny;i++) 
        {   
            for (int j=0;j<nx;j++)
            {
                spf[i][j] = img[i][j]-c;
                M = std::max(std::abs(spf[i][j]),M);
            }
        }
        //can be optimised to find the max value?
        
        for (int i=0;i<ny;i++) 
        {   
            for (int j=0;j<nx;j++)
            {
                spf[i][j] = spf[i][j]/M;    
            }
        }
       
        for (int i=0;i<ny;i++) 
        {   
            for (int j=0;j<nx;j++)
            {
                u[i][j] += valueinter*spf[i][j]*std::sqrt(gradimg.gradx[i][j]*gradimg.gradx[i][j]+gradimg.grady[i][j]*gradimg.grady[i][j]); 
            }
        }

        
        for (int i=0;i<ny;i++) 
        {
            for (int j=0;j<nx;j++)
            {
                if (u[i][j]>0) u[i][j] = 1.;
                if (u[i][j]<0) u[i][j] = -1.;
            } 
        }
        
        Conv(u,nx,ny,gaus);
    //convolution
    }
    
    return u;
}

void Savepgm(double** img,int nx,int ny, int maxval=3)
{ 
    std::ofstream f ("toto.pgm");
    if (!f.is_open())
        std::cout << "Impossible d'ouvrir le fichier en écriture !" << std::endl;
    else
    {
        f<<"P2"<<std::endl;
        f<<nx<<' '<<ny<<std::endl;
        f<<maxval<<std::endl;
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

Image Loadpgm(const char* filename)
{
    std::ifstream file(filename);
    
    char inter;
    int maxval;
    Image image;
    
     
    if(file.is_open()) 
    {
        file>>inter>>inter>>image.nx>>image.ny>>maxval;
        
        image.image_ =Make_img(image.nx, image.ny);
        
        int nx=image.nx;
        int ny=image.ny;
        double** img=image.image_;
        
        for (int i=0;i<ny;i++)
        {
            for (int j=0;j<nx;j++) 
            {
                file>>img[i][j];
            }
        }
        file.close();
     /*
        for (int i=0;i<ny;i++)
        {
            for (int j=0;j<nx;j++) 
            {
                if(img[i][j]<60) img[i][j]=0;
                if(img[i][j]>160) img[i][j]=255;
                
            }
        }
      */  
        
        return image;  
    }
    else
        exit(0);
        
}

    
int main(int argc,char* argv[])
{
    char* filename;
    if(argc==2)
        filename=argv[1];
    else 
    {
        char fname[100]="rayureA.pgm";
        filename=fname;
    }
        
    int shapeGaus;
    int shapeMedian;
    Image u;   //countour of image
    int Iter;
    double mu;
    double delt;
    
    shapeMedian = 3;
    delt = 0.5;    
    mu = 1;
    shapeGaus = 5;
    Iter = 1;
    
    Image image;
    image=Loadpgm(filename);
    u.nx=image.nx;
    u.ny=image.ny;
    
    
    Median(image.image_,image.nx,image.ny,shapeMedian);
    //Conv(image.image_,image.nx,image.ny,5);
    u.image_=Countour(image,Iter,shapeGaus,delt,mu);
    //Median(u.image_,image.nx,image.ny,shapeMedian);
    u.image_=Countour(u,1,shapeGaus,5,mu);
    
    
    Savepgm(u.image_,image.nx,image.ny,3);
    
    
    return 0;
}

