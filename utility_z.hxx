#ifndef __UTILITY_Z_HXX
#define __UTILITY_Z_HXX

#include <queue>
#include "utilities.hxx"
#include "numerictraits.hxx"
#include "stdimage.hxx"
#include "convolution.hxx"
#include "multi_shape.hxx"
#include <vigra/flatmorphology.hxx>
#include <vigra/matrix.hxx>
#include <vigra/linear_algebra.hxx>
#include <vigra/linear_solve.hxx>
#include <vigra/multi_math.hxx>
//#include <vigra/multi_fft.hxx>



using namespace std;
using namespace vigra;

namespace vigra_mod {

const int nl6[2][6][2] = { { {0,-1},{1,0},{0,1},{-1,1},{-1,0},{-1,-1}},
                            {{1,-1},{1,0},{1,1},{0,1},{-1,0},{0,-1}} };
const int nl8[2][8][2] = { { {1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1} },
                           { {1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1} } };


double diffclock(clock_t t1, clock_t t2){
    double diffticks=t1-t2;
    double diffms = (diffticks*1000)/CLOCKS_PER_SEC;
    return diffms;
}

template <class T1> int round(T1 number){
    double a;
    if (number>=0)
        return (float)modf(number,&a)>0.5?ceil(number):floor(number);
    else
        return (float)abs(modf(number,&a))>0.5?floor(number):ceil(number);
}

    
template <class T1, class S1>
void
imInf(MultiArrayView<2, T1, S1> const & imin1,
      MultiArrayView<2, T1, S1> const & imin2,
      MultiArrayView<2, T1, S1> imout)
{
    for (int k=0; k<imin1.size(); ++k){
        imout[k] = min(imin1[k], imin2[k]);
    }
}// end of function


template <class T1, class S1>
void
imSup(MultiArrayView<2, T1, S1> const & imin1,
      MultiArrayView<2, T1, S1> const & imin2,
      MultiArrayView<2, T1, S1> imout)
{
    for (int k=0; k<imin1.size(); ++k){
        imout[k] = max(imin1[k], imin2[k]);
    }
}// end of function
    
    
    
template <class T1, class S1,
	 class T2, class S2>
void
extendImBorder(MultiArrayView<2, T1, S1> const & src,
		MultiArrayView<2, T2, S2> dist)
{
    dist.init(0);
    int w = src.shape()[0];
    int h = src.shape()[1];
    int wext = w+2;
    int hext = h+2;

    for (int i=0; i<w; ++i){
    	for (int j=0; j<h; ++j){
		dist(i+1,j+1) = src(i,j);
	}
    }
    for (int i=0; i<wext; ++i){
	    dist(i,0) = dist(i,1);
	    dist(i,hext-1) = dist(i,hext-2);
    }
    for (int j=0; j<hext; ++j){
	    dist(0,j) = dist(1,j);
	    dist(wext-1,j) = dist(wext-2,j);
    }
} // end of function

template <class T1, class S1,
	 class T2, class S2>
void
removeImBorder(MultiArrayView<2, T1, S1> const & src,
		MultiArrayView<2, T2, S2> dist)
{
    dist.init(0);
    int w = dist.shape()[0];
    int h = dist.shape()[1];

    for (int i=0; i<w; ++i){
    	for (int j=0; j<h; ++j){
		dist(i,j) = src(i+1,j+1);
	}
    }
} // end of function

template <class T1, class S1>
void
drawImBorder(MultiArrayView<2, T1, S1> src, int v)
{
    int w = src.shape()[0];
    int h = src.shape()[1];

    for (int i=0; i<w; ++i){
	    src(i,0) = v;
	    src(i,h-1) = v;
    }
    for (int j=0; j<h; ++j){
	    src(0,j) = v;
	    src(w-1,j) = v;
    }
} // end of function


template <class T1, class S1>
int *histogram( const MultiArrayView<2, T1, S1> imin ) {
    int size[2] = { imin.shape()[0], imin.shape()[1]};
    int *hist = new int[258];
    memset( hist, 0, sizeof(int) * 258);

    for (int i=0; i<size[0]; ++i){
        for (int j=0; j<size[1]; ++j){
            ++ hist[(int) imin(i,j)];
        }
    }

    int vmin(0), vmax(255);
    bool fmin(true), fmax(true);
    for (int i=0; i<255; ++i){
        if (hist[i]==0 && fmin) ++vmin;
        else fmin = false;
        if (hist[255 - i]==0 && fmax) --vmax;
        else fmax = false;
    }

    hist[256] = vmin;
    hist[257] = vmax;

    return hist;

} // end of function


template <class T1, class S1>
void BIm2MultAr(vigra::FImage src, MultiArrayView<2, T1, S1> dest){
	FImage::Iterator itEnd = src.lowerRight();
	FImage::Iterator itCurrent = src.upperLeft();

    double maxV = 0;
	for (int y=0; itCurrent.y < itEnd.y; ++itCurrent.y, ++y){
		vigra::FImage::Iterator itxC = itCurrent;
		for (int x=0; itxC.x < itEnd.x; ++itxC.x, ++x){
			dest(x,y) = *itxC;
            if (*itxC > maxV) maxV = *itxC;
		}
	}
}




template <class T1, class S1>
void BIm2MultAr(vigra::BImage src, MultiArrayView<2, T1, S1> dest){
	BImage::Iterator itEnd = src.lowerRight();
	BImage::Iterator itCurrent = src.upperLeft();

    double maxV = 0;
	for (int y=0; itCurrent.y < itEnd.y; ++itCurrent.y, ++y){
		vigra::BImage::Iterator itxC = itCurrent;
		for (int x=0; itxC.x < itEnd.x; ++itxC.x, ++x){
			dest(x,y) = *itxC;
            if (*itxC > maxV) maxV = *itxC;
		}
	}
}


template <class T1, class S1>
void MultAr2BIm(MultiArrayView<2, T1, S1> src, vigra::FImage & dest ){
	FImage::Iterator itEnd = dest.lowerRight();
	FImage::Iterator itCurrent = dest.upperLeft();

    double maxV = 0;
	for (int y=0; itCurrent.y < itEnd.y; ++itCurrent.y, ++y){
		vigra::FImage::Iterator itxC = itCurrent;
		for (int x=0; itxC.x < itEnd.x; ++itxC.x, ++x){
			*itxC = src(x,y);
            if (*itxC > maxV) maxV = *itxC;
		}
	}
    // cout<<"D: "<<maxV<<endl;
}



template <class T1, class S1>
void MultAr2BIm(MultiArrayView<2, T1, S1> src, vigra::BImage & dest ){
	BImage::Iterator itEnd = dest.lowerRight();
	BImage::Iterator itCurrent = dest.upperLeft();

    double maxV = 0;
	for (int y=0; itCurrent.y < itEnd.y; ++itCurrent.y, ++y){
		vigra::BImage::Iterator itxC = itCurrent;
		for (int x=0; itxC.x < itEnd.x; ++itxC.x, ++x){
			*itxC = src(x,y);
            if (*itxC > maxV) maxV = *itxC;
		}
	}
    // cout<<"D: "<<maxV<<endl;
}



template <class T1, class S1>
void RecOverBuild( const MultiArrayView<2, T1, S1> immask, const MultiArrayView<2, T1, S1> immark, MultiArrayView<2, T1, S1> imout, int se){
	// Initialize HQ and state image
    int *** nl = new int ** [2];
    if (se == 6){
        for (int k=0; k<2; ++k){
            nl[k] = new int * [6];
            for (int l=0; l<6; ++l){
                nl[k][l] = new int [2];
                nl[k][l][0] = nl6[k][l][0];
                nl[k][l][1] = nl6[k][l][1];
            }
        }
    }
    else if (se == 8){
        for (int k=0; k<2; ++k){
            nl[k] = new int * [8];
            for (int l=0; l<8; ++l){
                nl[k][l] = new int [2];
                nl[k][l][0] = nl8[k][l][0];
                nl[k][l][1] = nl8[k][l][1];
            }
        }
    }
   
    int size[2] = {int(immask.shape()[0]), int(immask.shape()[1])};

	uint16_t x,y,xx,yy; 
    int prio;
    for (int k=0; k<immask.size(); ++k)
        imout[k] = max(immask[k], immark[k]);
	MultiArray<2, UInt8> imstate(immask.shape());
    imstate.init(0);
	queue<uint16_t> HQ[256][2];
    for (int i=0; i<size[0]; ++i){
        for (int j=0; j<size[1]; ++j){
            HQ[(int)imout(i,j)][0].push(i);
            HQ[(int)imout(i,j)][1].push(j);
        }
    }
	
	// Start treat HQ
	for (int i=0; i<256; i++){
		while (!HQ[i][0].empty()){
			x = HQ[i][0].front();
			y = HQ[i][1].front();
			HQ[i][0].pop();
			HQ[i][1].pop();
			imstate(x,y)=2;
			// find neighbor
			for (int k=0; k<se; ++k){
                xx = x + nl[y%2][k][0];
                yy = y + nl[y%2][k][1];
				if (xx<0 || xx>=size[0] || yy<0 || yy>=size[1]) continue;
				if (imstate(xx,yy)==2 ) continue;
				if (imstate(xx,yy)==0){
					imstate(xx,yy) = 1;  // update state image
					prio = max(imout(x,y),immask(xx,yy));
					HQ[prio][0].push(xx);  // put neighbors in the HQ
					HQ[prio][1].push(yy);
					imout(xx,yy) = prio;
				}
			}
		}
	}
	for (int k=0; k<2; ++k){
        for (int l=0; l<se; ++l){
		    delete[] nl[k][l];
        }
	}
	delete[] nl;

}


template <class T1, class S1>
void RecUnderBuild( const MultiArrayView<2, T1, S1> immask, const MultiArrayView<2, T1, S1> immark, MultiArrayView<2, T1, S1> imout, int se){
	// Initialize HQ and state image
    int *** nl = new int ** [2];
    if (se == 6){
        for (int k=0; k<2; ++k){
            nl[k] = new int * [6];
            for (int l=0; l<6; ++l){
                nl[k][l] = new int [2];
                nl[k][l][0] = nl6[k][l][0];
                nl[k][l][1] = nl6[k][l][1];
            }
        }
    }
    else if (se == 8){
        for (int k=0; k<2; ++k){
            nl[k] = new int * [8];
            for (int l=0; l<8; ++l){
                nl[k][l] = new int [2];
                nl[k][l][0] = nl8[k][l][0];
                nl[k][l][1] = nl8[k][l][1];
            }
        }
    }
   
	int size[2] = {int(immask.shape()[0]), int(immask.shape()[1])};
	uint16_t x,y,xx,yy; 
    int prio;
    for (int k=0; k<immask.size(); ++k)
        imout[k] = min(immask[k], immark[k]);
	MultiArray<2, UInt8> imstate(immask.shape());
    imstate.init(0);
	queue<uint16_t> HQ[256][2];
    for (int i=0; i<size[0]; ++i){
        for (int j=0; j<size[1]; ++j){
            HQ[(int)imout(i,j)][0].push(i);
            HQ[(int)imout(i,j)][1].push(j);
        }
    }
	
	// Start treat HQ
	for (int i=255; i>=0; --i){
		while (!HQ[i][0].empty()){
			x = HQ[i][0].front();
			y = HQ[i][1].front();
			HQ[i][0].pop();
			HQ[i][1].pop();
			imstate(x,y)=2;
			// find neighbor
			for (int k=0; k<se; ++k){
                xx = x + nl[y%2][k][0];
                yy = y + nl[y%2][k][1];
				if (xx<0 || xx>=size[0] || yy<0 || yy>=size[1]) continue;
				if (imstate(xx,yy)==2 ) continue;
				if (imstate(xx,yy)==0){
					imstate(xx,yy) = 1;  // update state image
					prio = min(imout(x,y),immask(xx,yy));
					HQ[prio][0].push(xx);  // put neighbors in the HQ
					HQ[prio][1].push(yy);
					imout(xx,yy) = prio;
				}
			}
		}
	}
	for (int k=0; k<2; ++k){
        for (int l=0; l<se; ++l){
		    delete[] nl[k][l];
        }
	}
	delete[] nl;

}



void Maxima( const MultiArrayView<2, UInt8> imin,  MultiArrayView<2, UInt8> imout, int se){
    using namespace vigra::multi_math;
	MultiArray<2, int> imtemp1(imin.shape());
	MultiArray<2, int> imtemp2(imin.shape());
	MultiArray<2, UInt8> imtemp3(imin.shape());
	MultiArray<2, UInt8> imtemp4(imin.shape());
    imtemp1 = imin;
    imtemp2 = imtemp1 -1;
    for (int k=0; k<imin.size(); ++k)
        if (imtemp2[k]<0) imtemp2[k] = 0;
    imtemp3 = imtemp2;

    RecUnderBuild(imin, imtemp3, imtemp4,6);
    imout = imin - imtemp4;
}

void Minima( const MultiArrayView<2, UInt8> imin,  MultiArrayView<2, UInt8> imout, int se){
    using namespace vigra::multi_math;
	MultiArray<2, int> imtemp1(imin.shape());
	MultiArray<2, int> imtemp2(imin.shape());
	MultiArray<2, UInt8> imtemp3(imin.shape());
	MultiArray<2, UInt8> imtemp4(imin.shape());
    imtemp1 = imin;
    imtemp2 = imtemp1 +1;
    for (int k=0; k<imin.size(); ++k)
        if (imtemp2[k]>255) imtemp2[k] = 255;
    imtemp3 = imtemp2;

    RecOverBuild(imin, imtemp3, imtemp4,6);
    imout = imtemp4 - imin;
}



void FillHoles( const MultiArrayView<2, UInt8> imin,  MultiArrayView<2, UInt8> imout, int se, int borderWidth=0){
	MultiArray<2, UInt8> imtemp1(imin.shape());
    MultiArray<2, UInt8> imtemp2(imin.shape());

    imtemp1.init(255);
    drawImBorder(imtemp1, 0);
    RecOverBuild(imin,imtemp1, imout,se);
    
    if (borderWidth > 0){
        imtemp1.init(255);
        MultiArrayView <2, UInt8> subarray = imtemp1.subarray(Shape2(borderWidth, borderWidth), Shape2(-borderWidth, -borderWidth));
        subarray.init(0);
        
        RecOverBuild(imin,imtemp1,imtemp2,se);
        
        imtemp1 = imout;
        imSup(imtemp2, imtemp1, imout);
    }
}


template <class T1, class S1>
void Threshold( const MultiArrayView<2, T1, S1> imin,  MultiArrayView<2, T1, S1> imout, int v_start, int v_end, int v_in, int v_out){
    imout.init(0);
    for (int k=0; k<imin.size(); ++k){
        if (imin[k]>=v_start && imin[k]<=v_end)
            imout[k] = v_in;
        else
            imout[k] = v_out;
    }
}




void Label(const MultiArrayView<2, UInt8> imin, MultiArrayView<2, int> imout, int se){
	/* Vincent's queue algo. scan once */
	// initialization
	// threshold(imin,imout,0,0,0);
    imout.init(0);
	MultiArray<2, UInt8> imflag(imin.shape());
	int size[2] = {int(imin.shape()[0]), int(imin.shape()[1]) };

	queue<int> Qx,Qy;
	int label(0),x,y,s,t;
    int *** nl = new int ** [2];
    if (se == 6){
        for (int k=0; k<2; ++k){
            nl[k] = new int * [6];
            for (int l=0; l<6; ++l){
                nl[k][l] = new int [2];
                nl[k][l][0] = nl6[k][l][0];
                nl[k][l][1] = nl6[k][l][1];
            }
        }
    }
    else if (se == 8){
        for (int k=0; k<2; ++k){
            nl[k] = new int * [8];
            for (int l=0; l<8; ++l){
                nl[k][l] = new int [2];
                nl[k][l][0] = nl8[k][l][0];
                nl[k][l][1] = nl8[k][l][1];
            }
        }
    }


	for(int j=0;j<size[1];++j){
		for(int i=0;i<size[0];++i){
			if (imflag(i,j)==1)
				continue;
			if (imin(i,j)==0 && imflag(i,j)==0){
				imflag(i,j) = 1;
				imout(i,j)= 0;
			}
			else {
				imout(i,j) = ++label;
				imflag(i,j)=1;
				Qx.push(i);
				Qy.push(j);
			}

			while (!Qx.empty()){
				s = Qx.front();
				t = Qy.front();
				Qx.pop();
				Qy.pop();

				for (int k=0; k<se; ++k){
                    x = s + nl[t%2][k][0];
                    y = t + nl[t%2][k][1];
					if (x<0 || x>=size[0] || y<0 || y>=size[1]) continue;
					if (imflag(x,y)==0){
						imflag(x,y)=1;
						if (imin(x,y)!=0){
							imout(x,y) = label;
							Qx.push(x);
							Qy.push(y);
						}
					}
				}
			}
		}
	}
	

	for (int k=0; k<2; ++k){
        for (int l=0; l<se; ++l){
		    delete[] nl[k][l];
        }
	}
	delete[] nl;

}

int labelCount(const MultiArrayView<2, int> imlabel){
	int n(0);
	for(int j=0;j<imlabel.shape()[1];++j){
		for(int i=0;i<imlabel.shape()[0];++i){
			if (imlabel(i,j)>n)
				n = imlabel(i,j);
		}
	}
	return n;
}



void drawEllipse(MultiArrayView<2, UInt8> imin, int x0, int y0, int a, int b, double theta, int fill = 0){
    int Npoint = a*b;
    double x,y;
    int x_, y_;
    double dx, dy;
    double interv = M_PI*2/Npoint;

    if (fill == -1){
        for (int k=0; k<Npoint; ++k){
           x = a*cos( k*interv); 
           y = b*sin( k*interv); 
           x_ = roundf( cos(theta)*x - sin(theta)*y ) + x0;
           y_ = roundf( sin(theta)*x + cos(theta)*y ) + y0;

           if (x_ < 0) x_=0;
           if (x_ >= imin.shape()[0]) x_=imin.shape()[0]-1;
           if (y_ < 0) y_=0;
           if (y_ >= imin.shape()[1]) y_=imin.shape()[1]-1;

           imin(x_,y_) = 255;
        }
        MultiArray<2, UInt8> imfill(imin.shape());
        MultiArray<2, UInt8> imtemp1(imin.shape());
        MultiArray<2, UInt8> imtemp2(imin.shape());
        imtemp1 = imin;
        vigra::discDilation(imtemp1, imtemp2, 1);
        FillHoles(imtemp2, imfill, 8);
        vigra::discErosion(imfill, imin, 1);
    }

    if (fill == 1){
        theta = M_PI - theta;
        Matrix<double> A(2,2);
        Matrix<double> Ainv(2,2);
        Matrix<double> X(2,1);
        Matrix<double> X_(2,1);
        A(0,0) = cos(theta);
        A(0,1) = -sin(theta);
        A(1,0) = sin(theta);
        A(1,1) = cos(theta);
        Ainv = inverse(A);

        for (int i=0; i<imin.shape()[0]; ++i){
            for (int j=0; j<imin.shape()[1]; ++j){
                X_(0,0) = i - x0;
                X_(1,0) = j - y0;
                X = A*X_;

                dx = X(0,0);
                dy = X(1,0);

                if ((dx / a)*(dx / a) + (dy / b)*(dy / b) <= 1)
                    imin(i,j) = 255;
            }
        }
    }


    else{
        for (int k=0; k<Npoint; ++k){
           x = a*cos( k*interv); 
           y = b*sin( k*interv); 

           x_ = roundf( cos(theta)*x - sin(theta)*y ) + x0;
           y_ = roundf( sin(theta)*x + cos(theta)*y ) + y0;

           if (x_ < 0) x_=0;
           if (x_ >= imin.shape()[0]) x_=imin.shape()[0]-1;
           if (y_ < 0) y_=0;
           if (y_ >= imin.shape()[1]) y_=imin.shape()[1]-1;

           imin(x_,y_) = 255;
        }
    }

} // end of function



void colorDeconv(const MultiArrayView<2, RGBValue<UInt8> > imin,  MultiArrayView<2, RGBValue<double> > imout, int stainM = 0){
    MultiArray<2, double> imout1(imin.shape());
    MultiArray<2, double> imout2(imin.shape());
    using namespace vigra::linalg;
    if (stainM == 0){
        double Adata[6] = {0.644211, 0.716556, 0.266844, 0.092789, 0.954111, 0.283111};

        Matrix<double> A(Shape2(3,2));
        for (int i=0; i<6; ++i){
            A[i] = Adata[i];
        }
        Matrix<double> C = inverse(A);

        for (int k=0; k<imin.size(); ++k){
            imout[k].red() = imin[k].red()*C(0,0) + imin[k].green()*C(0,1) + imin[k].blue()*C(0,2) ;
            imout[k].green() = imin[k].red()*C(1,0) + imin[k].green()*C(1,1) + imin[k].blue()*C(1,2) ;
            // imout1[k] = imin[k].red()*C(0,0) + imin[k].green()*C(0,1) + imin[k].blue()*C(0,2) ;
            // imout2[k] = imin[k].red()*C(1,0) + imin[k].green()*C(1,1) + imin[k].blue()*C(1,2) ;
        }

        imout1 = imout.bindElementChannel(0);
        imout2 = imout.bindElementChannel(1);

//        exportImage(imout1, ImageExportInfo("imDeconv1.png"));
//        exportImage(imout2, ImageExportInfo("imDeconv2.png"));
    }
} // end of function


template <class T1, class S1>
void
im2uint8(MultiArrayView<2, T1, S1> const & imin,
		MultiArrayView<2, UInt8> imout)
{
    using namespace vigra::multi_math;

    double vmax(imin[0]), vmin(vmax);
    for (int k=0; k<imin.size(); ++k){
        if (imin[k] > vmax) vmax = imin[k];
        if (imin[k] < vmin) vmin = imin[k];
    }
    if (vmax == vmin)
        imout.init(0);
    else{
        imout = (imin - vmin) / (vmax - vmin) * 255;
    }
}// end of function

    
template <class T1, class S1>
int
meanValue (MultiArrayView<2, T1, S1> const & imin, int onset = 0)
{
    float meanV = 0;
    int nbPixel = 0;
    for (int k=0; k<imin.size(); ++k){
        if (onset){
            if (imin[k] < onset)
                continue;
            else{
                meanV += imin[k];
                nbPixel ++;
            }
        }
        else{
            meanV += imin[k];
        }
    }
    
    if (onset){
        if (nbPixel==0)
            nbPixel = 1;
        meanV /= nbPixel;
    }
    else{
        meanV /= imin.size();
    }
    
    return int(round(meanV));
}// end of function
    
    

    
    
//template <class T1, class S1>
//void
//fastMeanFilter(MultiArrayView<2, T1, S1> const & imin,
//         MultiArrayView<2, T1, S1> imout, int w)
//{
//    // convolve real array with a Gaussian (sigma=1) defined in the spatial domain
//    // (implicitly uses padding by at least 4 pixels)
//    MultiArray<2, double> src(imin.shape()), dest(imin.shape());
//    for (int k=0; k < imin.size(); ++k)
//        src[k] = double(imin[k]);
//    //// create spatial kernel
//    MultiArray<2, double> spatial_kernel(Shape2(w, w));
////    Gaussian<double> gauss(1.0);
//    double kernelValue = 1.0 / (w * w);
//    for(int y=0; y<w; ++y)
//        for(int x=0; x<w; ++x)
//            spatial_kernel(x, y) = kernelValue;
//    convolveFFT(src, spatial_kernel, dest);
//    for (int k=0; k < imin.size(); ++k)
//        imout[k] = UInt8(dest[k]);
//    // convolve real array with a Gaussian (sigma=1) defined in the Fourier domain
//    // (uses no padding, because the kernel size corresponds to the input size)
////    MultiArray<2, FFTWComplex<double> > fourier_kernel(fftwCorrespondingShapeR2C(src.shape()));
////    int y0 = h / 2;
////    
////    for(int y=0; y<fourier_kernel.shape(1); ++y)
////        for(int x=0; x<fourier_kernel.shape(0); ++x)
////            fourier_kernel(x, y) = exp(-0.5*sq(x / double(w))) * exp(-0.5*sq((y-y0)/double(h)));
////    convolveFFT(src, fourier_kernel, dest);
//}// end of function
//
} // namespace vigra_mod



#endif /* FRST_VIGRA_MOD */
