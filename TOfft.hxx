#ifndef __TOfft_hxx
#define __TOfft_hxx

#include <complex>
#include <math.h>

#include "utility_z.hxx"

using namespace std;

#define PI 3.14159

void FFT(complex<double> *TD, complex<double> *FD, int r){
	long count;
	int bfsize,p;
	double angle;
	complex<double> *W,*X1,*X2,*X;

	count = 1 << r;

	W = new complex<double>[count/2];
	X1 = new complex<double>[count];
	X2 = new complex<double>[count];

	// cal pow coef
	for (int i=0; i<count/2; i++){
		angle = -i * PI * 2 / count;
		W[i] = complex<double> (cos(angle), sin(angle));
	}

	// real
	memcpy(X1, TD, sizeof(complex<double>) * count);

	// butterfly 
	for (int k=0; k<r;k++){
		for (int j=0; j<1<<k; j++){
			bfsize = 1<<(r-k);
			for (int i=0; i<bfsize/2; i++){
				p = j * bfsize;
				X2[i+p] = X1[i+p] + X1[i+p+bfsize/2];
				X2[i+p+bfsize/2] = (X1[i+p]-X1[i+p+bfsize/2]) * W[i*(1<<k)];
			}
		}
		X = X1;
		X1 = X2;
		X2 = X;
	}

	// order
	for (int j=0; j< count; j++){
		p = 0;
		for(int i=0; i<r; i++){
			if (j&(1<<i)){
				p+=1<<(r-i-1);
			}
		}
		FD[j]=X1[p];
	}

	delete W;
	delete X1;
	delete X2;

}

void IFFT(complex<double> *FD, complex<double> *TD, int r){
	long count;
	complex<double> *X;
	count = 1<<r;

	X = new complex<double>[count];
	memcpy(X,FD,sizeof(complex<double>)*count);

	//conjugate
	for(int i=0; i<count; i++){
		X[i] = complex<double> (X[i].real(), -X[i].imag());
	}

	FFT(X, TD, r);

	// conjugate in time
	for (int i=0; i<count; i++){
		TD[i] = complex<double> (TD[i].real()/count, -TD[i].imag()/count);
	}

	delete X;
}

//void FFT2D(Mat imin){
//	long lWidth = imin.cols;
//	long lHeight = imin.rows;
//
//	long w,h; //fft width and hight. should be pow of 2
//	int wp,hp;
//
//	w=1; h=1; wp=0; hp=0;
//	// calculate the width and height of FFT trans
//	while (w*2 <= lWidth){
//		w *=2;
//		wp ++;
//	}
//	while (h*2 <= lHeight){
//		h *=2;
//		hp ++;
//	}
//
//	// padding image
//	Mat imPad = Mat::zeros(h,w,CV_8U);
//	Mat imFFT = Mat::zeros(h,w,CV_32F);
//	for (int j=0; j<imin.rows; j++)
//		for (int i=0; i<imin.cols; i++)
//			imPad.at<uchar>(j,i) = imin.at<uchar>(j,i);
//
//	complex<double> *TD = new complex<double>[w*h];
//	complex<double> *FD = new complex<double>[w*h];
//
//	for (int j=0; j<h; j++){
//		for (int i=0; i<w; i++){
//			TD[i+j*w]= complex<double>(imPad.at<uchar>(j,i),0);
//		}
//	}
//
//	// y direction FFT
//	for(int j=0; j<h; j++) FFT(&TD[w*j], &FD[w*j], wp);
//	// transform
//	for(int j=0; j<h; j++)
//		for(int i=0; i<w; i++)
//			TD[j+h*i] = FD[i+w*j];
//	// x direction FFT
//	for (int i=0; i<w; i++) FFT(&TD[h*i], &FD[h*i], hp);
//
//
//	// IFFT
//	for (int i=0; i<w; i++) IFFT(&FD[h*i], &TD[h*i], hp);
//	// transform
//	for(int j=0; j<h; j++)
//		for(int i=0; i<w; i++)
//			FD[i+w*j] = TD[j+h*i];
//	for(int j=0; j<h; j++) IFFT(&FD[w*j], &TD[w*j], wp);
//
//	for (int j=0; j<h; j++){
//		for (int i=0; i<w; i++){
//			imPad.at<uchar>(j,i) = round(TD[i+j*w].real());
//		}
//	}
//
//
//	delete[] TD;
//	delete[] FD;
//}

template <class T1, class S1>
void
fastMeanFilter(MultiArrayView<2, T1, S1> const & imin,
               MultiArrayView<2, T1, S1> imout, int size){
//void fastMeanFilter(Mat imin, int size, Mat imout, int depth= 8){
	// depth 8->uchar    33->float
	// initialization /////////////////////////////////////////////////
	long lWidth = imin.shape()[0];
	long lHeight = imin.shape()[1];
	long w,h; //fft width and hight. should be pow of 2
	int wp,hp;

	w=1; h=1; wp=0; hp=0;
	// calculate the width and height of FFT trans
	while (w <= lWidth){
		w *=2;
		wp ++;
	}
	while (h <= lHeight){
		h *=2;
		hp ++;
	}
	/////////////////////////////////////////////////



	// padding image/////////////////////////////////////////////////
    MultiArray<2, double> imPad(w, h);
    MultiArray<2, double> imKernalPad(w, h);

    int meanV = vigra_mod::meanValue(imin,1);
    for (int k=0; k<imPad.size(); ++k){
        imPad[k] = meanV;
    }
	int xx0 = lWidth/2, yy0 = lHeight/2, x0_=w/2, y0_=h/2;
	for (int j=0; j<lHeight; j++){
		for (int i=0; i<lWidth; i++){
            imPad(i-xx0+x0_, j-yy0+y0_) = imin(i, j);
//			if (depth==8) imPad.at<float>(j-yy0+y0_,i-xx0+x0_) = imin.at<uchar>(j,i);
//			if (depth==33) imPad.at<float>(j-yy0+y0_,i-xx0+x0_) = imin.at<float>(j,i);
		}
	}
	int dwKernal = size/2;
	for (int j=-dwKernal; j<=dwKernal; j++){
		for (int i=-dwKernal; i<=dwKernal; i++){
//			imKernalPad.at<float>(h/2+j,w/2+i) = 1.0f/(size*size);
            imKernalPad(w/2+i, h/2+j) = 1.0 / (size * size);
		}
	}


	/////////////////////////////////////////////////



	// FFT/////////////////////////////////////////////////
	complex<double> *TD = new complex<double>[w*h];
	complex<double> *FD = new complex<double>[w*h];
	complex<double> *kTD = new complex<double>[w*h];
	complex<double> *kFD = new complex<double>[w*h];

	for (int j=0; j<h; j++){
		for (int i=0; i<w; i++){
//			TD[i+j*w]= complex<double>(imPad.at<float>(j,i),0);
//			kTD[i+j*w]= complex<double>(imKernalPad.at<float>(j,i),0);
            TD[i+j*w]= complex<double>(imPad(i,j),0);
            kTD[i+j*w]= complex<double>(imKernalPad(i,j),0);
		}
	}
	// y direction FFT
	for(int j=0; j<h; j++) FFT(&TD[w*j], &FD[w*j], wp);
	for(int j=0; j<h; j++) FFT(&kTD[w*j], &kFD[w*j], wp);
	// transform
	for(int j=0; j<h; j++){
		for(int i=0; i<w; i++){
			TD[j+h*i] = FD[i+w*j];
			kTD[j+h*i] = kFD[i+w*j];
		}
	}
	// x direction FFT
	for (int i=0; i<w; i++) FFT(&TD[h*i], &FD[h*i], hp);
	for (int i=0; i<w; i++) FFT(&kTD[h*i], &kFD[h*i], hp);
	/////////////////////////////////////////////////



	// Product in frequency domain//////////////////
	for (int i=0; i<h*w; i++)
		FD[i]= FD[i]*kFD[i];
	/////////////////////////////////////////////////


	///IFFT ////////////////////////////////////////
	// IFFT x
	for (int i=0; i<w; i++) IFFT(&FD[h*i], &TD[h*i], hp);
	// transform
	for(int j=0; j<h; j++)
		for(int i=0; i<w; i++)
			FD[i+w*j] = TD[j+h*i];
	// IFFT y
	for(int j=0; j<h; j++) IFFT(&FD[w*j], &TD[w*j], wp);


	for (int j=0; j<h; j++){
		for (int i=0; i<w; i++){
//            imPad.at<float>(j,i) = round(TD[(i<w/2?w/2+i:i-w/2)+(j<h/2?h/2+j:j-h/2)*w].real());
            imPad(i,j) = vigra_mod::round(TD[(i<w/2?w/2+i:i-w/2)+(j<h/2?h/2+j:j-h/2)*w].real());
		}
	}
	
	for (int j=0; j<lHeight; j++){
		for (int i=0; i<lWidth; i++){
            imout(i,j) = imPad(i-xx0+x0_, j-yy0+y0_);
//			if (depth==8) imout.at<uchar>(j,i) = imPad.at<float>(j-yy0+y0_,i-xx0+x0_);
//			if (depth==33) imout.at<float>(j,i) = imPad.at<float>(j-yy0+y0_,i-xx0+x0_);
		}
	}
	/////////////////////////////////////////////////


	delete[] TD;
	delete[] FD;
	delete[] kTD;
	delete[] kFD;
}


//void getVAR(Mat imin, Mat imout, int l){
//
//	Mat imtemp1 = imin.clone();
//	int w(0);
//	w = ((l-1)/2)*2+1;
//	blur(imin,imtemp1,Size(w,w));
//	imwrite("zz1.png",imtemp1);
//
//	imout.setTo(0);
//	w = (l-1)/2;
//	int x,y,np;
//	np = (w*2+1)*(w*2+1);
//	float meanv, sumv;
//	for (int j=0; j<imin.rows; j++){
//		for (int i=0; i<imin.cols; i++){
//			meanv = imtemp1.at<uchar>(j,i);
//			sumv = 0;
//			for (int m=-w; m<=w; m++){
//				for (int n=-w; n<=w; n++){
//					x = i+m;
//					y = j+n;
//					if (x<0 || x>=imin.cols || y<0 || y>=imin.rows) continue;
//					sumv += ((float)imin.at<uchar>(y,x) - meanv)*((float)imin.at<uchar>(y,x) - meanv);
//				}
//			}
//			sumv /= np;
//			if (sumv>255) sumv=255;
//			imout.at<uchar>(j,i)=round(sumv);
//		}
//	}
//
//}

#endif