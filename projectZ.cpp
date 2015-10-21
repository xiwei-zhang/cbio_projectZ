#include <iostream>
#include <assert.h>
#include <time.h>
#include <algorithm>

#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/impex.hxx>

#include <vigra/flatmorphology.hxx>


#include "utility_z.hxx"
#include "TOfft.hxx"

using namespace vigra;

int main (int argc, char ** argv)
{
    clock_t t1 = clock();

    using namespace vigra::multi_math;

    if (argc < 2){
        std::cout<<"Error: Give an image!!"<<std::endl;
        return 1;
    }
    
    //// Open an image, read info
    vigra::ImageImportInfo imageInfo(argv[1]);
    
    std::cout<<imageInfo.getFileType() <<std::endl;
    std::cout << "  pixel type:  " << imageInfo.getPixelType() << std::endl;
    std::cout << "  color image: ";
    if (imageInfo.isColor())    std::cout << "yes (";
    else                        std::cout << "no  (";
    std::cout << "number of channels: " << imageInfo.numBands() << ")"<< std::endl;
    
    assert(imageInfo.isColor());

    MultiArray<2, vigra::RGBValue<UInt8> > iminRGB(imageInfo.shape());
    MultiArray<2, vigra::RGBValue<double> > imtempRGB(imageInfo.shape());

    importImage(imageInfo, iminRGB);
    
    
    //// Color deconvolution
    MultiArray<2, double> imtempD1(iminRGB.shape());
    MultiArray<2, double> imtempD2(iminRGB.shape());
    MultiArray<2, UInt8> iminH(iminRGB.shape());
    MultiArray<2, UInt8> iminE(iminRGB.shape());
    
    vigra_mod::colorDeconv(iminRGB, imtempRGB);

    imtempD1 = imtempRGB.bindElementChannel(0);
    imtempD2 = imtempRGB.bindElementChannel(1);
    
    vigra_mod::im2uint8(imtempD1, iminH);
    vigra_mod::im2uint8(imtempD2, iminE);

    exportImage(iminH, ImageExportInfo("output/iminH.png"));
    exportImage(iminE, ImageExportInfo("output/iminE.png"));

    
    //// 1. Fill holes:
    MultiArray<2, UInt8> imtemp1(iminRGB.shape());
    MultiArray<2, UInt8> imtemp2(iminRGB.shape());
    MultiArray<2, UInt8> imtemp3(iminRGB.shape());
    MultiArray<2, UInt8> imRes(iminRGB.shape());


    // vigra::discMedian(iminH, imtemp3, 2); // PARA
    imtemp3 = iminH;
    exportImage(imtemp3, ImageExportInfo("output/imMedian3.png"));

    vigra_mod::FillHoles(imtemp3, imtemp1, 8);
    imRes = imtemp1 - imtemp3;
    exportImage(imRes, ImageExportInfo("output/imResFillholes.png"));

    

    //// 2. Threshold
    MultiArray<2, UInt8> imThd(iminRGB.shape());
    
    fastMeanFilter(imRes, imtemp1, 50);  // PARA
    exportImage(imtemp1, ImageExportInfo("output/imResMean.png"));
    UInt8 globalMeanV = vigra_mod::meanValue(imRes, 1);
    UInt8 tempV(0);
    for (int k=0; k < iminRGB.size(); ++k){
        tempV = std::max<UInt8>(globalMeanV, imtemp1[k]); // minimum value of thd
        imThd[k] = imRes[k] > tempV ? 255:0 ;
    }
    
    exportImage(imThd, ImageExportInfo("output/imThd.png"));

    float meanImV = 0;
    int nbPixel = 0;
    for (int k=0; k < iminRGB.size(); ++k){
        if (imRes[k]==0) continue;
        meanImV += imRes[k];
        nbPixel ++;
    }
    meanImV /= nbPixel;
    std::cout<<"Global mean value: "<< meanImV<<" "<<vigra_mod::meanValue(imRes, 1)<<std::endl;
    
    
    clock_t t2 = clock();

    cout<<"Total time: "<<double(vigra_mod::diffclock(t2,t1))<<"ms"<<endl;

    std::cout<<"Hello world"<<std::endl;
    
    return 0;
}
