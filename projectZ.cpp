/***
 ProjectZ, coded by Xiwei Zhang in Oct. 2015, CBIO, Paris
 to use the binary file:
 ./projectz [input color image] [output image name] (csv file, for annotation analysis) (features file name)
 ***/

#include <iostream>
#include <assert.h>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>

#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/impex.hxx>
#include <vigra/morpho/morpho_criteria.hxx>

#include <vigra/flatmorphology.hxx>


#include "utility_z.hxx"
#include "TOfft.hxx"

using namespace vigra;

class features{
    int areaRC;
    int areaMid;
    int meanIntensityResRC;
    int maxIntensityResRC;
    int minIntensityResRC;
    int meanIntensityResMid;
    int maxIntensityResMid;
    int minIntensityResMid;
    
    int x_coord;
    int y_coord;
    int x_center;
    int y_center;
    int n_label;
 
public:
    features(){
        areaRC = 0;
        areaMid = 0;
        meanIntensityResRC = 0;
        maxIntensityResRC = 0;
        minIntensityResRC = 255;
        meanIntensityResMid = 0;
        maxIntensityResMid = 0;
        minIntensityResMid = 255;
        
        x_coord = -1;
        y_coord = -1;
        x_center = 0;
        y_center = 0;
        n_label = -1;
    }
    
    void getFeatures(){
        std::cout<<"areaRC "<< areaRC<<std::endl
        <<"areaMid "<< areaMid << std::endl
        <<"meanIntensityResRC "<< meanIntensityResRC << std::endl
        <<"maxIntensityResRC " << maxIntensityResRC << std::endl
        <<"minIntensityResRC " << minIntensityResRC << std::endl
        <<"meanIntensityResMid "<< meanIntensityResMid << std::endl
        <<"maxIntensityResMid " << maxIntensityResMid << std::endl
        <<"minIntensityResMid " << minIntensityResMid << std::endl
        <<"x_coord " << x_coord <<std::endl
        <<"y_coord " << y_coord << std::endl
        <<"x_center " << x_center << std::endl
        <<"y_center " << y_center << std::endl;
    }
    
    void setLabel(int n){
        n_label = n;
    }
    
    void setCoord(int x, int y){
        x_coord = x;
        y_coord = y;
    }
    
    const int & _areaRC() const { return areaRC; }
    const int & _areaMid() const { return areaMid; }
    const int & _meanIntensityResRC() const { return meanIntensityResRC; }
    const int & _maxIntensityResRC() const { return maxIntensityResRC; }
    const int & _minIntensityResRC() const { return minIntensityResRC; }
    const int & _meanIntensityResMid() const { return meanIntensityResMid; }
    const int & _maxIntensityResMid() const { return maxIntensityResMid; }
    const int & _minIntensityResMid() const { return minIntensityResMid; }
    const int & _x_coord() const { return x_coord; }
    const int & _y_coord() const { return y_coord; }
    const int & _x_center() const { return x_center; }
    const int & _y_center() const { return y_center; }
    const int & _n_label() const { return n_label; }

    friend void computerFeats_1( vector<features> & feats, MultiArray<2, UInt8> const iminH, MultiArray<2, UInt8> const imRes, MultiArray<2, UInt8> const imCandi, MultiArray<2, UInt8> const imCandiMid, MultiArray<2, int> const imlabel, bool debugMode );
};


void computerFeats_1( vector<features> & feats, MultiArray<2, UInt8> const iminH, MultiArray<2, UInt8> const imRes, MultiArray<2, UInt8> const imCandi, MultiArray<2, UInt8> const imCandiMid, MultiArray<2, int> const imlabel, bool debugMode ){
    
    int n_candi = vigra_mod::labelCount(imlabel);
    if ( (n_candi + 1) != feats.size()) std::cout<<"SOMETHING IS WRONG!!"<<std::endl;
    
    for (int i=0; i<iminH.shape()[0]; ++i){
        for (int j=0; j<iminH.shape()[1]; ++j){
            //// For Bottom layer
            if (imlabel(i,j)==0) continue;
            int n = imlabel(i,j);
            feats[n].areaRC ++;
            feats[n].meanIntensityResRC += imRes(i,j);
            if (feats[n].maxIntensityResRC < imRes(i,j)) feats[n].maxIntensityResRC = imRes(i,j);
            if (feats[n].minIntensityResRC > imRes(i,j)) feats[n].minIntensityResRC = imRes(i,j);
            if (feats[n].x_coord == -1) feats[n].x_coord = i;
            if (feats[n].y_coord == -1) feats[n].y_coord = j;
            feats[n].x_center += i;
            feats[n].y_center += j;
            
            //// For mid layer
            if (imCandiMid(i,j)==0) continue;
            feats[n].areaMid ++;
            feats[n].meanIntensityResMid += imRes(i,j);
            if (feats[n].maxIntensityResMid < imRes(i,j)) feats[n].maxIntensityResMid = imRes(i,j);
            if (feats[n].minIntensityResMid > imRes(i,j)) feats[n].minIntensityResMid = imRes(i,j);
        }
    }
    
    for (int k=0; k < n_candi; ++k){
        float areaf = float( feats[k].areaRC );
        feats[k].meanIntensityResRC = int( vigra_mod::round( feats[k].meanIntensityResRC / areaf ));
        feats[k].x_center = int( vigra_mod::round( feats[k].x_center / areaf ));
        feats[k].y_center = int( vigra_mod::round( feats[k].y_center / areaf ));
        
        if (feats[k].areaMid != 0){
            areaf = float( feats[k].areaMid );
            feats[k].meanIntensityResMid = int( vigra_mod::round( feats[k].meanIntensityResMid / areaf ));
        }
    }
}


void findMitosis(vector<int> const (& mitosPos)[2],  MultiArray<2, int> const imlabel, vector<int> & mitoLabel, int R){
    const int nb[8][2] = { {-1,-1}, {-1,0}, {-1,1},
        {0,-1}, {0,1}, {1,-1}, {1,0}, {1,1} };
    for (int k = 0; k<mitosPos[0].size(); ++k){
        int x0 = mitosPos[0][k];
        int y0 = mitosPos[1][k];
        
        //// find the mitosis in the candidates
        bool findP = false;
        MultiArray<2, bool> imstate (imlabel.shape());
        imstate.init(false);
        int interCount = 0;
        queue<int> Qxy[2];
        Qxy[0].push(-1);
        Qxy[1].push(-1);
        Qxy[0].push(x0);
        Qxy[1].push(y0);
        imstate(x0, y0) = true;
        
        while (Qxy[0].size() > 0 && interCount < R){
            int x = Qxy[0].front();
            int y = Qxy[1].front();
            Qxy[0].pop();
            Qxy[1].pop();
            
            if (x == -1) {
                interCount ++;
                Qxy[0].push(-1);
                Qxy[1].push(-1);
            }
            else{
                if (imlabel(x,y) > 0){
                    mitoLabel.push_back(imlabel(x,y));
                    findP = true;
                    break;
                }
                else{
                    for (int coord = 0; coord < 8; ++coord){
                        int xx = x + nb[coord][0];
                        int yy = y + nb[coord][1];
                        if (xx < 0 || xx >= imlabel.shape()[0] ||
                            yy < 0 || yy >= imlabel.shape()[1] )
                            continue;
                        if (! imstate(xx, yy)){
                            imstate(xx, yy) = true;
                            Qxy[0].push(xx);
                            Qxy[1].push(yy);
                        }
                    }
                }
            }
        }
        
        if (!findP){
            std::cout<< "Not find mitosis "<<x0<<" "<<y0<<std::endl;
        }
    }
    
//    for (int k = 0; k<mitoLabel.size(); ++k)
//        std::cout<< mitoLabel[k] <<" "<< std::endl;
    
    
}


void imageProc( MultiArray<2, vigra::RGBValue<UInt8> > const & iminRGB, bool const debugMode, int argc, char ** const argv){
    
    using namespace vigra::multi_math;
    
    bool analyseMode = false;
    
    if (argc == 5) analyseMode = true;
    
    MultiArray<2, vigra::RGBValue<double> > imtempRGB(iminRGB.shape());

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
    
    if (debugMode){
        exportImage(iminH, ImageExportInfo("output/iminH.png"));
        exportImage(iminE, ImageExportInfo("output/iminE.png"));
    }
    
    
    //// 1. Fill holes:
    MultiArray<2, UInt8> imtemp1(iminRGB.shape());
    MultiArray<2, UInt8> imtemp2(iminRGB.shape());
    MultiArray<2, UInt8> imtemp3(iminRGB.shape());
    MultiArray<2, UInt8> imRes(iminRGB.shape());
    MultiArray<2, UInt8> imDiv(iminRGB.shape());
    MultiArray<2, UInt8> imDivArea(iminRGB.shape());

    
    // vigra::discMedian(iminH, imtemp3, 2); // PARA
    imtemp3 = iminH;
    if (debugMode)
        exportImage(imtemp3, ImageExportInfo("output/imMedian3.png"));
    
    vigra_mod::FillHoles(imtemp3, imtemp1, 8, 50);
    imRes = imtemp1 - imtemp3;
    if (debugMode)
        exportImage(imRes, ImageExportInfo("output/imResFillholes.png"));
    
    
    //// 1.b Divided by 2 and Red
    imtemp1 = imRes / 2; // PARA
    imtemp2.init(10);  // PARA
    vigra_mod::imSup(imtemp1, imtemp2, imtemp3);
    vigra_mod::RecUnderBuild(imRes, imtemp3, imtemp2, 8);
    imDiv = imRes - imtemp2;
    
    //// 1.c Area opening
    vigra_mod::Threshold(imDiv, imtemp2, 1, 255, 255, 0);
    BasicImage<UInt8> imb1(iminRGB.shape()[0], iminRGB.shape()[1]);
    BasicImage<UInt8> imb2(iminRGB.shape()[0], iminRGB.shape()[1]);
    Diff2D image_size = Diff2D(iminRGB.shape()[0], iminRGB.shape()[1]);
    morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
    
    vigra_mod::MultAr2BIm(imtemp2, imb1);
    morpho::morphoAreaOpening(imb1, imb2, 25, nb);  // PARA
    vigra_mod::BIm2MultAr(imb2, imDivArea);
    if (debugMode){
        exportImage(imDiv, ImageExportInfo("output/imDividedRec.png"));
        exportImage(imDivArea, ImageExportInfo("output/imDividedRecArea.png"));
    }
    
    
    //// 2. Dynamic threshold
    ////     -> mean filter
    ////     -> compare to the filtered image
    MultiArray<2, UInt8> imThd(iminRGB.shape());
    
    fastMeanFilter(imRes, imtemp1, 30);  // PARA
    if (debugMode)
        exportImage(imtemp1, ImageExportInfo("output/imResMean.png"));
    UInt8 globalMeanV = vigra_mod::meanValue(imRes, 5);
    UInt8 tempV(0);
    for (int k=0; k < iminRGB.size(); ++k){
        tempV = std::max<UInt8>(globalMeanV, imtemp1[k]); // minimum value of thd
        imThd[k] = imRes[k] > tempV ? 255:0 ;
    }
    
    exportImage(imThd, ImageExportInfo(argv[2]));
//    exportImage(imDivArea, ImageExportInfo(argv[2]));

    
    
    //// Analyse annotated mitosis
    if (analyseMode){
        
        //// Read csv file to get mitosis center coordinate
        std::ifstream infile(argv[3]);
        std::string line;
        vector<int> mitosPos[2];
        while (std::getline(infile, line))
        {
            stringstream ss( line );
            string field;
            
            for (int i = 0; i<2; ++i){
                std::getline( ss, field, ',' );
                stringstream fs( field );
                int value(-1);
                fs >> value;
                mitosPos[i].push_back(value);
            }
        }
        infile.close();
    
        //// Get coresponding candidate
        MultiArray<2, int> imlabel(iminRGB.shape());

        vigra_mod::Label(imThd, imlabel, 8);
        int n_candi = vigra_mod::labelCount(imlabel);
        
        vector <features> feats(n_candi + 1); // "+ 1" is because in imlabel, start from 1, not zero. Thus feats[0] is Null.

        computerFeats_1( feats, iminH, imRes, imThd, imDivArea, imlabel, debugMode);
        
        vector<int> mitoLabel;
        findMitosis(mitosPos, imlabel, mitoLabel, 20);
        
        // feats[0].getFeatures();
        
        ofstream myfile;
        myfile.open (argv[4]);
        for (int k=1; k < feats.size(); ++k){
            bool isMito(false);
            for (int l=0; l<mitoLabel.size(); ++l){
                if (k == mitoLabel[l]) isMito = true;
            }
            myfile << feats[k]._x_center() <<" "<< feats[k]._y_center()<<" "
            <<feats[k]._x_coord()<<" "<<feats[k]._y_coord()<<" "
            <<feats[k]._areaRC()<<" "<<feats[k]._areaMid()<<" "<<feats[k]._meanIntensityResRC()<<" "
            <<feats[k]._maxIntensityResRC()<<" "<<feats[k]._minIntensityResRC()<<" "
            <<feats[k]._meanIntensityResMid()<<" "<<feats[k]._maxIntensityResMid()<<" "
            <<feats[k]._minIntensityResMid();
            
            if (isMito) myfile << " 1\n";
            else myfile << " 0\n";
        }
        myfile.close();
    }
}



int main (int argc, char ** argv)
{
    bool debugMode = false;
    
    clock_t t1 = clock();


    if (argc < 3){
        std::cout<<"Error: Give an image!! and an output image name!!"<<std::endl;
        return 1;
    }
    else if (argc > 3 && argc < 5){
        std::cout<<"If you want to computer the features, give csv file and output name!!"<<std::endl;
        return 1;
    }
    
    
    //// Open an image, read info
    vigra::ImageImportInfo imageInfo(argv[1]);
    
    if (debugMode){
        std::cout<<imageInfo.getFileType() <<std::endl;
        std::cout << "  pixel type:  " << imageInfo.getPixelType() << std::endl;
        std::cout << "  color image: ";
        if (imageInfo.isColor())    std::cout << "yes (";
        else                        std::cout << "no  (";
        std::cout << "number of channels: " << imageInfo.numBands() << ")"<< std::endl;
    }
    
    assert(imageInfo.isColor());

    MultiArray<2, vigra::RGBValue<UInt8> > iminRGB(imageInfo.shape());

    importImage(imageInfo, iminRGB);

    //// start image processing
    imageProc(iminRGB, debugMode, argc, argv);
    
    
    
    clock_t t2 = clock();

    
    if (debugMode)
        cout<<"Total time: "<<double(vigra_mod::diffclock(t2,t1))<<"ms"<<endl;

    
    return 0;
}
