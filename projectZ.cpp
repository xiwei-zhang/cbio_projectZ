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
#include <list>
#include <set>

#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/impex.hxx>
#include <vigra/morpho/morpho_criteria.hxx>

#include <vigra/flatmorphology.hxx>


#include "utility_z.hxx"
#include "TOfft.hxx"
// #include "maxTree.hxx"

const int nl6[2][6][2] = { { {0,-1},{1,0},{0,1},{-1,1},{-1,0},{-1,-1}},
    {{1,-1},{1,0},{1,1},{0,1},{-1,0},{0,-1}} };
const int nl8[2][8][2] = { { {1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1} },
    { {1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1} } };
const float pi = 3.14159;


using namespace vigra;

float getStd( vector<float> const & vec ){
    float mean = 0, variance(0), M2 = 0;
    size_t n = vec.size();
    for(size_t i=0; i<n; ++i) {
        float delta = vec[i] - mean;
        mean += delta/n;
        M2 += delta*(vec[i] - mean);
        variance = M2/(n - 1);
    }
    return sqrt(variance);
}

class features{
    //// basic features
    int areaRC;  // bottom layer
    int areaMid; // middle layer
    int meanIntensityResRC; // residu image
    int maxIntensityResRC;
    int minIntensityResRC;
    int meanIntensityResMid;
    int maxIntensityResMid;
    int minIntensityResMid;
    
    int meanIRC, maxIRC, minIRC, meanIMid, maxIMid, minIMid; // the same on original image
    
    int volumeRC, volumeMid, pics, borderDepth, StdRC, StdMid, UORC, UOMid;
    
    //// geometric features
    int geoLength;  // geodesic length
    int geoLength2; // L2 length between two furthest points
    int perimeter;
    float circularity; // (4*A) / (pi * D^2)  0 ~ 1
    int geoLengthMid;  // in case of serveral CC in the mid layer, take the largest one
    int geoLength2Mid;
    int perimeterMid;
    float circularityMid;
    
    float ratioLen2; // geoLength2Mid / geoLegnth2, for "U" shape structures
    float ratioFillHolesArea;  // volume of the holes over areaRC
    
    
    //// texture and contexture
    int nb_mid_obj;
    
    float nPic_h, nPic_l, area_h, area_l, volume_h, volume_l;

    
    //// others
    int x_coord;
    int y_coord;
    int x_center;
    int y_center;
    int n_label;
    
    bool iskept;
    
    int x_start;
    int x_end;
    int y_start;
    int y_end;
    
    int p1[2]; // start point
    int p2[2]; // end point

    
public:
    features(int width, int height){
        areaRC = 0;
        areaMid = 0;
        meanIntensityResRC = 0;
        maxIntensityResRC = 0;
        minIntensityResRC = 255;
        meanIntensityResMid = 0;
        maxIntensityResMid = 0;
        minIntensityResMid = 255;
        
        meanIRC = 0; maxIRC = 0; minIRC = 255; meanIMid = 0; maxIMid = 0; minIMid = 255;
        
        volumeRC = 0; volumeMid = 0; pics = 0; borderDepth= 0;
        StdRC = 0; StdMid = 0; UORC = 0; UOMid = 0;
        
        geoLength = 0;
        geoLength2 = 0;
        perimeter = 0;
        circularity = 0.0f;
        geoLengthMid = 0;
        geoLength2Mid = 0;
        perimeterMid = 0;
        circularityMid = 0.0f;
        
        nb_mid_obj = 0;
        nPic_h = 0; nPic_l = 0; area_h = 0; area_l = 0; volume_h = 0; volume_l = 0;
        
        x_coord = -1;
        y_coord = -1;
        x_center = 0;
        y_center = 0;
        n_label = -1;
        
        iskept = true;
        
        x_start = width;
        x_end = 0;
        y_start = height;
        y_end = 0;
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
        <<"y_center " << y_center << std::endl
        <<"x_start " << x_start << std::endl
        <<"x_end " << x_end << std::endl
        <<"y_start " << y_start << std::endl
        <<"y_end " << y_end << std::endl;
        
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
    const bool & _iskept() const { return iskept; }
    const int & _geoLength() const { return geoLength; }
    const int & _geoLength2() const { return geoLength2; }
    const int & _perimeter() const { return perimeter; }
    const float & _circularity() const { return circularity; }
    const int & _nb_mid_obj() const { return nb_mid_obj; }




    friend void computerFeats_1( vector<features> & feats, MultiArrayView<2, UInt8> const iminH, MultiArrayView<2, UInt8> const imRes, MultiArrayView<2, UInt8> const imCandi, MultiArrayView<2, UInt8> const imCandiMid, MultiArrayView<2, int> const imlabel, bool debugMode );
    friend void computerFeats_2( vector<features> & feats, MultiArrayView<2, UInt8> const iminH, MultiArrayView<2, UInt8> const imRes, MultiArrayView<2, UInt8> const imCandi, MultiArrayView<2, UInt8> const imCandiMid, MultiArrayView<2, int> const imlabel, MultiArrayView<2, UInt8> const imMaxima, bool debugMode );
    
    friend void selection_1( vector<features> & feats, MultiArrayView<2, int> const imlabel, MultiArrayView<2, UInt8> imout);
    
    friend int cropGeoLength (features const & feat, MultiArrayView<2, int> const & imlabel, int (& p1)[2], int (& p2)[2], int const label, int const se, int & perimeter);
    friend int cropGeoLengthMid (features & feat, MultiArrayView<2, int> const & imlabel, MultiArrayView<2, int> const & imlabel2,  MultiArrayView<2, UInt8> const & iminH, int (& p1)[2], int (& p2)[2], int const * label2Area, int const label, int const se);
    
    friend void imageProc( MultiArrayView<2, vigra::RGBValue<UInt8> > const & iminRGB, bool const debugMode, int argc, char ** const argv);
    
    friend float getRatioFillHolesArea (features & feat, MultiArrayView<2, int> const & imlabel, MultiArrayView<2, UInt8> const & iminH, MultiArrayView<2, int> const & imCandiFHLabel, MultiArrayView<2, UInt8> const & imMaxima, int const label);
    
    friend void getContextualFeats(features & feat, MultiArrayView<2, int> const & imlabel, MultiArrayView<2, UInt8> const & iminH, MultiArrayView<2, UInt8> const & imRes, MultiArrayView<2, UInt8> const & imMaxima, int const label);

    friend void writeFile(char const * output_name, vector<features> const & feats, vector<int> const & mitoLabel);

};  // end of class definition

int cropGeoLength(features const & feat, MultiArrayView<2, int> const & imlabel, int (& p1)[2], int (& p2)[2], int const label, int const se, int & perimeter){
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
    
    MultiArrayView<2, int> imcrop = imlabel.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));
    MultiArray<2, int> imstate(imcrop.shape());
    imstate.init(-1);
    
    queue<int> Q[2];
    float dist, maxDist(0);
    int mx, my, f(0), len(0),sumX(0),sumY(0);
//    std::list<int> temppp[2];
    
    //// get border pixels push into Q
    for (int i=0; i<imcrop.shape()[0]; ++i) {
        for (int j=0; j<imcrop.shape()[1]; ++j) {
            if (imcrop(i,j) != label) continue;
            bool onEdge = false;
            
            if (i==0 || j==0 || i==imcrop.shape()[0]-1 || j==imcrop.shape()[1]-1){
                onEdge = true;
            }
            for (int k=0; k<se; ++k) {
                int px = i + nl[j%2][k][0];
                int py = j + nl[j%2][k][1];
                if (imcrop(px,py) != label) {
                    onEdge = true;
                    break;
                }
            }
            
            if (onEdge) {
                Q[0].push(i);
                Q[1].push(j);
                sumX += i;
                sumY += j;
            }
            
            imstate(i,j) = 0;
        }
    }
    
    perimeter = int(Q[0].size());
    int x0 = roundf(sumX / float(perimeter));
    int y0 = roundf(sumY / float(perimeter));

    //// get the pixel most far
    while(!Q[0].empty()){
        int px = Q[0].front();
        int py = Q[1].front();
        dist = sqrt((float)(x0-px)*(x0-px) + (y0-py)*(y0-py));
        if (dist>=maxDist){
            maxDist = dist;
            mx = px;
            my = py;
        }
        Q[0].pop();
        Q[1].pop();
    }
    
    Q[0].push(mx);
    Q[1].push(my);
    
    while(!Q[0].empty()){
        mx = Q[0].front();
        my = Q[1].front();
        if (imstate(mx,my)!=0){
            Q[0].pop();
            Q[1].pop();
            continue;
        }
        
        for (int k=0; k<se; ++k){
            int px = mx + nl[my%2][k][0];
            int py = my + nl[my%2][k][1];
            if (px<0 || px>=imcrop.shape()[0] || py<0 || py>=imcrop.shape()[1]) continue;
            if (imcrop(px, py)==label && imstate(px,py)==0){  // see if it's on the edge;
                Q[0].push(px);
                Q[1].push(py);
            }
        }
        imstate(mx, my)= 1;
        Q[0].pop();
        Q[1].pop();
    }
    Q[0].push(mx);
    Q[1].push(my);
    Q[0].push(-1); // -1 is a mark point
    Q[1].push(-1);
    imstate(mx, my) = 2;
    p1[0] = mx;
    p1[1] = my;

    
    // 4. Second propagation
    while(!Q[0].empty()){
        mx = Q[0].front();
        my = Q[1].front();
        
        if (mx == -1) {  // if the mark point pop out, one iteration is done, len ++
            ++len;
            Q[0].pop();
            Q[1].pop();
            if (Q[0].empty()) break;
            Q[0].push(-1);
            Q[1].push(-1);
            mx = Q[0].front();
            my = Q[1].front();
        }
        p2[0] = mx;
        p2[1] = my;
        
        // f = 0;
        for (int k=0; k<se; ++k){
            int px = mx + nl[my%2][k][0];
            int py = my + nl[my%2][k][1];
            if (px<0 || px>=imcrop.shape()[0] || py<0 || py>=imcrop.shape()[1]) continue;
            if (imcrop(px, py)==label && imstate(px, py)==1){  // see if it's on the edge;
                Q[0].push(px);
                Q[1].push(py);
                imstate(px, py) = 2;
            }
        }
        
        Q[0].pop();
        Q[1].pop();
    }
    
    
    for (int k=0; k<2; ++k){
        for (int l=0; l<se; ++l){
            delete[] nl[k][l];
        }
    }
    delete[] nl;
    
    
    return len;

}

int cropGeoLengthMid(features & feat, MultiArrayView<2, int> const & imlabel, MultiArrayView<2, int> const & imlabel2, MultiArrayView<2, UInt8> const & iminH, int (& p1)[2], int (& p2)[2], int const * label2Area, int const label, int const se){
    
    MultiArrayView<2, int> imcrop = imlabel.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));
    MultiArrayView<2, int> imcrop2 = imlabel2.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));
    MultiArrayView<2, UInt8> imcropH = iminH.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));

    set<int> midLabels;
    for (int k=0; k<imcrop.size(); ++k) {
        if (imcrop[k] != label) continue;
        if (imcrop2[k] == 0) continue;
        midLabels.insert( imcrop2[k] );
    }
    
    feat.nb_mid_obj = midLabels.size();
    
    int maxArea = 0;
    int maxLabel = 0;
    std::set<int>::iterator it;
    for (it=midLabels.begin(); it!=midLabels.end(); ++it){
        if (label2Area[ *it ] > maxArea){
            maxArea = label2Area[ *it ];
            maxLabel = *it;
        }
    }

    int len = cropGeoLength(feat, imlabel2, p1, p2, maxLabel, 8, feat.perimeterMid);

    //// circularity
    feat.circularityMid = ( 4 * maxArea ) / ( pi * len * len );

    
    //// variance
    vector<float> pixelValue1, pixelValue2;
    for (int k=0; k<imcrop.size(); ++k) {
        if (imcrop2[k] == maxLabel)
            pixelValue2.push_back( imcropH[k] );
        if (imcrop[k] == label)
            pixelValue1.push_back( imcropH[k] );
    }
    
    feat.StdMid = getStd(pixelValue2);
    feat.StdRC = getStd(pixelValue1);

    return len;

}


float getRatioFillHolesArea (features & feat, MultiArrayView<2, int> const & imlabel, MultiArrayView<2, UInt8> const & iminH, MultiArrayView<2, int> const & imCandiFHLabel, MultiArrayView<2, UInt8> const & imMaxima, int const label){
    
    MultiArrayView<2, int> imCropLabel = imlabel.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));
    MultiArrayView<2, UInt8> imCropH = iminH.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));
    MultiArrayView<2, int> imCropFHLabel = imCandiFHLabel.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));
    MultiArrayView<2, UInt8> imCropMaxima = imMaxima.subarray(Shape2(feat.x_start, feat.y_start), Shape2(feat.x_end, feat.y_end));

    float volume(0.0);
    int newLabel(0);
    //// get fillhole image label
    for (int k=0; k<imCropLabel.size(); ++k) {
        if (imCropLabel[k] != label) continue;
        if (imCropFHLabel[k] != 0){
            newLabel = imCropFHLabel[k];
            break;
        }
    }
    
    for (int k=0; k<imCropLabel.size(); ++k) {
        if (imCropFHLabel[k] != newLabel) continue;
        if (imCropLabel[k] == 0){
            volume += imCropH[k];
        }
    }
    
    //// get pics
    for (int k=0; k<imCropLabel.size(); ++k) {
        if (imCropLabel[k] != label) continue;
        if (imCropMaxima[k] > 0)
            feat.pics ++;
    }
    
    return ( volume / feat.areaRC ) ;
}


void getContextualFeats(features & feat, MultiArrayView<2, int> const & imlabel, MultiArrayView<2, UInt8> const & iminH, MultiArrayView<2, UInt8> const & imRes, MultiArrayView<2, UInt8> const & imMaxima, int const label){
    int W = 50; // PARA  // window size 2*W x 2*W
    int th_high = feat.maxIRC; // it's dark object. th_high will include more noise, it correspond to feat_l
    int th_mid = feat.meanIRC;
    
    int x_start = std::max<int>(0, (feat.x_center - W));
    int x_end = std::min<int>((imlabel.shape()[0]-1), (feat.x_center + W) );
    int y_start = std::max<int>(0, (feat.y_center - W) );
    int y_end = std::min<int>((imlabel.shape()[1]-1), (feat.y_center + W) );

    
    MultiArrayView<2, int> imCropLabel = imlabel.subarray(Shape2(x_start, y_start), Shape2(x_end, y_end));
    MultiArrayView<2, UInt8> imCropH = iminH.subarray(Shape2(x_start, y_start), Shape2(x_end, y_end));
    MultiArrayView<2, UInt8> imCropRes = imRes.subarray(Shape2(x_start, y_start), Shape2(x_end, y_end));
    MultiArrayView<2, UInt8> imCropMaxima = imMaxima.subarray(Shape2(x_start, y_start), Shape2(x_end, y_end));

    MultiArray<2, UInt8> imtemp1(imCropH.shape());
    MultiArray<2, UInt8> imtemp2(imCropH.shape());
    vigra_mod::Threshold(imCropH, imtemp1, 0, th_high, 255, 0);
    vigra_mod::Threshold(imCropH, imtemp2, 0, th_mid, 255, 0);

    for (int k=0; k<imCropLabel.size(); ++k) {
        if (imCropLabel[k] == label) continue;
        if (imtemp1[k] !=0){
            feat.area_l ++ ;
            feat.volume_l += imCropRes[k];
            if (imCropMaxima[k] > 0)
                feat.nPic_l ++;
        }
        
        if (imtemp2[k] !=0){
            feat.area_h ++ ;
            feat.volume_h += imCropRes[k];
            if (imCropMaxima[k] > 0)
                feat.nPic_h ++;
        }
    }
    
    feat.area_l /= feat.areaRC;
    feat.volume_l /= feat.areaRC;
    feat.area_h /= feat.areaRC;
    feat.volume_h /= feat.areaRC;
}


void computerFeats_1( vector<features> & feats, MultiArrayView<2, UInt8> const iminH, MultiArrayView<2, UInt8> const imRes, MultiArrayView<2, UInt8> const imCandi, MultiArrayView<2, UInt8> const imCandiMid, MultiArrayView<2, int> const imlabel, bool debugMode ){
    
    int n_candi = vigra_mod::labelCount(imlabel);
    if ( (n_candi + 1) != feats.size()) std::cout<<"SOMETHING IS WRONG!!"<<std::endl;
    
    for (int i=0; i<iminH.shape()[0]; ++i){
        for (int j=0; j<iminH.shape()[1]; ++j){
            //// For Bottom layer
            if (imlabel(i,j)==0) continue;
            int n = imlabel(i,j);
            feats[n].areaRC ++;
            feats[n].meanIntensityResRC += imRes(i,j);
            feats[n].meanIRC += iminH(i,j);
            if (feats[n].maxIntensityResRC < imRes(i,j)) feats[n].maxIntensityResRC = imRes(i,j);
            if (feats[n].minIntensityResRC > imRes(i,j)) feats[n].minIntensityResRC = imRes(i,j);
            if (feats[n].maxIRC < iminH(i,j)) feats[n].maxIRC = iminH(i,j);
            if (feats[n].minIRC > iminH(i,j)) feats[n].minIRC = iminH(i,j);
            if (feats[n].x_coord == -1) feats[n].x_coord = i;
            if (feats[n].y_coord == -1) feats[n].y_coord = j;
            if (feats[n].x_start > i) feats[n].x_start = i;
            if (feats[n].x_end < i) feats[n].x_end = i;
            if (feats[n].y_start > j) feats[n].y_start = j;
            if (feats[n].y_end < j) feats[n].y_end = j;
            feats[n].x_center += i;
            feats[n].y_center += j;
            
            
            //// For mid layer
            if (imCandiMid(i,j)==0) continue;
            feats[n].areaMid ++;
            feats[n].meanIntensityResMid += imRes(i,j);
            feats[n].meanIMid += iminH(i,j);
            if (feats[n].maxIntensityResMid < imRes(i,j)) feats[n].maxIntensityResMid = imRes(i,j);
            if (feats[n].minIntensityResMid > imRes(i,j)) feats[n].minIntensityResMid = imRes(i,j);
            if (feats[n].maxIMid < iminH(i,j)) feats[n].maxIMid = iminH(i,j);
            if (feats[n].minIMid > iminH(i,j)) feats[n].minIMid = iminH(i,j);
        }
    }
    
    for (int k=0; k < n_candi; ++k){
        feats[k].volumeRC = feats[k].meanIntensityResRC;
        feats[k].volumeMid = feats[k].meanIMid;
        float areaf = float( feats[k].areaRC );
        float areaMidf = float( feats[k].areaMid );
        feats[k].meanIntensityResRC = int( vigra_mod::round( feats[k].meanIntensityResRC / areaf ));
        feats[k].meanIRC = int( vigra_mod::round( feats[k].meanIRC / areaf ));
        feats[k].meanIMid = int( vigra_mod::round( feats[k].meanIMid / areaMidf ));
        feats[k].x_center = int( vigra_mod::round( feats[k].x_center / areaf ));
        feats[k].y_center = int( vigra_mod::round( feats[k].y_center / areaf ));
        
        if (feats[k].areaMid != 0){
            areaf = float( feats[k].areaMid );
            feats[k].meanIntensityResMid = int( vigra_mod::round( feats[k].meanIntensityResMid / areaf ));
        }
    }
}

void computerFeats_2( vector<features> & feats, MultiArrayView<2, UInt8> const iminH, MultiArrayView<2, UInt8> const imRes, MultiArrayView<2, UInt8> const imCandi, MultiArrayView<2, UInt8> const imCandiMid, MultiArrayView<2, int> const imlabel, MultiArrayView<2, UInt8> const imMaxima, bool debugMode ){

    using namespace vigra::multi_math;

    MultiArray<2, int> imlabel2(iminH.shape());
    vigra_mod::Label(imCandiMid, imlabel2, 8);
    int n_candi2 = vigra_mod::labelCount(imlabel2) + 1;
    
    //// Preparation: get area of each label2
    int * label2Area = new int [ n_candi2];
    std::fill_n(label2Area, n_candi2, 0);
    for (int k=0; k<imlabel2.size(); ++k) {
        if (imlabel2[k]==0) continue;
        label2Area[imlabel2[k]] ++;
    }
    
    //// Preparation: get fillholes of candidates and label image
    MultiArray<2, UInt8> imCandiFH(iminH.shape());
    MultiArray<2, int> imCandiFHLabel(iminH.shape());

    vigra_mod::FillHoles(imCandi, imCandiFH, 8, 50);
    vigra_mod::Label(imCandiFH, imCandiFHLabel, 8);

    
    //// Preparation: get VAR
    // MultiArray<2, UInt8> imVAR(iminH.shape());
    /// getVAR(imRes, imVAR, 10);

    
    if (debugMode){
        exportImage(imCandiFH, ImageExportInfo("output/imCandiFillHoles.png"));
    }
    
    
    //// computer geo features and other features
    for (int k=1; k<feats.size(); ++k) {
        int p1[2], p2[2];
        if (!feats[k]._iskept()) continue;
        feats[k].geoLength = cropGeoLength(feats[k], imlabel, p1, p2, k, 8, feats[k].perimeter);
        feats[k].p1[0] = p1[0];
        feats[k].p1[1] = p1[1];
        feats[k].p2[0] = p2[0];
        feats[k].p2[1] = p2[1];
        feats[k].geoLength2 = sqrt(pow(float(p1[0]-p2[0]),2) + pow(float(p1[1]-p2[1]),2));
        feats[k].geoLengthMid = cropGeoLengthMid(feats[k], imlabel, imlabel2, iminH, p1, p2, label2Area, k, 8);
        feats[k].geoLength2Mid = sqrt(pow(float(p1[0]-p2[0]),2) + pow(float(p1[1]-p2[1]),2));
        feats[k].circularity = ( 4 * feats[k].areaRC ) / ( pi * feats[k].geoLength * feats[k].geoLength );
        feats[k].ratioLen2 = float(feats[k].geoLength2Mid) / feats[k].geoLength2;
        feats[k].ratioFillHolesArea = getRatioFillHolesArea(feats[k], imlabel, iminH, imCandiFHLabel, imMaxima, k);
        
        getContextualFeats(feats[k], imlabel, iminH, imRes, imMaxima, k);
    }
    
    delete[] label2Area;
}


void selection_1( vector<features> & feats, MultiArrayView<2, int> const imlabel, MultiArrayView<2, UInt8> imout){
    
    int * mapLabelIntensity = new int [ feats.size() ];
    std::fill_n(mapLabelIntensity, feats.size(), 0);

    for (int k=1; k < feats.size(); ++k){
        // 47 70 30 2000
        if (feats[k].meanIntensityResRC < 36 || feats[k].areaRC < 70 || feats[k].areaMid < 30 || feats[k].areaRC > 2000) {            feats[k].iskept = false;
        }
        else {
            mapLabelIntensity[k] = 255;
        }
    }
    
    for (int k=0; k < imlabel.size(); ++k){
        imout[k] = mapLabelIntensity[ imlabel[k] ];
    }
}


void findMitosis(vector<int> const (& mitosPos)[2],  MultiArrayView<2, int> const imlabel, vector<int> & mitoLabel, int R){
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


void writeFile(char const * output_name, vector<features> const & feats, vector<int> const & mitoLabel){
    ofstream myfile;
    myfile.open (output_name);
    for (int k=1; k < feats.size(); ++k){
        if (!feats[k]._iskept()) continue;
        bool isMito(false);
        for (int l=0; l<mitoLabel.size(); ++l){
            if (k == mitoLabel[l]) isMito = true;
        }
        myfile << feats[k]._x_center() <<" "<< feats[k]._y_center()<<" "
        <<feats[k]._x_coord()<<" "<<feats[k]._y_coord()<<" "
        <<feats[k]._areaRC()<<" "<<feats[k]._areaMid()<<" "<<feats[k]._meanIntensityResRC()<<" "
        <<feats[k]._maxIntensityResRC()<<" "<<feats[k]._minIntensityResRC()<<" "
        <<feats[k]._meanIntensityResMid()<<" "<<feats[k]._maxIntensityResMid()<<" "
        <<feats[k]._minIntensityResMid()<<" "<<feats[k]._geoLength()<<" "<<feats[k]._geoLength2()<<" "
        <<feats[k]._perimeter()<<" "<< feats[k].geoLengthMid <<" "<< feats[k].geoLength2Mid <<" "
        <<feats[k].perimeterMid <<" "<<feats[k]._circularity()<<" "<<feats[k].circularityMid<<" "
        <<feats[k]._nb_mid_obj()<<" "<<feats[k].ratioFillHolesArea;
        
        if (isMito) myfile << " 1\n";
        else myfile << " 0\n";
        
        if (isMito) std::cout<<"MITOS : "<<k<<std::endl;
    }
    myfile.close();
}


void imageProc( MultiArrayView<2, vigra::RGBValue<UInt8> > const & iminRGB, bool const debugMode, int argc, char ** const argv){
    
    using namespace vigra::multi_math;
    
    bool analyseMode = false;
    
    if (argc == 5) analyseMode = true;
    
    int width = int(iminRGB.shape()[0]);
    int height = int(iminRGB.shape()[1]);

    
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
    
    if (debugMode){
        exportImage(imThd, ImageExportInfo("output/imThreshold.png"));
    }
//    exportImage(imThd, ImageExportInfo(argv[2]));
//    exportImage(imDivArea, ImageExportInfo(argv[2]));


    //// Get coresponding candidate
    MultiArray<2, int> imlabel(iminRGB.shape());
    MultiArray<2, int> imlabel2(iminRGB.shape());
    MultiArray<2, int> imlabelMaxi(iminRGB.shape());
    MultiArray<2, UInt8> imSelection(iminRGB.shape());
    MultiArray<2, UInt8> imMaxima(iminRGB.shape());

    //// Computer Basic features
    vigra_mod::Label(imThd, imlabel, 8);
    int n_candi = vigra_mod::labelCount(imlabel);
    vector <features> feats(n_candi + 1, features(width, height)); // "+ 1" is because in imlabel, start from 1, not zero. Thus feats[0] is Null.
    computerFeats_1( feats, iminH, imRes, imThd, imDivArea, imlabel, debugMode);
    selection_1( feats, imlabel, imSelection );
    
    
    //// get maxima image (one pixel for each maxima)
    vigra_mod::Maxima(imRes, imMaxima, 8);
    vigra_mod::Label(imMaxima, imlabelMaxi, 8);
    int N_maxi = vigra_mod::labelCount(imlabelMaxi);
    int * maxiLabel = new int [N_maxi + 1];
    std::fill_n(maxiLabel, (N_maxi + 1), 0);
    for (int k = 0; k<imlabelMaxi.size(); ++k) {
        if (imlabelMaxi[k] == 0) continue;
        if (maxiLabel[imlabelMaxi[k]] == 0){
            maxiLabel[imlabelMaxi[k]] = 1;
            imMaxima[k] = 255;
        }
        else
            imMaxima[k] = 0;
    }

    if (debugMode) {
        exportImage(imMaxima, ImageExportInfo("output/imMaxima.png"));
    }
    
    
    //// Computer advanced features
    computerFeats_2( feats, iminH, imRes, imThd, imDivArea, imlabel, imMaxima, debugMode);

    
    if (debugMode){
        exportImage(imSelection, ImageExportInfo("output/imSelection.png"));
    }
    
    exportImage(imSelection, ImageExportInfo(argv[2]));
    

    
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
    
        for (int k=0; k<imlabel.size(); ++k) {
            imlabel2[k] = imSelection[k] > 0 ? imlabel[k] : 0;
        }
        vector<int> mitoLabel;
        findMitosis(mitosPos, imlabel2, mitoLabel, 20);
        
        // feats[0].getFeatures();
        
        writeFile(argv[4], feats, mitoLabel);
    }
    
    delete [] maxiLabel;
}



int main (int argc, char ** argv)
{
    bool debugMode = true;
    
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
