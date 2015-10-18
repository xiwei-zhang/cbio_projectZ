#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/impex.hxx>

using namespace vigra;

int main (int argc, char ** argv)
{
    vigra::ImageImportInfo imageInfo(argv[1]);
    std::cout<<imageInfo.getFileType() <<std::endl;
    std::cout<<"hello world"<<std::endl;

    return 0;
}
