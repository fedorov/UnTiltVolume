#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"
#include "itkVector.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"

#include "UnTiltVolumeCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  typedef    T PixelType;

  typedef itk::Image<PixelType, 3> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typedef itk::Vector<float, 3> VectorType;
  typedef itk::AffineTransform<float,3> TransformType;
  typedef itk::TransformFileWriter TransformWriterType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  reader->Update();  

  typename ImageType::Pointer image = reader->GetOutput();
  typename ImageType::DirectionType dir = image->GetDirection();
  std::cout << "Direction: " << dir << std::endl;

  //
  // find the closest direction for each of the major axes
  //
 
  // define major directions
  double directions [6][3] = {
                   {  1,  0,  0 },   // right
                   { -1,  0,  0 },   // left
                   {  0,  1,  0 },   // anterior
                   {  0, -1,  0 },   // posterior
                   {  0,  0,  1 },   // superior
                   {  0,  0, -1 } }; // inferior
  
  int closestAxis[3] = {0., 0., 0.};
  double closestDot[3] = {-1., -1., -1.};

  int direction;
  for (int direction = 0; direction < 6; direction++)
    {
    double dot[3];
    for (int col = 0; col < 3; col++)
      {
      dot[col] = 0;
      int i;
      for (int i = 0; i < 3; i++)
        {
        dot[col] += dir[col][i] * directions[direction][i];
        }
      if (dot[col] > closestDot[col]) 
        {
        closestDot[col] = dot[col];
        closestAxis[col] = direction;
        }
      }
    }

  std::cout << "Closest axis: " << closestAxis[0] << ", " << closestAxis[1] << ", " << closestAxis[2] << std::endl;
  typename ImageType::DirectionType newDir;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      newDir[i][j] = directions[closestAxis[i]][j];
    }
  }

  image->SetDirection(newDir);

  typename TransformType::Pointer tiltTfm = TransformType::New();
  typename TransformWriterType::Pointer tfmWriter = TransformWriterType::New();
  typename TransformType::MatrixType tiltMatrix, origMatrix, dirMatrix;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dirMatrix[i][j] = newDir[i][j];
      origMatrix[i][j] = dir[i][j];
    }
  }

  tiltMatrix = origMatrix*dirMatrix.GetInverse();
  tiltMatrix = tiltMatrix.GetInverse();
  typename ImageType::PointType origin, origOrigin = image->GetOrigin();
  itk::Point<float,3> origOriginPt, originPt;
  origOriginPt = origOrigin;
  originPt = tiltMatrix*origOriginPt;
  origin = originPt;
  image->SetOrigin(origin);

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput(image);
  writer->Update();  

  tiltTfm->SetMatrix(tiltMatrix);
  tfmWriter->SetInput(tiltTfm);
  tfmWriter->SetFileName(outputTransform.c_str());
  tfmWriter->Update();

  return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
