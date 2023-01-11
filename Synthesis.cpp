//Mathew Bernard
//Image Synthesis
//12/9/22

#include <OpenImageIO/imageio.h>
#include <iostream>
#include <GL/glut.h>

using namespace std;
OIIO_NAMESPACE_USING

struct Pixel{ // defines a pixel structure
	float r,g,b,a;
}; 

//
// Global variables and constants
//

// everything with a "2" preceeding it is referring to the target image

const int DEFAULTWIDTH = 600;	// default window dimensions if no image
const int DEFAULTHEIGHT = 600;

int WinWidth, WinHeight;	// window width and height
int WinWidth2, WinHeight2;
int ImWidth, ImHeight;		// image width and height
int ImChannels;           // number of channels per image pixel

int ImWidth2, ImHeight2;
int ImChannels2;

float avg_red_source, avg_blue_source, avg_green_source;
float avg_red_target, avg_blue_target, avg_green_target;

float stdr, stdg, stdb;
float stdred_target, stdgreen_target, stdblue_target;

float FLT_TRUE_MIN = numeric_limits<float>::denorm_min(); 

float L, M, S; //variables to hold LMS values for source and target
float L2, M2, S2;

int VpWidth, VpHeight;		// viewport width and height
int Xoffset, Yoffset;     // viewport offset from lower left corner of window

Pixel **source = NULL;  // the image source used for OpenGL display
Pixel **target = NULL; // the image source for the second image
int pixformat; 			// the pixel format used to correctly  draw the image
int pixformat2;

string saveimage = "None";

//
//  Routine to cleanup the memory.   
//
void destroy(){
 if (source){
     delete source[0];
	 delete source;  
  }
}


//
//  Routine to read an image file and store in a source
//  returns the size of the image in pixels if correctly read, or 0 if failure
//
int readImage(string infilename){
// Create the oiio file handler for the image, and open the file for reading the image.
  // Once open, the file spec will indicate the width, height and number of channels.
  std::unique_ptr<ImageInput> infile = ImageInput::open(infilename);
  if(!infile){
    cerr << "Could not input image file " << infilename << ", error = " << geterror() << endl;
    return 0;
  }

  // Record image width, height and number of channels in global variables
  ImWidth = infile->spec().width;
  ImHeight = infile->spec().height;
  ImChannels = infile->spec().nchannels;
 
  // allocate temporary structure to read the image 
  vector<float> tmp_pixels(ImWidth * ImHeight * ImChannels);

  // read the image into the tmp_pixels from the input file, flipping it upside down using negative y-stride,
  // since OpenGL sources have the bottom scanline first, and 
  // oiio expects the top scanline first in the image file.
  if(!infile->read_image(TypeDesc::FLOAT, &tmp_pixels[0])){
    cerr << "Could not read image from " << infilename << ", error = " << geterror() << endl;
    return 0;
  }
 // get rid of the old source and make a new one of the new size
  destroy();
 // allocate space for the source (contiguous approach, 2d style access)
  source = new Pixel*[ImHeight];
  if(source != NULL)
	source[0] = new Pixel[ImWidth * ImHeight];
  for(int i = 1; i < ImHeight; i++)
	source[i] = source[i - 1] + ImWidth;

 
 //  assign the read pixels to the the data structure
 int index;
  //Work backwards to flip vertically
  for(int row = ImHeight-1, x=0; row >= 0; row--, x++) {
    //No need to flip horizontally
    for(int col = 0; col < ImWidth; col++) {
        index = (x*ImWidth+col)*ImChannels;
        if (ImChannels==1){ 
          source[row][col].r = tmp_pixels[index];
			    source[row][col].g = tmp_pixels[index];
			    source[row][col].b = tmp_pixels[index];
			    source[row][col].a = 255;
		    }
        else{
		      source[row][col].r = tmp_pixels[index];
			    source[row][col].g = tmp_pixels[index+1];
			    source[row][col].b = tmp_pixels[index+2];			
			  if (ImChannels <4) // no alpha value is present so set it to 255
				  source[row][col].a = 255; 
			  else // read the alpha value
				  source[row][col].a = tmp_pixels[index+3];			
      }
    }
  }
 
  // close the image file after reading, and free up space for the oiio file handler
  infile->close();
  
  // set the pixel format to GL_RGBA and fix the # channels to 4  
  pixformat = GL_RGBA;  
  ImChannels = 4;

  // return image size in pixels
  return ImWidth * ImHeight;
}



//
// Routine to display a source in the current window
//
void displayImage(){
  // if the window is smaller than the image, scale it down, otherwise do not scale
  if(WinWidth < ImWidth2  || WinHeight < ImHeight2)
    glPixelZoom(float(VpWidth) / ImWidth2, float(VpHeight) / ImHeight2);
  else
    glPixelZoom(1.0, 1.0);
  
  // display starting at the lower lefthand corner of the viewport
  glRasterPos2i(0, 0);
  

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glDrawPixels(ImWidth2, ImHeight2, pixformat, GL_FLOAT, target[0]);
}


//
// Routine to write the current framebuffer to an image file
//
void writeImage(string outfilename){
  // make a source that is the size of the window and grab OpenGL framebuffer into it
  // alternatively, you can read the source into a 1d array and export this 
   unsigned char local_source[WinWidth * WinHeight * ImChannels];
   glReadPixels(0, 0, WinWidth, WinHeight, pixformat, GL_UNSIGNED_BYTE, local_source);
  
  // create the oiio file handler for the image
  std::unique_ptr<ImageOutput> outfile = ImageOutput::create(outfilename);
  if(!outfile){
    cerr << "Could not create output image for " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // Open a file for writing the image. The file header will indicate an image of
  // width WinWidth, height WinHeight, and ImChannels channels per pixel.
  // All channels will be of type unsigned char
  ImageSpec spec(WinWidth, WinHeight, ImChannels, TypeDesc::UINT8);
  if(!outfile->open(outfilename, spec)){
    cerr << "Could not open " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // Write the image to the file. All channel values in the source are taken to be
  // unsigned chars. While writing, flip the image upside down by using negative y stride, 
  // since OpenGL sources have the bottom scanline first, and oiio writes the top scanline first in the image file.
  int scanlinesize = WinWidth * ImChannels * sizeof(unsigned char);
  if(!outfile->write_image(TypeDesc::UINT8, local_source + (WinHeight - 1) * scanlinesize, AutoStride, -scanlinesize)){
    cerr << "Could not write image to " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // close the image file after the image is written and free up space for the
  // ooio file handler
  outfile->close();
}

//take the image and copy it to the second source so we
//can reuse the readimage function
void copyImage(){
  
  target = new Pixel*[ImHeight];
  if(target != NULL)
	target[0] = new Pixel[ImWidth * ImHeight];
  for(int i = 1; i < ImHeight; i++)
	target[i] = target[i - 1] + ImWidth;

  for(int row = 0; row < ImHeight; ++row) {
    for(int col = 0; col < ImWidth; ++col) {
    
		target[row][col].r = source[row][col].r;
		target[row][col].g = source[row][col].g;
		target[row][col].b = source[row][col].b;
	}
  }
  ImHeight2 = ImHeight;
  ImWidth2 = ImWidth;
}


//
//   Display Callback Routine: clear the screen and draw the current image
//
void handleDisplay(){
  
  // specify window clear (background) color to be opaque black
  glClearColor(0, 0, 0, 1);
  // clear window to background color
  glClear(GL_COLOR_BUFFER_BIT);  
  
  // only draw the image if it is of a valid size
  if(ImWidth > 0 && ImHeight > 0)
    displayImage();
  
  // flush the OpenGL pipeline to the viewport
  glFlush();
}

//
//RGB TO LMS conversion for the first image (source)
//
Pixel rgb_lms(Pixel rgbPixmap){

  float l, m, s;
  //create a new Pixel array to store the lms pixmap
  Pixel LMSPixel;

  //prevents log(0) from happening
  rgbPixmap.r = max(FLT_TRUE_MIN, rgbPixmap.r);
  rgbPixmap.g = max(FLT_TRUE_MIN, rgbPixmap.g);
  rgbPixmap.b = max(FLT_TRUE_MIN, rgbPixmap.b);
  LMSPixel.a = rgbPixmap.a;
  
  //apply the matrix to convert to lms
  l = 0.3811*rgbPixmap.r + 0.5783*rgbPixmap.g + 0.0402*rgbPixmap.b;
  m = 0.1967*rgbPixmap.r + 0.7244*rgbPixmap.g + 0.0782*rgbPixmap.b;
  s = 0.0241*rgbPixmap.r + 0.1288*rgbPixmap.g + 0.8444*rgbPixmap.b;

  //convert to logrithmic space
  LMSPixel.r = log10(l);
  LMSPixel.g = log10(m);
  LMSPixel.b = log10(s);

  return LMSPixel;

}

//
// LMS TO LAB conversion
//
Pixel lms_lab(Pixel LMSPixel){

  //variables
  float l, a, b;
  Pixel LABPixel;

    //apply 1st matrix
  l = 1.0*LMSPixel.r + 1.0*LMSPixel.g + 1.0*LMSPixel.b;
  a = 1.0*LMSPixel.r + 1.0*LMSPixel.g + -2.0*LMSPixel.b;
  b = 1.0*LMSPixel.r + -1.0*LMSPixel.g + 0.0*LMSPixel.b;
  LABPixel.a = LMSPixel.a;
  //apply 2nd matrix
  LABPixel.r = l / pow(3.f, 0.5);
  LABPixel.g = a / pow(6.f, 0.5);
  LABPixel.b = b / pow(2.f, 0.5);
      
  return LABPixel;
}

//
//conversion of image from lab to lms colorspace
//
Pixel lab_lms(Pixel LABPixel){

  //create variables
  float l, m, s;
  Pixel LMSPixel;

 //apply 1st matrix
  l = LABPixel.r * (pow(3.f, 0.5)/3.f);
  m = LABPixel.g * (pow(6.f, 0.5)/6.f);
  s = LABPixel.b * (pow(2.f, 0.5)/2.f);

  //apply 2nd matrix
  LMSPixel.r = 1.0*l + 1.0*m + 1.0*s;
  LMSPixel.g = 1.0*l + 1.0*m + -1.0*s;
  LMSPixel.b = 1.0*l + -2.0*m + 0.0*s;

  LMSPixel.a = LABPixel.a;

  return LMSPixel;
}

//
//Conversion of image from lms to rgb colorspace
//
Pixel lms_rgb(Pixel LMSPixel){

  float r, g, b;
  Pixel RGBPixel;

  //reduce from log space  
  r = pow(10, LMSPixel.r);
  g = pow(10, LMSPixel.g);
  b = pow(10, LMSPixel.b);

  //apply matrix
  RGBPixel.r = 4.4679*r + -3.5873*g + 0.1193*b;
  RGBPixel.g = -1.2186*r + 2.3809*g + -0.1624*b;
  RGBPixel.b = 0.0497*r + -0.2439*g + 1.2045*b;

  //bound check to sure final RGB value is between 0 and 1
  RGBPixel.r = max(0.f, min(1.f, RGBPixel.r));
  RGBPixel.g = max(0.f, min(1.f, RGBPixel.g));
  RGBPixel.b = max(0.f, min(1.f, RGBPixel.b));
  RGBPixel.a = LMSPixel.a;

  return RGBPixel;

}

void synthesis(){

    //
    // Stddev calculations for target image
    //
    double lAvg1 = 0, aAvg1 = 0, bAvg1 = 0;
    double lStddev1 = 0, aStddev1 = 0, bStddev1 = 0;

    //converstion of source image from rgb to lab
    for(int i=0; i<ImHeight; i++){
      for(int j=0; j<ImWidth; j++){
        source[i][j] = rgb_lms(source[i][j]);
        source[i][j] = lms_lab(source[i][j]);

        //calculate the sums for each channel
        lAvg1 += source[i][j].r;
        aAvg1 += source[i][j].g;
        bAvg1 += source[i][j].b;

      }
    }

    cout << lAvg1 << endl;
    cout << "check" << endl;

    //calculate the averages for stddev
    lAvg1 /= (ImHeight*ImWidth);
    aAvg1 /= (ImHeight*ImWidth);
    bAvg1 /= (ImHeight*ImWidth);

    //calculate the stddev for LAB channels
    for(int i=0; i<ImHeight; i++){
        for(int j=0; j<ImWidth; j++){
            lStddev1 += pow(source[i][j].r - lAvg1, 2);
            aStddev1 += pow(source[i][j].g - aAvg1, 2);
            bStddev1 += pow(source[i][j].b - bAvg1, 2);
        }
    }
    lStddev1 = pow(lStddev1 / (ImHeight*ImWidth), 0.5);
    aStddev1 = pow(aStddev1 / (ImHeight*ImWidth), 0.5);
    bStddev1 = pow(bStddev1 / (ImHeight*ImWidth), 0.5);
    
    cout << "check2" << endl;
    //
    // Stddev calculations for source image
    //
    double lAvg2 = 0, aAvg2 = 0, bAvg2 = 0;
    double lStddev2 = 0, aStddev2 = 0, bStddev2 = 0;

    //convert RGBtoLAB
    for(int i=0; i<ImHeight2; i++){
        for(int j=0; j<ImWidth2; j++){
            target[i][j] = rgb_lms(target[i][j]);
            target[i][j] = lms_lab(target[i][j]);

            //calc sums
            lAvg2 += target[i][j].r;
            aAvg2 += target[i][j].g;
            bAvg2 += target[i][j].b;
        }
    }

cout << lAvg2 << endl;
cout << "check3" << endl;

    //calculate averages for stddev
    lAvg2 /= (ImHeight2*ImWidth2);
    aAvg2 /= (ImHeight2*ImWidth2);
    bAvg2 /= (ImHeight2*ImWidth2);

    //calculate stddev for LAB channels
    for(int i=0; i<ImHeight2; i++){
        for(int j=0; j< ImWidth2; j++){
            lStddev2 += pow(target[i][j].r - lAvg2, 2);
            aStddev2 += pow(target[i][j].g - aAvg2, 2);
            bStddev2 += pow(target[i][j].b - bAvg2, 2);
        }
    }
    lStddev2 = pow(lStddev2 / (ImHeight2*ImWidth2), 0.5);
    aStddev2 = pow(aStddev2 / (ImHeight2*ImWidth2), 0.5);
    bStddev2 = pow(bStddev2 / (ImHeight2*ImWidth2), 0.5);
    
    cout << "check4" << endl;
    //
    // perform transfer
    //
    for(int i=0; i<ImHeight; i++){
        for(int j=0; j<ImWidth; j++){
            //subtract target averages
            source[i][j].r -= lAvg1;   
            source[i][j].g -= aAvg1;
            source[i][j].b -= bAvg1;

            //multipy by target/source stddev
            source[i][j].r *= (lStddev1/lStddev2);
            source[i][j].g *= (aStddev1/aStddev2);
            source[i][j].b *= (bStddev1/bStddev2);

            //add source averages
            source[i][j].r += lAvg2;
            source[i][j].g += aAvg2;
            source[i][j].b += bAvg2;
        }
    }

  cout << "check5" << endl;
    //
    // convert LAB to RGB
    //
    for(int i=0; i<ImHeight; i++){
        for(int j=0; j<ImWidth; j++){
            source[i][j] = lab_lms(source[i][j]);
            source[i][j] = lms_rgb(source[i][j]);
        }
    }
  cout << "check6" << endl;
  return;

}

//
//  Keyboard Callback Routine: 'r' - read and display a new image,
//  'w' - write the current window to an image file, 'q' or ESC - quit
//
void handleKey(unsigned char key, int x, int y){
  string infilename, outfilename;
  int ok;
  
  switch(key){
    case 'w':		// 'w' - write the image to a file
    case 'W':
      if(saveimage == "None"){
      cout << "Output image filename? ";  // prompt user for output filename
      cin >> outfilename;
      writeImage(outfilename);
      }
      else{
        writeImage(saveimage);
      }
      break;
    
    // case 'l': // 'l' - convert RGB space to LMS
    // case 'L':
    //   rgb_lms_image1(); //for the source image
    //   rgb_lms_image2(); //for the target image
    //   lms_lab(); //handles the conversion for both images
    //   standard_deviation();
    //   lab_lms(); //conversion of the target image
    //   lms_rgb(); //conversion of the target image
    //   glutPostRedisplay();

    //   break;

    case 's': //'s' - Synthesis main code
    case 'S':
      synthesis();
      glutPostRedisplay();
      break;
	
	case 'q':		// q or ESC - quit
    case 'Q':
    case 27:
      destroy();
      exit(0);
      
    default:		// not a valid key -- just ignore it
      return;
  }
}


//
//  Reshape Callback Routine: If the window is too small to fit the image,
//  make a viewport of the maximum size that maintains the image proportions.
//  Otherwise, size the viewport to match the image size. In either case, the
//  viewport is centered in the window.
//
void handleReshape(int w, int h){

  float width = min(ImWidth, ImWidth2);
  float height = min(ImHeight, ImHeight2);

  float imageaspect = (float)width / (float)height;	// aspect ratio of image
  float newaspect = (float)w / (float)h; // new aspect ratio of window
  
  // record the new window size in global variables for easy access
  width = w;
  height = h;
  
  // if the image fits in the window, viewport is the same size as the image
  if(w >= width && h >= height){
    Xoffset = (w - width) / 2;
    Yoffset = (h - height) / 2;
    VpWidth = width;
    VpHeight = height;
  }
  // if the window is wider than the image, use the full window height
  // and size the width to match the image aspect ratio
  else if(newaspect > imageaspect){
    VpHeight = h;
    VpWidth = int(imageaspect * VpHeight);
    Xoffset = int((w - VpWidth) / 2);
    Yoffset = 0;
  }
  // if the window is narrower than the image, use the full window width
  // and size the height to match the image aspect ratio
  else{
    VpWidth = w;
    VpHeight = int(VpWidth / imageaspect);
    Yoffset = int((h - VpHeight) / 2);
    Xoffset = 0;
  }
  
  // center the viewport in the window
  glViewport(Xoffset, Yoffset, VpWidth, VpHeight);
  
  // viewport coordinates are simply pixel coordinates
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, VpWidth, 0, VpHeight);
  glMatrixMode(GL_MODELVIEW);
}


//
// Main program to scan the commandline, set up GLUT and OpenGL, and start Main Loop
//
int main(int argc, char* argv[]){
  // scan command line and process
  // only one parameter allowed, an optional image filename and extension
  if(argc > 4){
    cout << "usage: ./final [inputimage1.ext][inputimage2.ext][saveimage.ext]" << endl;
    exit(1);
  }
  
  // set up the default window and empty source if no image or image fails to load
  WinWidth = DEFAULTWIDTH;
  WinHeight = DEFAULTHEIGHT;
  ImWidth = 0;
  ImHeight = 0;
  
  WinWidth2 = DEFAULTWIDTH;
  WinHeight2 = DEFAULTHEIGHT;
  ImWidth2 = 0;
  ImHeight2 = 0;

  // load the image if present, and size the window to match
  if(argc == 3 ){
    if(readImage(argv[1])){
      WinWidth = ImWidth;
      WinHeight = ImHeight;
    }
    copyImage();
    if(readImage(argv[2])){
      WinWidth = ImWidth;
      WinHeight = ImHeight;
    }
    //saveimage = argv[3];
  }

  //load the image if present, and size the window to matc
  //as well as create the name to save the image out to when finished
  if(argc == 4 ){
    if(readImage(argv[1])){
      WinWidth = ImWidth;
      WinHeight = ImHeight;
    }
    copyImage();
    if(readImage(argv[2])){
      WinWidth = ImWidth;
      WinHeight = ImHeight;
    }
    saveimage = argv[3];
  }
  
  //reshape the output to the smaller image when calling glutInit
  float width = min(ImWidth, ImWidth2);
  float height = min(ImHeight, ImHeight2);

  // start up GLUT
  glutInit(&argc, argv);
  
  // create the graphics window, giving width, height, and title text
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
  glutInitWindowSize(width, height);
  glutCreateWindow("Final Project");
  
  // set up the callback routines
  glutDisplayFunc(handleDisplay); // display update callback
  glutKeyboardFunc(handleKey);	  // keyboard key press callback
  glutReshapeFunc(handleReshape); // window resize callback
  
  
  // Enter GLUT's event loop
  glutMainLoop();
  return 0;
}
