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
//RGB TO LMS conversion for each pixel, returns the pixel in LMS form
//
Pixel rbg_to_lms(Pixel rgbPixmap){

  float r, g, b;
  float l, m, s;
  Pixel lmsPixmap;


  r = max(FLT_TRUE_MIN, float(0.3811*rgbPixmap.r + 0.5783*rgbPixmap.g + 0.0402*rgbPixmap.b));
  g = max(FLT_TRUE_MIN, float(0.1967*rgbPixmap.r + 0.7244*rgbPixmap.g + 0.0782*rgbPixmap.b));
  b = max(FLT_TRUE_MIN, float(0.0241*rgbPixmap.r + 0.1288*rgbPixmap.g + 0.8444*rgbPixmap.b));
  
  if(r == 0.0)
    r = 1.0;
  if(g == 0.0)
    g = 1.0;
  if(b == 0.0)
    b = 1.0;

  l = log10(r);
  m = log10(g);
  s = log10(b);

  lmsPixmap.r = l;
  lmsPixmap.g = m;
  lmsPixmap.b = s;
      
  return lmsPixmap;

}


//
// LMS TO LAB conversion
//
void lms_lab(){

  //variables for source
  float r, g, b;
  float l, a, beta;
  //variables for target
  float red_target, green_target, blue_target;
  float l2, a2, beta2;

  //
  //setting all the average values to 0
  //
  avg_red_source = 0;
  avg_blue_source = 0;
  avg_green_source = 0;

  avg_red_target = 0;
  avg_blue_target = 0;
  avg_green_target = 0;


  for(int row = 0; row < ImHeight; ++row) {
      for(int col = 0; col < ImWidth; ++col) {
        
        r = 1.0*source[row][col].r + 1.0*source[row][col].g + 1.0*source[row][col].b;
        g = 1.0*source[row][col].r + 1.0*source[row][col].g + -2.0*source[row][col].b;
        b = 1.0*source[row][col].r + -1.0*source[row][col].g + 0.0*source[row][col].b;

        l = (1.0/sqrt(3.0))*r;
        a = (1.0/sqrt(6.0))*g;
        beta =(1.0/sqrt(2.0))*b;

        source[row][col].r = l;
        source[row][col].g = a;
        source[row][col].b = beta;

        avg_red_source += source[row][col].r;
        avg_green_source += source[row][col].g;
        avg_blue_source += source[row][col].b;
      }
  }

  avg_red_source = avg_red_source/(ImHeight*ImWidth);
  avg_green_source = avg_green_source/(ImHeight*ImWidth);
  avg_blue_source = avg_blue_source/(ImHeight*ImWidth);

    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
      
        red_target = 1.0*target[row][col].r + 1.0*target[row][col].g + 1.0*target[row][col].b;
        green_target = 1.0*target[row][col].r + 1.0*target[row][col].g + -2.0*target[row][col].b;
        blue_target = 1.0*target[row][col].r + -1.0*target[row][col].g + 0.0*target[row][col].b;

        l2 = (1.0/sqrt(3.0))*red_target;
        a2 = (1.0/sqrt(6.0))*green_target;
        beta2 = (1.0/sqrt(2.0))*blue_target;

        target[row][col].r = l2;
        target[row][col].g = a2;
        target[row][col].b = beta2;

        avg_red_target += target[row][col].r;
        avg_green_target += target[row][col].g;
        avg_blue_target += target[row][col].b;

   

      }
  }


  avg_red_target = avg_red_target/(ImHeight2*ImWidth2);
  avg_green_target = avg_green_target/(ImHeight2*ImWidth2);
  avg_blue_target = avg_blue_target/(ImHeight2*ImWidth2);


}

//
//Calculation of the standard deviation
//
void standard_deviation(){

float sumr = 0.0, sumg = 0.0, sumb = 0.0; 
float mean_red_of_source, mean_green_of_source, mean_blue_of_source;
float std_red_source = 0.0, std_green_source = 0.0, std_blue_source = 0.0;
float sumred_target = 0.0, sumgreen_target = 0.0, sumblue_target = 0.0;
float mean_red_of_target, mean_green_of_target, mean_blue_of_target;
float std_red_target = 0.0, std_green_target = 0.0, std_blue_target = 0.0;


    for(int row = 0; row < ImHeight; ++row) {
      for(int col = 0; col < ImWidth; ++col) {
        sumr += source[row][col].r;
        sumg += source[row][col].g;
        sumb += source[row][col].b;
      }
    }

  //mean for the second image
  mean_red_of_source = sumr / (ImHeight * ImWidth);
  mean_green_of_source = sumg / (ImHeight * ImWidth);
  mean_blue_of_source = sumb / (ImHeight * ImWidth);

    //subtract the mean from the source image
    for(int row = 0; row < ImHeight; ++row) {
      for(int col = 0; col < ImWidth; ++col) {
        std_red_source += pow(source[row][col].r - mean_red_of_source, 2);
        std_green_source += pow(source[row][col].g - mean_green_of_source, 2);
        std_blue_source += pow(source[row][col].b - mean_blue_of_source, 2);
      }
    }

  std_red_source = sqrt(std_red_source/(ImHeight*ImWidth));
  std_green_source = sqrt(std_green_source/(ImHeight*ImWidth));
  std_blue_source = sqrt(std_blue_source/(ImHeight*ImWidth));

    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
        sumred_target += target[row][col].r;
        sumgreen_target += target[row][col].g;
        sumblue_target += target[row][col].b;
      }
    }


  //mean for the first image loaded
  mean_red_of_target = sumred_target / (ImHeight2 * ImWidth2);
  mean_green_of_target = sumgreen_target / (ImHeight2 * ImWidth2);
  mean_blue_of_target = sumblue_target / (ImHeight2 * ImWidth2);

    //subtract the mean from the target image
    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
        std_red_target += pow(target[row][col].r - mean_red_of_target, 2);
        std_green_target += pow(target[row][col].g - mean_green_of_target, 2);
        std_blue_target += pow(target[row][col].b - mean_blue_of_target, 2);
      }
    }
  
  std_red_target = sqrt(std_red_target/(ImHeight2*ImWidth2));
  std_green_target = sqrt(std_green_target/(ImHeight2*ImWidth2));
  std_blue_target = sqrt(std_blue_target/(ImHeight2*ImWidth2));

    //subtract the average of image 2 from the target image
    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
        target[row][col].r = target[row][col].r - avg_red_target;
        target[row][col].g = target[row][col].g - avg_green_target;
        target[row][col].b = target[row][col].b - avg_blue_target;

      }
    }

    //multiply by the standard deviation ratio of source / target
    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
        target[row][col].r = (std_red_source/std_red_target) * target[row][col].r;
        target[row][col].g = (std_green_source/std_green_target) * target[row][col].g;
        target[row][col].b = (std_blue_source/std_blue_target) * target[row][col].b;
      }
    }

    //add the average of image 1 to the target image
    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
        target[row][col].r = target[row][col].r + avg_red_source;
        target[row][col].g = target[row][col].g + avg_green_source;
        target[row][col].b = target[row][col].b + avg_blue_source;

      }
    }

}

//
//conversion of the target image from lab to lms colorspace
//
void lab_lms(){

  float r = 0.0, g = 0.0, b = 0.0;

    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {

        r = (sqrt(3)/3)*target[row][col].r;
        g = (sqrt(6)/6)*target[row][col].g;
        b = (sqrt(2)/2)*target[row][col].b;

        L2 = 1.0*r + 1.0*g + 1.0*b;
        M2 = 1.0*r + 1.0*g + -1.0*b;
        S2 = 1.0*r + -2.0*g + 0.0*b;

        target[row][col].r = L2;
        target[row][col].g = M2;
        target[row][col].b = S2;

      }
  }
}

//
//Conversion of the target image from lms to rgb colorspace
//
void lms_rgb(){

  float r = 0.0, g = 0.0, b = 0.0;

    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
      

        target[row][col].r = pow(10, target[row][col].r);
        target[row][col].g = pow(10, target[row][col].g);
        target[row][col].b = pow(10, target[row][col].b);

        r = 4.4679*target[row][col].r + -3.5873*target[row][col].g + 0.1193*target[row][col].b;
        g = -1.2186*target[row][col].r + 2.3809*target[row][col].g + -0.1624*target[row][col].b;
        b = 0.0497*target[row][col].r + -0.2439*target[row][col].g + 1.2045*target[row][col].b;

        target[row][col].r = r;
        target[row][col].g = g;
        target[row][col].b = b;

      }
    }
}

void synthesis(){

  for(int row = 0; row < ImHeight; ++row) {
    for(int col = 0; col < ImWidth; ++col) {
        source[row][col] = rbg_to_lms(source[row][col]);
    }
  }
  for(int row = 0; row < ImHeight2; ++row) {
    for(int col = 0; col < ImWidth2; ++col) {
        target[row][col] = rbg_to_lms(target[row][col]);
    }
  }

  lms_lab(); //handles the conversion for both images
  standard_deviation();
  lab_lms(); //conversion of the target image
  lms_rgb(); //conversion of the target image


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
    
    case 'l': // 'l' - convert RGB space to LMS
    case 'L':
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
  if(argc == 1 || argc == 2 || argc > 4){
    cout << "usage: ./synthesis [target.png][source.png][saveimage.png]" << endl;
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
  glutCreateWindow("Image Synthesis Project");
  
  // set up the callback routines
  glutDisplayFunc(handleDisplay); // display update callback
  glutKeyboardFunc(handleKey);	  // keyboard key press callback
  glutReshapeFunc(handleReshape); // window resize callback
  
  
  // Enter GLUT's event loop
  glutMainLoop();
  return 0;
}
