# Image-Synthesis
This program was created to apply the contents of the “Color Transfer between Images” by Reinhard et al.This is a tool for image synthesis, where one image can apply its color and characteristics to another image. 

This program uses the OpenImageIO API to read and write images and uses the 	OpenGL/GLUT to display them to the screen.The user is able to read in two images, a source image and a target image. The program gathers information about the source image and then applies it to the target image. The user should run the program by the following command:

./synthesis target-image.png source-image.png [save-image.png]

Where [save-image.png] is an optional parameter for a file name to save the program out to.

Within the program the user can operate different features by:

W - this will write the image out to a file if specified in the third parameter, otherwise the program will ask the user for the name of a file to save to.

L - this is where the image synthesis happens, first converting the images to LMS and then lαβ colorspace. Transfer the source images lαβ colorspace onto target images lαβ colorspace. Convert the target back to LMS and then to RGB to display the new image.

Q - used to exit the program
