convert image*.ppm image.gif
mogrify -format png *.ppm
mv image*.png image_png
mv image*.ppm image_ppm
