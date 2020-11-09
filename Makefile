#para ver las advertencias
#--cxxflags
#BINNAME?=XXXXXXXXXXXXXXX

CC=g++ -std=c++11 
CFLAGS=`Magick++-config --cppflags --ldflags --libs` `pkg-config gtkmm-3.0 --cflags --libs`  -Wl,-V
CFLAGS2=`pkg-config gtkmm-3.0 --cflags --libs` -Wl,-V 
CFLAGS3=`-O2 -L/usr/X11R6/lib -lm -lpthread -lX11` 


IMG = mglp
$(IMG):	
	$(CC) ./src/main.cpp -o ./dist/$@ -O2 -L/usr/X11R6/lib -lm -lpthread -lX11

clean:
	rm -f $(IMG)
