all: c_lib cpp

c_lib:
	gcc -c soft20_wrap.c -I. -I./soft20/include -I/usr/include -L./soft20 -lsoft1 -lfftw3 -lm
	ar -crs libsoft20_api.a soft20_wrap.o

cpp:
	g++ main.cpp -std=c++17 -L./ -L./soft20 -lsoft20_api -lsoft1 -lfftw3 -lm -o main
