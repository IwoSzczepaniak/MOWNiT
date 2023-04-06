all:
	gcc -O0 main.c `gsl-config --cflags --libs` -o main
test:
	./main
clean:
	rm -rf main C_results.csv