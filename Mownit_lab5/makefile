all:
	gcc -o main main.c `gsl-config --cflags --libs`
test:
	./main
clean:
	rm -rf main C_results.csv