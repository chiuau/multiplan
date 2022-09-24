all: multiplan

multiplan:
	g++ -std=c++2a -O3 *.cpp -o multiplan
