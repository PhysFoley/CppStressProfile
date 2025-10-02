stresscalc: main.cpp mat3.cpp
	g++ -pthread -O3 -o stresscalc main.cpp mat3.cpp -Wall -Wpedantic

clean:
	rm -f stresscalc

.PHONY: clean
