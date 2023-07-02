BettiMatching: BettiMatching.o inter_dimensions.o top_dimension.o enumerators.o data_structures.o utils.o
	g++ --std=c++17 BettiMatching.o inter_dimensions.o top_dimension.o enumerators.o data_structures.o utils.o -o BettiMatching

BettiMatching.o: BettiMatching.cpp inter_dimensions.h top_dimension.h utils.h config.h
	g++ --std=c++17 -c BettiMatching.cpp

inter_dimensions.o: inter_dimensions.cpp inter_dimensions.h enumerators.h
	g++ --std=c++17 -c inter_dimensions.cpp

top_dimension.o: top_dimension.cpp top_dimension.h enumerators.h
	g++ --std=c++17 -c top_dimension.cpp

enumerators.o: enumerators.cpp enumerators.h
	g++ --std=c++17 -c enumerators.cpp

data_structures.o: data_structures.cpp data_structures.h template_functions.h
	g++ --std=c++17 -c data_structures.cpp

utils.o: utils.cpp utils.h
	g++ --std=c++17 -c utils.cpp

clean:
	rm -rf *o BettiMatching

