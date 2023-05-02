BettiMatching: BettiMatching.o cube.o cubical_grid_complex.o write_pair.o union_find.o boundary_enumerator.o top_dimension.o dimension_1.o dimension_0.o match.o
	g++ BettiMatching.o cube.o cubical_grid_complex.o write_pair.o union_find.o boundary_enumerator.o top_dimension.o  dimension_1.o dimension_0.o match.o -o BettiMatching

BettiMatching.o: BettiMatching.cpp top_dimension.h dimension_1.h dimension_0.h match.h
	g++ --std=c++17 -c BettiMatching.cpp

boundary_enumerator.o: boundary_enumerator.cpp boundary_enumerator.h
	g++ --std=c++17 -c boundary_enumerator.cpp

cube.o: cube.cpp cube.h
	g++ --std=c++17 -c cube.cpp

cubical_grid_complex.o: cubical_grid_complex.cpp cubical_grid_complex.h
	g++ --std=c++17 -c cubical_grid_complex.cpp

dimension_0.o: dimension_0.cpp dimension_0.h union_find.h
	g++ --std=c++17 -c dimension_0.cpp

dimension_1.o: dimension_1.cpp dimension_1.h boundary_enumerator.h
	g++ --std=c++17 -c dimension_1.cpp

match.o: match.cpp match.h
	g++ --std=c++17 -c match.cpp

top_dimension.o: top_dimension.cpp top_dimension.h union_find.h
	g++ --std=c++17 -c top_dimension.cpp

union_find.o: union_find.cpp union_find.h
	g++ --std=c++17 -c union_find.cpp

write_pair.o: write_pair.cpp write_pair.h
	g++ --std=c++17 -c write_pair.cpp

clean:
	rm -rf *o BettiMatching

