CXX = g++
CXXFLAGS=-std=gnu++11 -W -Wall -Wextra -Wno-strict-aliasing -Wno-address-of-packed-member -Wno-format -g -O3 -DDEBUG

../mspredict: main.o encoder.o decoder.o msprint.o compare.o tinyxml2.o
	$(CXX) $(CXXFLAGS) main.o encoder.o decoder.o msprint.o compare.o tinyxml2.o -o ../mspredict

main.o: main.cpp mspredict.h base64.h utils.h head.h lossy.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

encoder.o: encoder.cpp mspredict.h utils.h head.h encoder.h lossy.h
	$(CXX) $(CXXFLAGS) -c encoder.cpp -o encoder.o

decoder.o: decoder.cpp mspredict.h utils.h head.h decoder.h lossy.h
	$(CXX) $(CXXFLAGS) -c decoder.cpp -o decoder.o

tinyxml2.o: tinyxml2.cpp tinyxml2.h
	$(CXX) $(CXXFLAGS) -c tinyxml2.cpp -o tinyxml2.o

msprint.o: msprint.cpp head.h
	$(CXX) $(CXXFLAGS) -c msprint.cpp -o msprint.o

compare.o: compare.cpp head.h
	$(CXX) $(CXXFLAGS) -c compare.cpp -o compare.o

clean:
	rm *.o ../mspredict
