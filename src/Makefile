CXX = g++
CXXFLAGS=-std=gnu++11 -W -Wall -Wextra -Wno-address-of-packed-member -Wno-strict-aliasing -Wno-array-bounds -Wno-unknown-warning-option -Wno-format -Wno-unknown-pragmas -Wno-unused-private-field -g -O3 -DDEBUG

../mspack: main.o encoder.o decoder.o msprint.o compare.o tinyxml2.o SHA1.o
	$(CXX) $(CXXFLAGS) main.o encoder.o decoder.o msprint.o compare.o tinyxml2.o SHA1.o -o ../mspack

main.o: main.cpp mspack.h base64.h utils.h head.h lossy.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

encoder.o: encoder.cpp mspack.h utils.h head.h encoder.h lossy.h
	$(CXX) $(CXXFLAGS) -c encoder.cpp -o encoder.o

decoder.o: decoder.cpp mspack.h utils.h head.h decoder.h lossy.h
	$(CXX) $(CXXFLAGS) -c decoder.cpp -o decoder.o

tinyxml2.o: tinyxml2.cpp tinyxml2.h
	$(CXX) $(CXXFLAGS) -c tinyxml2.cpp -o tinyxml2.o

msprint.o: msprint.cpp head.h
	$(CXX) $(CXXFLAGS) -c msprint.cpp -o msprint.o

compare.o: compare.cpp head.h
	$(CXX) $(CXXFLAGS) -c compare.cpp -o compare.o

SHA1.o: SHA1.cpp SHA1.h
	$(CXX) $(CXXFLAGS) -c SHA1.cpp -o SHA1.o

clean:
	rm *.o ../mspack
