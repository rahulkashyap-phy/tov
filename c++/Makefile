CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -Wextra
TARGET = tov
TARGET_FIXED = tov_fixed
SOURCE = tov.cpp
SOURCE_FIXED = tov_fixed.cpp

all: $(TARGET) $(TARGET_FIXED)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

$(TARGET_FIXED): $(SOURCE_FIXED)
	$(CXX) $(CXXFLAGS) -o $(TARGET_FIXED) $(SOURCE_FIXED)

clean:
	rm -f $(TARGET) $(TARGET_FIXED)

.PHONY: clean all

