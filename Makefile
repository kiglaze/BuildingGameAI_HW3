# Makefile was generated with ChatGPT, based on the compile commands listed in the README.md file.
# Compiler settings - Can be customized.
CXX = g++ 
CXXFLAGS = -Wall -g

# Linker settings
LDFLAGS = -lsfml-graphics -lsfml-window -lsfml-system

# Add all cpp files to your project here
SOURCES = main.cpp Sprite.cpp SpriteCollection.cpp SteeringBehavior.cpp Arrive.cpp SteeringData.cpp Align.cpp Face.cpp LookWhereYoureGoing.cpp Kinematic.cpp Crumb.cpp

# Object files
OBJS = $(SOURCES:.cpp=.o)

# Executable name
EXEC = sfml-app

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $<

clean:
	rm -f $(OBJS) $(EXEC)
