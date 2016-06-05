CXX= mpicxx
CXXFLAGS= -std=c++11 -fopenmp 
#-DDCO_AMPI -I dco_cpp -I AdjointMPI_example/AdjointMPI/include/ -I AdjointMPI_example/AdjointMPI/interfaces/dco -L dco_cpp -ldco AdjointMPI_example/AdjointMPI/libAMPI.a
OBJ=main.o matrix_template.o matrix_vector_product.o functions.o
EXE=Test

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf $(OBJ) $(EXE) main.o
