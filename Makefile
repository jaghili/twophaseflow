# DIRECTORIES
MKLROOT := /home/${USER}/intel
RLSDIR := .
SRC_DIR := ./src
OBJ_DIR := ./obj

# FILES
HEADER_FILES := $(wildcard $(SRC_DIR)/*.h)
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

# COMPILERS 
#CC := ${MKLROOT}/bin/icc -qopenmp
CC := g++ -fopenmp

# LIBRARIES
LDFLAGS := -L./lib/
#LDFLAGS += -L${MKLROOT}/mkl/lib/intel64_lin/ # Intel's MKL library

# Using SuiteSparse's UMFPACK 
#LDFLAGS += -L/home/${USER}/codes/SuiteSparse/UMFPACK/Lib/
#LDFLAGS += -L/home/${USER}/codes/SuiteSparse/AMD/Lib/
#LDFLAGS += -L/home/${USER}/codes/SuiteSparse/SuiteSparse_config/
#LDFLAGS += -L/home/${USER}/codes/SuiteSparse/CHOLMOD/Lib/
#LDFLAGS += -L/home/${USER}/codes/SuiteSparse/COLAMD/Lib/

LDFLAGS += -lMeshReader # FVCA5 Mesh reader library
LDFLAGS += -lsfml-system -lsfml-window -lsfml-graphics
#LDFLAGS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LDFLAGS += -lpthread -lm -ldl
#LDFLAGS += -lsuperlu
#LDFLAGS += -lumfpack -lamd -lsuitesparseconfig -lcholmod -lcolamd

INCLUDES := -I/usr/include/eigen3/
#INCLUDES += -I/usr/include/eigen3/unsupported/   
INCLUDES += -I./src/
#INCLUDES += -I/home/${USER}/codes/superlu/SRC/
#INCLUDES += -I/home/${USER}/intel/compilers_and_libraries_2018.1.163/linux/mkl/include/

#INCLUDES += -I/home/${USER}/codes/SuiteSparse/UMFPACK/Include/
#INCLUDES += -I/home/${USER}/codes/SuiteSparse/SuiteSparse_config/
#INCLUDES += -I/home/${USER}/codes/SuiteSparse/AMD/Include/
#INCLUDES += -I/home/${USER}/codes/SuiteSparse/CHOLMOD/Include/
#INCLUDES += -I/home/${USER}/codes/SuiteSparse/COLAMD/Include/

CXXFLAGS := -Ofast  $(DEBUG)
#CXXFLAGS += -DEIGEN_USE_MKL_ALL -DMKL_LP64
#CXXFLAGS += -DRECORDMOVIE 			#Uncomment to save screenshots
#CXXFLAGS += -DBF_DEBUG

OUTDIR := ./output

all: $(RLSDIR)/frac.exe

$(RLSDIR)/frac.exe: obj/CoreyConstant.o obj/Column.o obj/Data.o obj/Corey.o obj/Constants.o obj/Physics.o obj/Corey.o obj/Flux.o obj/Jacobian.o obj/Residual.o obj/frac.o
	$(CC)  $^ $(LDFLAGS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADER_FILES)
	$(CC) $(INCLUDES) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR)/*.o
	rm -rf $(RLSDIR)/*.exe
	rm -rf $(OUTDIR)/movie/*.jpg
	rm -rf $(OUTDIR)/*.dat
	rm -rf $(OUTDIR)/*.tmp

./output/movie/constant.avi: ./output/movie/constant*.jpg
	ffmpeg -start_number 0 -i constant_%d.jpg -vcodec mjpeg ./output/movie/constant.avi
