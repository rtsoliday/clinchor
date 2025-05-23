SDDS_REPO = $(firstword $(wildcard ../../SDDS ../../../../../epics/extensions/src/SDDS))
ifeq ($(SDDS_REPO),)
  $(error SDDS source code not found. Run 'git clone https://github.com/rtsoliday/SDDS.git' next to the clinchor repository)
endif

include ../Makefile.rules

PROD = clinchor
clinchor_SRC = clinchor.c
PROD_DEPS += $(OBJ_DIR)/bunchTrain.h \
			$(OBJ_DIR)/clinchor.h \
			$(OBJ_DIR)/generalBunchPattern.h \
			$(OBJ_DIR)/symmetricBunchPattern.h \
			$(OBJ_DIR)/bunchPatternFromFile.h

CFLAGS += -I$(SDDS_REPO)/include -I$(OBJ_DIR) -DUSE_GSL

ifeq ($(OS), Linux)
  CFLAGS += -fopenmp $(MKL_CFLAG) $(MKL_INCLUDE) $(LAPACK_CFLAG) $(LAPACK_INCLUDE)
  LDFLAGS := -fopenmp -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) $(LDFLAGS)
  PROD_SYS_LIBS := $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS) $(MKL_LIB) $(LAPACK_LIB)
  PROD_LIBS_SDDS = -lmdbcommon -lmatlib -lSDDS1 -lnamelist -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Darwin)
  CFLAGS += -DLAPACK
  LDFLAGS := -L$(SDDS_REPO)/lib/$(OS)-$(ARCH) $(LDFLAGS)
  PROD_SYS_LIBS := -framework Accelerate $(LZMA_LIB) $(GSL_LIB) $(GSLCBLAS_LIB) $(Z_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS_SDDS = -lmdbcommon -lmatlib -lSDDS1 -lnamelist -lrpnlib -lmdbmth -lmdblib
endif

ifeq ($(OS), Windows)
  CFLAGS += -I$(SDDS_REPO)/lzma -wd4101 -wd4244 -openmp $(MKL_CFLAG) $(MKL_INCLUDE)
  PROD_SYS_LIBS := gsl.lib gslcblas.lib $(MKL_LIB) $(PROD_SYS_LIBS)
  PROD_LIBS_SDDS = mdbcommon.lib matlib.lib SDDS1.lib namelist.lib rpnlib.lib mdbmth.lib mdblib.lib lzma.lib z.lib
  LIB_LINK_DIRS += -LIBPATH:$(SDDS_REPO)/lib/$(OS)-$(ARCH)
endif

include ../Makefile.build

$(OBJ_DIR)/clinchor.$(OBJEXT): $(OBJ_DIR)/bunchTrain.h $(OBJ_DIR)/clinchor.h $(OBJ_DIR)/generalBunchPattern.h $(OBJ_DIR)/symmetricBunchPattern.h $(OBJ_DIR)/bunchPatternFromFile.h

$(OBJ_DIR)/bunchTrain.h: bunchTrain.nl
	$(SDDS_REPO)/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables bunchTrain.nl $(OBJ_DIR)/bunchTrain.h

$(OBJ_DIR)/clinchor.h: clinchor.nl
	$(SDDS_REPO)/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables clinchor.nl $(OBJ_DIR)/clinchor.h

$(OBJ_DIR)/generalBunchPattern.h: generalBunchPattern.nl
	$(SDDS_REPO)/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables generalBunchPattern.nl $(OBJ_DIR)/generalBunchPattern.h

$(OBJ_DIR)/symmetricBunchPattern.h: symmetricBunchPattern.nl
	$(SDDS_REPO)/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables symmetricBunchPattern.nl $(OBJ_DIR)/symmetricBunchPattern.h

$(OBJ_DIR)/bunchPatternFromFile.h: bunchPatternFromFile.nl
	$(SDDS_REPO)/bin/$(OS)-$(ARCH)/nlpp -suppressSummaryVariables bunchPatternFromFile.nl $(OBJ_DIR)/bunchPatternFromFile.h
