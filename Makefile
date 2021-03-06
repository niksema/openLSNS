# Location of the CUDA Toolkit
CUDA_PATH     ?= /usr/local/cuda-7.0

# architecture & os
HOST_ARCH     := $(shell uname -m)
TARGET_ARCH   ?= $(HOST_ARCH)
TARGET_SIZE   := 
HOST_OS       := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
TARGET_OS     ?= $(HOST_OS)
include Makefile.host

# host compiler & flags
HOST_COMPILER ?= g++
NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)
NVCCFLAGS     := -m${TARGET_SIZE} --use_fast_math
CCFLAGS       := 
LDFLAGS       :=
GENCODE_FLAGS := -gencode arch=compute_20,code=compute_20
BUILD_TYPE    :=
include Makefile.flags

ALL_CCFLAGS   :=
ALL_CCFLAGS   += $(NVCCFLAGS)
ALL_CCFLAGS   += $(EXTRA_NVCCFLAGS)
ALL_CCFLAGS   += $(addprefix -Xcompiler ,$(CCFLAGS))
ALL_CCFLAGS   += $(addprefix -Xcompiler ,$(EXTRA_CCFLAGS))

ALL_LDFLAGS   :=
ALL_LDFLAGS   += $(ALL_CCFLAGS)
ALL_LDFLAGS   += $(addprefix -Xlinker ,$(LDFLAGS))
ALL_LDFLAGS   += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))

# Common includes and paths for CUDA
BUILD_ENABLED := 1
OUTDIR        := bin
COREDIR       := core
CUDADIR       := cuda
LIBRARIES     :=
INCLUDES      := -I./src/include -I./src/include/cuda -I./src/core -I./src/code/cuda
DEPS          := info.h kernel.h config.h cuda_timer.h lsnsmath.h gates.h gateproc.h
_OBJ          = info.o kernel.o gates.o
OBJ_          = $(patsubst %,$(OUTDIR)/%,$(_OBJ))
include Makefile.cuda
ifeq ($(BUILD_ENABLED),0)
	EXEC ?= @echo "[@]"
endif

# Target rules
all: build

build: lsns

check.deps:
$(OUTDIR)/%.o: %.cpp $(DEPS)
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<
$(OUTDIR)/%.o: $(COREDIR)/%.cpp $(DEPS)
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<
$(OUTDIR)/%.o:%.cu $(DEPS) 
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<
$(OUTDIR)/%.o: $(CUDADIR)/%.cpp $(DEPS)
	$(EXEC) $(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

lsns: $(OBJ_)
	$(EXEC) $(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -o $@ $+ $(LIBRARIES)

run: build
	$(EXEC) ./lsns

clean:
	rm -f $(OBJ_)

clobber: clean
