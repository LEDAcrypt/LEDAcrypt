ifndef SL
$(warning -------------------------------------------------------------------- )
$(warning Missing definition for the security level, defaulting to SL=1. Define)
$(warning it as an extra make parameter SL=<desired level> where the desired   )
$(warning level is an integer between 1 and 5)
$(warning -------------------------------------------------------------------- )
SL = 1
endif


ifneq ($(SL),1)
ifneq ($(SL),2)
ifneq ($(SL),3)
ifneq ($(SL),4)
ifneq ($(SL),5)
$(warning -------------------------------------------------------------------- )
$(warning Invalid definition for the security level,defaulting to SL=1. )
$(warning Define it as an extra make parameter as ) 
$(warning SL=<desired level> N0=<desired number of blocks> make )
$(warning where the desired security level --category-- is an integer between )
$(warning 1 and 5)
$(warning -------------------------------------------------------------------- )
SL = 1
endif
endif
endif
endif
endif

ifndef DFR_SL_LEVEL
$(warning -------------------------------------------------------------------- )
$(warning Missing definition for the desired DFR , defaulting )
$(warning to 2^-64. Define it as an extra make parameter as ) 
$(warning SL=<desired level> DFR_SL_LEVEL=<0 or 1> make )
$(warning where 0 indicates 2^-64 DFR and 1 2^-SL DFR)
$(warning -------------------------------------------------------------------- )
DFR_SL_LEVEL = 0
endif

ifneq ($(DFR_SL_LEVEL),0)
ifneq ($(DFR_SL_LEVEL),1)
$(warning -------------------------------------------------------------------- )
$(warning Invalid definition for the desired DFR , defaulting )
$(warning to 2^-64. Define it as an extra make parameter as ) 
$(warning SL=<desired level> DFR_SL_LEVEL=<0 or 1> make )
$(warning where 0 indicates 2^-64 DFR and 1 2^-SL DFR)
$(warning -------------------------------------------------------------------- )
DFR_SL_LEVEL = 0
endif
endif

N0 = 2
  CC = gcc 

  CFLAGS = -DCATEGORY=$(SL) -DN0=$(N0) -DDFR_SL_LEVEL=$(DFR_SL_LEVEL) -DCPU_WORD_BITS=64 \
           -std=c11 -Wall -pedantic -Wmaybe-uninitialized -Wuninitialized \
	   -march=native -O3 -g3
  LDFLAGS = -lm -lkeccak
  INCLUDES = -I./include
  SRCDIR = library
  OBJDIR = bin

# Gathers the names of all C files 
    CSRC = $(wildcard $(SRCDIR)/*.c) 
# Produces in $(COBJS) the names of .o object files for all C files
   COBJS = $(CSRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

$(OBJDIR)/%.o: $(SRCDIR)/%.c 
	$(CC) -c $(CFLAGS) $(INCLUDES)  $< -o $@ 

.PHONY : all clean

all: $(COBJS) 
	ar rcs $(OBJDIR)/libLEDAkem_sl$(SL)_N0$(N0).a $(COBJS) 
clean:
	$(RM) $(OBJDIR)/*.o  $(OBJDIR)/libLEDAkem_sl$(SL)_N0$(N0).a
