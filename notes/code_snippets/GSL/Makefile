CC     = gcc
CFLAGS = -Wall -g
LIBS   = -lz -lm -lrt -ldl
SRCS   = gamma.c gamma_inc.c erfc.c expint.c pow_int.c beta_inc.c beta.c \
	error.c exp.c  log.c psi.c trig.c stream.c math.c zeta.c infnan.c elementary.c \
	fdiv.c main.c
OBJS   = $(SRCS:.c=.o)

EXECUTABLE := gamma

# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'

.PHONY: clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXECUTABLE)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	-rm -f $(EXECUTABLE) *.o *~ core

