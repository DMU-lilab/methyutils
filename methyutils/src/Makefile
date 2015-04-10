CC = gcc

CFLAGES = -Wall -fPIC
LFLAGES = -lm

OBJS = chisq.o statistics.o
DEPS = %.h
TARGET = statistics.so

all : $(OBJS)
	R CMD SHLIB -o $(TARGET) $(OBJS) $(LFLAGS)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $<

clean:
	-rm *.o *.so
