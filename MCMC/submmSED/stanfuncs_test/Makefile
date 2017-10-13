

CXX = clang++

PYSTAN = /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/pystan
STAN_MATH = $(PYSTAN)/stan/lib/stan_math

INC = -I/Users/jaffe/home/dev/stan/stan/src/ \
-I${STAN_MATH} -I${STAN_MATH}/lib/boost_1.64.0 \
-I${STAN_MATH}/lib/eigen_3.3.3/ \
-I${STAN_MATH}/lib/cvodes_2.9.0/include

CXXFLAGS = $(INC) -fPIC # -Wall -Wextra -O2
LDFLAGS = -dynamiclib
RM = rm -f

TARGET = ./blackbody.dylib
SOURCE = ./blackbody.cpp

OBJS = $(SOURCE:.cpp=.o)

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

$(OBJS): $(SOURCE)
	$(CXX) $(CXXFLAGS) -c -o $@ $<


.PHONY: clean
clean:
	-${RM} ${OBJS} #$(SOURCE:.cpp=.d)

clean-all:
	-${RM} ${TARGET} ${OBJS}


### see https://gist.github.com/xuhdev/1873316 for a version that makes
###     intermediate .d files to track includes
