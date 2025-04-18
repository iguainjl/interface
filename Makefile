CXX = nvcc

TAU?=1.0

INCLUDES = -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.7/math_libs/12.2/include 
FLAGS = --expt-extended-lambda -lcufft -std=c++14 -arch=sm_75 
PARAMS = -DC2=1.0 -DC4=0.0 -DTAU=$(TAU) #-DDOUBLE  

LDFLAGS = -L/opt/nvidia/hpc_sdk/Linux_x86_64/23.7/math_libs/12.2/lib64 


activeinterface: ew.cu
	$(CXX) $(FLAGS) $(PARAMS) ew.cu -o activeinterface $(LDFLAGS) $(INCLUDES) 


update_git:
	git add *.cu Makefile *.h *.sh README.md ; git commit -m "program update"; git push

clean:
	rm activeinterface
