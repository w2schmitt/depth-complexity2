OS=$(shell uname -s)

GLFW_CFLAGS = $(shell pkg-config --cflags libglfw)
GLFW_LDFLAGS = $(shell pkg-config --libs libglfw)

ATB_CFLAGS = 
ATB_LDFLAGS = -lAntTweakBar

ifeq "$(OS)" "Darwin"
	GLUT_LDFLAGS += -framework GLUT -framework OpenGL
	GLEW_CFLAGS = -I/opt/local/include
	GLEW_LDFLAGS = -L/opt/local/lib -lGLEW
else
    ifeq "$(OS)" "Linux"
	    GLUT_LDFLAGS = -lglut -lGLU -lGL
    else
	    GLUT_LDFLAGS = -lglut
    endif
	GLEW_LDFLAGS = -lGLEW
endif

CPPFLAGS = $(GLFW_CFLAGS) $(ATB_CFLAGS) $(GLEW_CFLAGS) -Wall -Werror -O3
#CPPFLAGS = $(GLFW_CFLAGS) $(ATB_CFLAGS) $(GLEW_CFLAGS) -Wall -Werror -g
LDFLAGS = $(ATB_LDFLAGS) $(GLEW_LDFLAGS) -Wall -Werror
#LDFLAGS = $(ATB_LDFLAGS) $(GLEW_LDFLAGS) -Wall -Werror -g


.PHONY: all clean

#all : sphereray main_joao depthcomplexity2d depthcomplexity3d_offline depthcomplexity3d \
	  #randomdepthcomplexity3d randomdepthcomplexity3d_offline newrandomdepthcomplexity3d mergeOFF

#all : depthcomplexity2d depthcomplexity3d_offline depthcomplexity3d \
      newrandomdepthcomplexity3d mergeOFF

all : newrandomdepthcomplexity3d depthcomplexity3d_offline thickness_2d mergeOFF
	
	
clean:
	rm -f *.o camera/*.o depthcomplexity3d thickness_2d depthcomplexity2d main_joao \
	test_intersection test depthcomplexity3d_offline randomdepthcomplexity3d \
	randomdepthcomplexity3d_offline newrandomdepthcomplexity3d mergeOFF
	
hard_clean: clean
	rm -f *.off hist*.txt *.svg toPlot.txt
	rm -rf Tests
	
#depthcomplexity3d: depthcomplexity3d.o camera/Camera.o camera/float3.o dc_2d.o dc_3d.o util.o ShaderMgr.o timer.o
#	g++ -g -o $@ $^ $(LDFLAGS) $(GLFW_LDFLAGS) $(GLUT_LDFLAGS)

#depthcomplexity2d: depthcomplexity2d.o dc_2d.o util.o ShaderMgr.o flags.o timer.o 
#	g++ -g -o $@ $^ $(LDFLAGS) $(GLUT_LDFLAGS)

depthcomplexity3d_offline: \
	RFDepthComplexity3D.o camera/Camera.o camera/float3.o util.o ShaderMgr.o Texture3D.o timer.o flags.o
	g++ -g -o $@ depthcomplexity3d.cpp -D "COMPUTE_OFFLINE" $^ $(LDFLAGS) \
	$(GLFW_LDFLAGS) $(GLUT_LDFLAGS) $(GLEW_CFLAGS) -lpthread

thickness_2d: thickness_2d.o
	g++ -g -o $@ $^ $(LDFLAGS) $(GLUT_LDFLAGS)
#main_joao: main_joao.o
#	g++ -g -o $@ $^ $(LDFLAGS) $(GLUT_LDFLAGS)

#sphereray: sphererays.o
#	g++ -g -o $@ $^ $(LDFLAGS) $(GLUT_LDFLAGS)
  
#randomdepthcomplexity3d: dc_3d_random.o camera/Camera.o camera/float3.o util.o timer.o flags.o
#	g++ -g -o $@ depthcomplexity3d.cpp -D "USE_RANDOM_DC3D" $^ $(LDFLAGS) \
	$(GLFW_LDFLAGS) $(GLUT_LDFLAGS) $(GLEW_CFLAGS)
	
newrandomdepthcomplexity3d: \
	RFDepthComplexity3D.o camera/Camera.o camera/float3.o util.o ShaderMgr.o Texture3D.o timer.o flags.o
	g++ -g -o $@ depthcomplexity3d.cpp $^ $(LDFLAGS) \
	$(GLFW_LDFLAGS) $(GLUT_LDFLAGS) $(GLEW_CFLAGS) -lpthread
	
#randomdepthcomplexity3d_offline: dc_3d_random.o flags.o util.o timer.o flags.o
#	g++ -g -o $@ depthcomplexity3d_offline.cpp -D "USE_RANDOM_DC3D" $^ \
	$(LDFLAGS) $(GLFW_LDFLAGS) $(GLUT_LDFLAGS) $(GLEW_CFLAGS)

mergeOFF: mergeOFF.o
	g++ -g -o $@ $^
	
