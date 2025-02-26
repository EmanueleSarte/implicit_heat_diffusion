Simulating the diffusion of heat in a 2D grid with the Euler implicit method.

(Work In Progess)

## Running the executable

In order to compile the source, you need (tested on Ubuntu 22):
```
sudo apt-get install libblas-dev liblapack-dev  #LAPACK
sudo apt-get install freeglut3-dev # OpenGL/GLUT
sudo apt-get install libglew-dev  # GLEW
```
Then, if everything is in the right place
```
g++ exam_opengl.cpp -o exam_opengl -llapack -lGL -lglut -lGLEW -Ofast
```

To run the executable just
```
./exam_opengl
```

Possible arguments:
TODO

## Preview
The white circle at the start is an infinite source of heat fixed at 300K, the boundary conditions are fixed at 0K.

![heat](./media/preview.gif)

## Commands

**S** - To toggle the start and stop of the animation

**R** - To restart the animation

**P** - To toggle between the pixelated mode

**N** - To step of one iteration (it makes sense if the animation is stopped)

**NUMPAD-N** To toggle the heat contour lines (N between 0 and 9, it makes sense when the pixelated mode is not active)

**1..6** Use the numbers to toggle between different color palettes

Clicking on a point will show the temperature in the terminal.

## Explanation and Algorithm:
Read [explanation.md](./explanation.md)
