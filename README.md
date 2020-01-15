# Sheaf stability

This program computes whether a given sheaf is semistable following the algorithm laid out in the paper (Computing stability of sheafs - Holger Brenner and Jonathan Steinbuch (TODO link).

# Usage
To compile the program you need to have the GNU Multiprecision Library (https://gmplib.org/) and the Boost Program Options (https://www.boost.org/) installed.

To compile just run the command 
```
make
```
in the main directory.

The basic usage is as follows. In the directory with the executable run:
```
./stability --input-file="input.txt" --output-file="output.txt" 
```
There are several command line options. In most cases going with the default should be fine.

There are the following accessible routines:
## Semistability
This is the main routine. To use it provide an input file of the following form:
```
semistability
characteristic: 0
variables: "x", "y", "z"
relations: z^5+y^5+x^5
matrix: {{x^4+y^2z^2,y^4,z^4,x^7}}
```
In this case we compute whether the sheaf given by Syz(x^4+y^2z^2,y^4,z^4,x^7) over \Q[x,y,z]/z^5+y^5+x^5 is semistable.

The output file will contain the following:
```
1
{ExteriorPower => 1, SymmetricPower => 4, Twist => 25}
```
The 1 in the first line means that the sheaf is not semistable. The next line describes where a destabilizing subsheaf has been found.

## Powers
With this we can directly get Sym^q(Ext^s(F)) of a kernel sheaf F. 

To use it provide an input file of the following form:
```
powers
characteristic: 3
variables: "x", "y", "z"
relations: z^5+x^5+y^5
matrix: {{x^4+y^2z^2,y^4,z^4,x^7}}
symmetricpower: 2
exteriorpower: 1
```
In this case we compute Sym^2(Ext^1(Syz(x^4+y^2z^2,y^4,z^4,x^7))) over \F_3[x,y,z]/z^5+y^5+x^5. The output will be 
```
    |         8       8       8          11   8  8          11   8          11        14 
----| --------- ------- ------- ----------- --- -- ----------- --- ----------- --------- 
   4| 2z2y2+2x4      y4      z4 2z5x2+2y5x2   0  0           0   0           0         0 
   4|         0 z2y2+x4       0           0 2y4 z4 2z5x2+2y5x2   0           0         0 
   4|         0       0 z2y2+x4           0   0 y4           0 2z4 2z5x2+2y5x2         0 
   7|         0       0       0     z2y2+x4   0  0          y4   0          z4 z5x2+y5x2 
```
The first line and row describe the negatives of the twists making up the sheafs in the kernel sequence. 
Consequently the symmetric power sheaf is the kernel of the map O(-8)^6+O(-11)^3+O(-14)^1-->O(-4)^3+O(-7)^1 given by the output matrix.

## Kernel
With this we can directly get global sections of Sym^q(Ext^s(F))(t) of a kernel sheaf F. 

To use it provide an input file of the following form:
```
kernel
characteristic: 0
variables: "x", "y", "z"
relations: z^5+x^5+y^5
matrix: {{x^4+y^2z^2,y^4,z^4,x^7}}
symmetricpower: 2
exteriorpower: 1
twist: 50
```
In this case we compute global sections of Sym^2(Ext^1(Syz(x^4+y^2z^2,y^4,z^4,x^7)))(50) over \F_3[x,y,z]/z^5+y^5+x^5. The output will be 
```
1
    |                                                                                                   50 
----| ---------------------------------------------------------------------------------------------------- 
   8|                                     -5z32y9x-6z31y8x3-4z30y12-2z29y11x2+3z28y10x4+2z26y13x3+z24y16x2 
   8| -10z37y5-12z36y4x2+6z34y7x+18z33y6x3-2z32y10-4z31y9x2+2z30y8x4-2z29y12x+2z28y11x3+2z26y14x2+2z24y17x 
   8|                                                                                                    0 
  11|                                                                                                    0 
   8|                            -6z41x+8z39y3+15z38y2x2+5z37yx4-11z36y5x-8z35y4x3+3z34y8+12z33y7x2+z24y18 
   8|                                                                                                    0 
  11|                                                                                                    0 
   8|                                                                                                    0 
  11|                                                                                                    0 
  14|                                                                                                    0 
```
The first line of the output file is 0 if there are no global sections. If you specifiy the program option to compute the full kernel (with the option -f) it will have the dimension of the space of global sections. By default it's just one if there is a global section. 

The matrix after that is in the same format as for the powers. The rows are linearly independent global sections.

# Notes
The polynomials are output in a short form where an integer directly after a variable means that that variable is taken to the power of the integer. So for example -5z32y9x stands for -5*z^32*y^9*x^1. If you prefer a longer form there is a command line option (-p) for that.
You can also enter polynomials in short form in the input file if you want. The parser is relatively robust in that regard.

You can input variables in many forms. The only restrictions are that the first character has to be a letter from the alphabet and that a variable descriptor can not be contained in full in another. So "z_{0,1}" is a totally fine variable descriptor.

