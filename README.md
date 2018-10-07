Copyright
The original package is available at http://infonet.gist.ac.kr/
COPYRIGHT (c) 2018 Heung-No Lee, and sangjun Park
E-mail: heungno@gist.ac.kr, sjpark1@gist.ac.kr

This package is distributed under the terms of the GNU General Public License 3.0. http://www.gnu.org/copyleft/gpl.html

Permission to use, copy, modify, and distribute this software for any purpose without fee is hereby granted, provided that this entire notice is included in all copies of any software which is or includes a copy or modification of this software and in all copies of the supporting documentation for such software. This software is being provided "as is", without any express or implied warranty. In particular, the authors do not make any representation or warranty of any kind concerning the merchantability of this software or its fitness for any particular purpose.

This pacakge provides a fast mixed integer quadratic programming algorithm for solving sparse signal estimation. 
The following codes can be used to reproduce some of simulation results in our paper titled as "Fast Mixed Integer Quadratic Programming for Sparse Signal Esitmation", IEEE Access.

Example1  : It aims to estimate n-dimensional sparse signal from its corresponding m-dimensional measurement vector.
Figure1_2 : It plots both the averaged mean square error and the averaged support set error of ADM-MIQP, respectively,
            depending on the number of iterations
Figure3   : It plots the phase transitions of the following methods such as ADM-MIQP, MDAL and YALL1.
Figure4_5 : It plots both the averaged mean square error and the averaged support set error of ADM-MIQP, respectively,
            by varying the sparsity level k.
Figure6   : It plots the averaged running time of the following methods such as ADM-MIQP, MDAL and YALL1,
            by varying the problem dimension n.
