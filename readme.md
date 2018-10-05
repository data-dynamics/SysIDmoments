This code is the implementation of experiment in paper 
Ozay, Necmiye, Constantino Lagoa, and Mario Sznaier. "Set membership identification of 
switched linear systems with known number of subsystems." Automatica 51 (2015): 180-191.

Parameter Explanation:
* Line 44---N         Number of time-steps
* Line 47---sw_seq    Label of test case
* Line 48---num_sys   Number of subsystems
* Line 53---p1,p2,p3  Coefficients of subsystem through the time horizon
* Line 102--w         Noise
* Line 106--u         Input
* Line 108--m         Order of system
* Line 112--y         Output of system

The default hybrid system is given by
* I.   y_t = 0.2y_(t-1)+0.24y_(t-2)+2u_(t-1)+ w_t, poles at -0.4, 0.6
* II.  y_t = -1.4y_(t-1)-0.53y_(t-2)+u_(t-1)+w_t, poles at -0.7+0.2i, -0.7-0.2i
* III. y_t = 1.7y_(t-1)-0.72y_(t-2)+0.5u_(t-1)+w_t, poles at 0.9, 0.8

7 test cases with different switching sequences are given in Line 52-59. And they can be
selected by setting the value of sw_seq from 1 to 7. The default sequence (used in paper)
corresponds to sw_seq=1.


Authors: Necmiye Ozay & Zhe Du, University of Michigan

Research partly supported by DARPA grant N66001-14-1-4045.