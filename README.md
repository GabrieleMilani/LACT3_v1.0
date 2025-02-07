# Limit Analysis Code Tilting Table Tests ($\rm LACT^3$)

A MATLAB application for fast tilt table test of 2D rigid block structures based on limit analysis.  

## General Info

$\rm LACT^3$ is a MATLAB®-based graphical user interface useful to determine the ultimate load-bearing capacity of masonry walls in-plane loaded. The latter are modelled within a heterogeneous approach where blocks are assumed infinitely resistant and joints are reduced to interfaces exhibiting an associated rigid-plastic behaviour, ruled by a Mohr-Coulomb failure criterion. At incipient collapse, the mechanical problem can be therefore described by means of the two classic limit analysis theorems. $\rm LACT^3$ extracts the geometry of the masonry wall directly from a dxf file, with a precise description of block dimension and shape. The upper bound theorem of limit analysis is used, and the collapse tilting angle is determined along with the corresponding failure mechanism, re-cursively solving a linear programming problem at progressively increased values of rotation of the tilting table. Additionally, the self-dual linear programming approach enables the evaluation of internal actions. The proposed tool is highly user-friendly, requiring only a basic knowledge of CAD software, and is easily manageable as it requires only two mechanical parameters for the joints: (i) friction angle and (ii) cohesion. $\rm LACT^3$ provides an efficient means for the rapid assessment of 2D masonry structures under horizontal loads.

## Setup instructions

To run the software, MATLAB 2019b or later is recommended. Up to now, the only platform supported is Windows 64-bit x86. Download the whole package and unzip it into a folder. Open your MATLAB and run the "LACT3.m" to initial the GUI interface.

## Getting started

After running the "GUI_example.m", an interface consisting of four plot windows and a functional panel will be launched. To start the analysis, press the "Select file" button on the "Model input" tab to indicate the path of a user-defined dxf file properly representing the geometric features of the target 2D masonry wall. Then, press the "Plot model" button aside to import and display the model for double check. The basic information of the model will also be reported in the message panel, such as the number of nodes, blocks, and interfaces. Basically, the analysis can be executed by pressing the "Run" button below since all the material parameters have been set as default. The results will be visualized at the right three windows. The summarizing report will be displayed in the "Output Message" panel, where a spinner is available, allowing users to review all the converged collapse frames recorded during the iteration process.

In the "Property" tab, users can set the mechanical properties of the materials, including friction angle and cohesion, according to the specific case. Program initiation settings can be assigned in the "Iteration" tab, in which the configuration of the initial tilting angle, the step increment, the converging tolerance, and the step control algorithm are accessible. In the folder "sample", we also provided several examples that have passed through our test. Enjoy!

For more details and documentation please also refer to: <>

## Release Date

February 06, 2025

## Developed By

Yiwei Hua, Martina Buzzetti, Natalia Pingaro, Luis C.M. da Silva, Gabriele Milani

## Cite As

Y. Hua, M. Buzzetti, N. Pingaro, L.C.M. da Silva, G. Milani. A computerized tool for the kinematic limit analysis of 2D masonry structures failing on a tilting table. <>
