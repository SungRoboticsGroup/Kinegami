# Kinegami
Given a Denavit-Hartenberg (D-H) specification of a kinematic chain robot, the program generates a crease pattern that folds into a kinematically equivalent robot with compliant joints. The program takes in the D-H specification and assigns the corresponding joint to a location that is sufficiently far from other joints while keeping the same kinematic properties. Then an origami link inspired by the Dubin's path method is created to connect every two conseccutive joints. 

To construct the D-H specification, follow the variable definitaion ("Link length (a)", "Link twist (alpha)", "Joint offset (d)", and "Joint angle (theta)") and numbering system of this kinematic chain mechanism schematic drawing (for details, please read "W. Khalil and E. Dombre, *Modeling, identification and control of robots*") to form a table.

<!-- ![DH](https://user-images.githubusercontent.com/50150425/161108095-26ed20c8-596d-4ba4-a642-271e6f2d4c32.png) -->
<!-- ![DHlight](https://user-images.githubusercontent.com/50150425/161122307-ad8ce29e-18ea-4b91-883e-ca186c5232fc.png) -->
<img src="https://user-images.githubusercontent.com/50150425/161122307-ad8ce29e-18ea-4b91-883e-ca186c5232fc.png" alt="DHlight" width="600"/>

Our Kinegami algorithm recruits a catalogue of parameterized modules. The folded state of the origami module, its spacial operator representation, and its crease pattern is shown here: (a) the origami prism tube, (b) the twist fitting, (c) the elbow fitting, (d) the prismatic joint, (e) the revolute joint, (f) partial close-up of an elbow fitting and (g) partial close-up of a revolute joint.

![OrigamiModuleNew](https://user-images.githubusercontent.com/50150425/161108362-0ee75174-3fda-47f5-8d4d-f7433aadc7c7.png)

In addition, it provides the additional contribution of automatically choosing the relevant modules (with design parameters that control their compliance) and composing them into a non-self-intersecting single sheet pattern, thus reducing the design problem into simply one of abstract specification.
The resulting pipeline does not require any additional human input beyond the D-H specification, though its algorithmic steps are sufficiently transparent as to facilitate the integration of designers' alternative modules or more suitably optimized compositions when desired.

## User Guide
Run scripts `Kinegami_******.m` for existing examples and change parameters if desired. 
To create your kinematic chain robot, fill out the `Kinegami_Template.m` file in this order:
1. Design the regular polygon shape (the number of sides "nside" and circumradius "r") as the tubular origami base
2. Determine the number of joints for your robot, fill in the number of joints plus one (the fingertip) in variable "n".
3. Fill in the DH parameters specifications for variable "D" following the joint design.
4. Specify joints information, all variables that contain "??", including TYPE, maximum joint range, initial joint configurations, etc.
5. In addition, modify the user specifications for crease pattern generation including DXF printing and segmentations.

For more comprehensive understanding, reference supporting functions. The algorithm requires "fSolve" from the MatLab Optimization Toolbox. Please install the add-on Optimization Toolbox.

## Laser cutting instructions
The crease pattern generated by the algorithm can be output as a `.dxf` file, its unit is meters.
The blue lines indicates the mountain folds, the red indicates the valley folds, and the black lines indicates the boarder edges.
We then cut our specimen with the "PLS 4.75" laser cutting machine.
For the blue and red line, we perforated the paper by setting the machine parameters to be (power 10, speed 20, ppi 25).
For the black line, we cut the paper by setting the machine parameters to be (power 10, speed 20, ppi 750).
The paper we use is the 8 mil thick Durilla synthetics paper with polyester finish (CTI Paper, USA).

<!-- ## Updates:
7/5/2021:
Edited `JointAssignment.m` to include correct value of rs for Prismatic Joints.
Edited `Kinegami.m` to support plotting for Proximal and Distal Frames. Added new function frameplot.m for frame plotting. Changed manner in which figures are closed in papercut files. -->

## Acknowledgement
This work was supported in part by the Army Research Office under the SLICE Multidisciplinary University Research Initiatives Program, award under Grant \#W911NF1810327 and in part by the National Science Foundation under Grant \#1845339.

## License
This code is released using the Penn Software License. Please refer to `LICENSE.txt` for details.
