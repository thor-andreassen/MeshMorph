Finite Element Morphing, Overclosure Removal, and Slicing Code

This code contains a package of functions to morph template meshes to other geometries using a GRNN based algorithm herein.

The code is thoroughly described in the full paper, available at: www.doi.org/10.1016/j.medengphy.2024.104136

The full paper also includes pseudocode, additional examples, test cases, and results demonstrated in the full manuscript.

This code was developed with another project to make freely available geometries of the lower limbs of the Visibile Human Female and Visible Human Male. Whereby, this code was used in the final stages to remove all overclosures between adjacent geometries. This project and all geometries are available at the following: https://digitalcommons.du.edu/visiblehuman/

The code has been implemented in MATLAB utilizing the Machine Learning Toolbox and is available free and open-source, but we ask that you cite the following two works:

T. E. Andreassen, D. R. Hume, L. D. Hamilton, K. E. Walker, S. E. Higinbotham, and K. B. Shelburne, “Three Dimensional Lower Extremity Musculoskeletal Geometry of the Visible Human Female and Male”, Scientific Data, vol. 10, no. 1, p. 34, Jan. 2023, doi: 10.1038/s41597-022-01905-2.

T. E. Andreassen, D. R. Hume, L. D. Hamilton, S. E. Higinbotham, and K. B. Shelburne, "Automated 2D and 3D Finite Element Overclosure Adjustment and Mesh Morphing Using Generalized Regression Neural Networks", Medical Engineering & Physics, Mar. 2024, doi: 10.1016/j.medengphy.2024.104136.

The main GRNN algorithm is contained in the function GRNNMorph.m, but only contains the GRNN steps of the proposed algorithm and not the mesh reduction steps as outline in the paper. The example file called "Example_2D_to_2D_Mesh_Morph.m" contains an example of the morphing algorithm being used on an example dataset. A user can simply use
the instructions provided to create a new pre-formatted folder structure that can be easily run by the example function and morph the given user's geometries instead of the example case.

Adding changes to the code is encouraged and can be added to the repository by contacting the author. The author will check new or revised content for accuracy and completeness and add to the repository.

Future/ongoing work aims to recreate the code using code that does not need the Machine Learning Toolbox, as well as implementing the code into a Python Toolbox for additional use to a larger audience.
