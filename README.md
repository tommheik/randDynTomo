# cylSh-inverse-learning-rates

## Inverse learning of randomly sampled sparse dynamic tomography
Numerical methods for estimating the Bregman distance decay rate using cylindrical shearlet regularization for dynamic tomography.

This algorithm is used to produce the numerical results in  
T. A. Bubba, T. Heikkilä, D. Labate, L. Ratti and ?J. P. Rodriguez Ayllon?,  
"Super snappy title",  
_SIAM_ (2023). 

## Credits

Original algorithm by **Tatiana A. Bubba** and **Luca Ratti**, used in  
- T. A. Bubba, M. Burger, T. Helin, and L. Ratti, "Convex regularization in statistical inverse learning problems", _Inverse Problems and Imaging_, 17, pp. 1193–1225 (2023).  
- T. A. Bubba and L. Ratti, "Shearlet-based regularization in statistical inverse learning with an application to x-ray tomography", 
_Inverse Problems_, p. 054001, (2022).

The minimization problem is solved with the Variable Metric Inexact Line-search Algorithm (**VMILA**) introduced in  
- S. Bonettini, I. Loris, F. Porta and M. Prato, "Variable metric inexact line-search based methods for nonsmooth optimization", _SIAM J. Optim._, 26 891–921, (2016).

The proximal operators require solving additional minimization problems. These are computed using the Scaled Gradient Projection (**SGP**) algorithm introduced in  
- S. Bonettini, R. Zanella and L. Zanni, "A scaled gradient projection method for constrained image deblurring", 
_Inverse Problems_ 25 015002, (2006).

The 3D Discrete Wavelet Transform is part of the **Wavelet Toolbox** and the **3D cylindrical shearlets** are based on
- G. R. Easley, K. Guo, D. Labate and B. R. Pahari, "Optimally sparse representations of cartoon-like cylindrical data", 
__The Journal of Geometric Analysis__ 31: 8926-8946, (2021).  
https://github.com/tommheik/3d_cylind_shear

The tomography forward operator is computed using the **ASTRA Toolbox** 
- W. Van Aarle, W. J. Palenstijn, J. Cant, E. Janssens, F. Bleichrodt, A. Dabravolski, J. De Beenhouwer,
K. J. Batenburg and J. Sijbers, "Fast and flexible x-ray tomography using the ASTRA toolbox", 
_Opt. Express_ 24 25129–47 (2016).
- W. Van Aarle, W. J. Palenstijn, J. De Beenhouwer, T. Altantzis, S. Bals, K. J. Batenburg and J. Sijbers, "The ASTRA toolbox: a platform for advanced algorithm development in electron tomography", _Ultramicroscopy_ 157 35–47, (2015).

The algorithm also uses **Spot- A Linear-Operator Toolbox** and **HelTomo**
- E. Van den Berg and M. P. Friedlander, "Spot-a linear-operator toolbox", v1.2, (2013) http://cs.ubc.ca/labs/scl/spot
- A. Meaney, "HelTomo - Helsinki Tomography Toolbox", v2.0.0, (2022) https://github.com/Diagonalizable/HelTomo