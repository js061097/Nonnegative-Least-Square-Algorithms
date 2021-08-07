# Nonnegative-Least-Square-Algorithms (NNLS) 
  Contains the Lawson &amp; Hanson active set, fnnls by Bro &amp; Jong, Projected Quasi-Newton and Random projections NNLS algorithms, all implemented in MATLAB.
  Thesis: url (https://drive.google.com/file/d/1VfRrUy50YOgNNjRS9-6904cNzgLc7aec/view)

## Applications of NNLS
 &nbsp;Can be used to solve any overdetermined system of equations where the solution vector is constrained to take non-negative values. &nbsp;For example, in image processing we require the pixel values of image stored as a vector to be positive.<br/>
  Overdetermined system: Ax=b, where A = mxn real valued system  with m>n.
  
## About the algorithms (References) 
  &nbsp;The theory behind the algorithms can be found here:
  1. <b>Lawson & Hanson active set algorithm:</b> [C. L. Lawson and R. J. Hanson, Solving least squares Problems, Prentice Hall, 1987.](https://en.wikipedia.org/wiki/Non-negative_least_squares) <br/>
  2. <b>fnnls by Bro & Jong:</b> [R. Bro, S. D. Jong, A fast non-negativity-constrained least squares algorithm, Journal of Chemometrics, Vol. 11, No. 5, (1997), pp. 393–401.](http://xrm.phys.northwestern.edu/research/pdf_papers/1997/bro_chemometrics_1997.pdf) <br/>
  3. <b>Projected Quasi-Newton:</b> [D. Kim, S. Sra, and I. S. Dhillon, A new projected quasi-Newton Approach for the non-negative least squares problem, Technical Report TR-06-54, Computer Sciences, The Universityof Texas at Austin, (2006).](https://www.cs.utexas.edu/ftp/techreports/tr06-54.pdf) <br/>
  4. <b>Random Projections:</b> [C. Boutsidis and P. Drineas, Random projections for nonnegative least squares, Linear Algebra Appl., 431 (2009), pp. 760–771.](https://www.sciencedirect.com/science/article/pii/S0024379509001633)

## Tests and Results
 &nbsp;The algorithms were tested for randomly generated matrix problems with sizes varying from 2800x2000 to 10400x6800. The computational time taken to solve the problem and the residual error of the results were compared. The computational time measured using cputime and the residual error ||Ax-b|| were measured in MATLAB were obtained as follows:
 
 ![alt text](https://github.com/js061097/Nonnegative-Least-Square-Algorithms/blob/main/Results/cputime.jpg)
 ![alt text](https://github.com/js061097/Nonnegative-Least-Square-Algorithms/blob/main/Results/RelativeNorms.jpg)
 
 ## Conclusion
  - Generally, the computational time taken to solve a large system would be in the order:</br> 
  &nbsp;&nbsp;&nbsp;Lawson & Hanson < fnnls < Quasi-Newton < Randomized.
  - Randomized algorithm using Quasi-Newton method will produce results that come with some relative errors in comparision to the others but they solve the problems 10 times faster than the Quasi-Newton method for matrix sizes upto 10400x6800.
