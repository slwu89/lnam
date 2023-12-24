# lnam

The Linear Network Autocorrelation Model is:

$$ y = W_{1} y + X \beta + \epsilon $$

$$ \epsilon = W_{2} \epsilon + \nu $$

and where the matrices $W_{1}$ and $W_{2}$ include the "AR" and "MA" parameters:

$W_{1} = \[ \sum_{i}^{p} \rho_{1,i} W_{1,i} \]$ and $W2 = sum( rho2_i W2_i, i=1..q)$

  * the code from sna: https://github.com/cran/sna/blob/master/R/models.R#L942
  * description for easy reading: https://rdrr.io/cran/sna/man/lnam.html
