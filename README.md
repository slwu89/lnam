# lnam

The Linear Network Autocorrelation Model is:

$$ y = W_{1} y + X \beta + \epsilon $$

$$ \epsilon = W_{2} \epsilon + \nu $$

and where the matrices $W_{1}$ and $W_{2}$ include the "AR" and "MA" parameters:

$W_{1} = \sum_{i=1}^{p} \rho_{1,i} W_{1,i}$ and $W_{2} = \sum_{i=1}^{q} \rho_{2,i}W_{2,i}$

  * the code from sna: https://github.com/cran/sna/blob/master/R/models.R#L942
  * description for easy reading: https://rdrr.io/cran/sna/man/lnam.html
  * my comparison with purely spatial regression approaches https://gist.github.com/slwu89/f4ed5a57f5490ec52648bbd687a6a3e7#file-netreg-r

To-do:

  * basic version (direct translation of R)
  * JuMP version for MLE
  * RXInfer https://rxinfer.ml/
  * Turing https://turinglang.org/stable/

Ideas:

  * Static arrays as much as possible. Not sure if we can do static sparse matrices, but that would be a huge speed improvement.