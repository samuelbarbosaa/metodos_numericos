# grid
tauchen_grid = function(N, m, sigma, rho) {
  thetaN = m * sigma / sqrt(1-rho^2)  
  seq(-thetaN, thetaN, length.out = N)
}

# matriz de transicao
tauchen_P = function(grid, rho, sigma, mu) {
  N = length(grid)
  delta = (max(grid) - min(grid)) / (N-1)
  PT = matrix(NA_real_, N, N)
  for(i in 1:N) {
    for(j in 1:N) {
      if(j == 1) {
        PT[i,j] = pnorm( (grid[1] - (1-rho)*mu - rho*grid[i] + delta/2) / sigma )
      } else if(j==N) {
        PT[i,j] = 1 - pnorm( (grid[N] - (1-rho)*mu - rho*grid[i] - delta/2) / sigma )
      } else {
        PT[i,j] = 
          pnorm((grid[j] + delta/2 - (1-rho)*mu - rho*grid[i]) / sigma) -
          pnorm((grid[j] - delta/2 - (1-rho)*mu - rho*grid[i]) / sigma)
      }
    }
  }
  return(PT)
}