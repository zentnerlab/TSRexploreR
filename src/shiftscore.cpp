#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double ShiftScoreFast(vec x, vec y, int k, int xn, int yn, uvec w){
  // Inputs are unnormalized histograms, evenly spaced with 0 as necessary.
  // This is a utility function to avoid repetitive calculations
  vec px = x / xn;
  vec py = y / yn;
  
  vec weighted = (px - py) % w;
  double ans = sum(weighted) / k;
  return(ans);
}


vec ShiftScore(sp_mat x, sp_mat y, int calcP, int nresamp){
  // Inputs are sparse vectors of the same length with counts, unnormalized

  int xn = accu(x);
  int yn = accu(y);
  
  x /= xn;
  y /= yn;
  
  sp_mat w = (x - y);
  vec dense_weights = nonzeros(w);
  uvec pos = find(w);
  pos += 1;
  
  vec out = dense_weights % pos;
  vec ans(2);
  ans.zeros();
  ans(0) += accu(out) / y.n_elem;
  
  if(calcP){
    double pp = 0.0;
    sp_mat pxy = (xn*x + yn*y) / (xn + yn);
    vec dense_probs = nonzeros(pxy);
    uvec pos_pxy = find(pxy);

    int k = dense_probs.n_elem;
    ivec simx(k);
    ivec simy(k);
    double sim = 0.0;
  
    RNGScope scope;
    for(int i = 0; i < nresamp; i++){
      rmultinom(xn, dense_probs.begin(), k, simx.begin());
      rmultinom(yn, dense_probs.begin(), k, simy.begin());
      vec sx = conv_to<vec>::from(simx);
      vec sy = conv_to<vec>::from(simy);
      // Rcout << "sx " << sx << std::endl;
      sim = ShiftScoreFast(sx, sy, k, xn, yn, pos_pxy);
      pp += (abs(sim) > abs(ans(0)));
    }
    ans(1) += pp / nresamp;
  }
  
  
  return(ans);
}


// [[Rcpp::export]]
mat allTheShiftScores(CharacterVector fhash, uvec dists, vec scores, vec sample, 
                      int calcP, int nresamp, int ntests){
  mat ans(2,ntests);
  
  // loop to find the sequence starts
  int startix = 0;
  int endix = 0;
  int n = fhash.length();
  int ntest = 0;
  for(int it = 1; it < n; it++){
    // Rcout << "it = " << it << std::endl;
    if(fhash(it-1) != fhash(it) || it == n-1){
      // Rcout << "startix = " << startix << std::endl;
      endix = it - 1;
      // Rcout << "endix = " << endix << std::endl;
      vec samps = sample.subvec(startix,endix);
      vec vals = scores.subvec(startix,endix);
      uvec samp1 = find(samps);
      uvec samp0 = find(samps<0.5);
      uvec dists_sub = dists.subvec(startix,endix);
      int nsp = dists_sub.max() + 1;
      int n1 = samp1.size();
      int n0 = samp0.size();
      // Rcout << "n1 = " << n1 << std::endl;
      // Rcout << "n0 = " << n0 << std::endl;
      umat locx(2,n1,fill::zeros);
      umat locy(2,n0,fill::zeros);
      locx.row(0) = dists_sub.elem(samp1).t();
      locy.row(0) = dists_sub.elem(samp0).t();
      // Rcout << locy << std::endl;
      // Rcout << vals.elem(samp0) << std::endl;
      // Rcout << nsp << std::endl;
      sp_mat x(locx, vals.elem(samp1), nsp, 1);
      sp_mat y(locy, vals.elem(samp0), nsp, 1);
      // Rcout << y << std::endl;
      
      ans.col(ntest) = ShiftScore(x, y, calcP, nresamp);
      ntest++;
      startix = it;
    }
  }
  return(ans);
}
