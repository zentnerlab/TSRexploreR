#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double ShiftScoreFast(arma::vec x, arma::vec y, int k, int xn, int yn, arma::uvec w){
  // Inputs are unnormalized histograms, evenly spaced with 0 as necessary.
  // This is a utility function to avoid repetitive calculations
  arma::vec px = x / xn;
  arma::vec py = y / yn;
  
  arma::vec weighted = (px - py) % w;
  double ans = arma::sum(weighted) / k;
  return(ans);
}


arma::vec ShiftScore(arma::sp_mat x, arma::sp_mat y, int calcP, int nresamp){
  // Inputs are sparse vectors of the same length with counts, unnormalized

  int xn = accu(x);
  int yn = accu(y);
  
  x /= xn;
  y /= yn;
  
  arma::sp_mat w = (x - y);
  arma::vec dense_weights = arma::nonzeros(w);
  arma::uvec pos = arma::find(w);
  pos += 1;
  
  arma::vec out = dense_weights % pos;
  arma::vec ans(2);
  ans.zeros();
  ans(0) += accu(out) / y.n_elem;
  
  if(calcP){
    double pp = 0.0;
    arma::sp_mat pxy = (xn*x + yn*y) / (xn + yn);
    arma::vec dense_probs = arma::nonzeros(pxy);
    arma::uvec pos_pxy = arma::find(pxy);

    int k = dense_probs.n_elem;
    arma::ivec simx(k);
    arma::ivec simy(k);
    double sim = 0.0;
  
    RNGScope scope;
    for(int i = 0; i < nresamp; i++){
      rmultinom(xn, dense_probs.begin(), k, simx.begin());
      rmultinom(yn, dense_probs.begin(), k, simy.begin());
      arma::vec sx = arma::conv_to<arma::vec>::from(simx);
      arma::vec sy = arma::conv_to<arma::vec>::from(simy);
      // Rcout << "sx " << sx << std::endl;
      sim = ShiftScoreFast(sx, sy, k, xn, yn, pos_pxy);
      pp += (abs(sim) > abs(ans(0)));
    }
    ans(1) += pp / nresamp;
  }
  
  
  return(ans);
}


// [[Rcpp::export]]
arma::mat allTheShiftScores(CharacterVector fhash, arma::uvec dists, arma::vec scores, 
                            arma::vec sample, int calcP, int nresamp, int ntests){
  arma::mat ans(2,ntests);
  
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
      arma::vec samps = sample.subvec(startix,endix);
      arma::vec vals = scores.subvec(startix,endix);
      arma::uvec samp1 = arma::find(samps);
      arma::uvec samp0 = arma::find(samps<0.5);
      arma::uvec dists_sub = dists.subvec(startix,endix);
      int nsp = dists_sub.max() + 1;
      int n1 = samp1.size();
      int n0 = samp0.size();
      // Rcout << "n1 = " << n1 << std::endl;
      // Rcout << "n0 = " << n0 << std::endl;
      arma::umat locx(2,n1,arma::fill::zeros);
      arma::umat locy(2,n0,arma::fill::zeros);
      locx.row(0) = dists_sub.elem(samp1).t();
      locy.row(0) = dists_sub.elem(samp0).t();
      // Rcout << locy << std::endl;
      // Rcout << vals.elem(samp0) << std::endl;
      // Rcout << nsp << std::endl;
      arma::sp_mat x(locx, vals.elem(samp1), nsp, 1);
      arma::sp_mat y(locy, vals.elem(samp0), nsp, 1);
      // Rcout << y << std::endl;
      
      ans.col(ntest) = ShiftScore(x, y, calcP, nresamp);
      ntest++;
      startix = it;
    }
  }
  return(ans);
}