namespace beam_pp {
  namespace quadrature
  {
    // gauss quadrature for 6th order integration
    constexpr int nquadrature  = 6;
    double weight[nquadrature] = {.171324492e0,.360761573e0,.467913935e0,
				  .467913935e0,.360761573e0,.171324492e0};
    double xloc[nquadrature] = {-.932469514e0,-.661209386e0,-.238619186e0,
				.238619186e0,.661209386e0,.932469514e0};  
  };
}
