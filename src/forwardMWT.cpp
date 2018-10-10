#include "forwardMWT.hpp"

#include <functional>

/*
template <typename T>
std::vector<T> forwardMWT(const int lev, const int deg, const T Lmin, const T
Lmax, std::function<T(T,PARAMS<T>)> foo, PARAMS<T> &params) {

    matrix2d<T> FMWT = operator_two_scale<T>(deg, static_cast<int>(pow(2, lev))
);

    // Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
    // [-1,+1] for performing quadrature.

    int quad_num = 10;

    std::vector<T> quad_x, quad_w;

    lgwt(quad_x,quad_w, quad_num, -1, 1);

    int N = quad_x.size();

    // Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
    // to order k.

    auto p_val1D = dlegendre2(quad_x,deg)[0];
    matrix2d<T> p_val(deg,N);
    for (auto i=0; i<N; ++i) {
            p_val(0,i) = p_val1D[i];
            p_val(1,i) = p_val1D[N+i];
    }

    // Get grid spacing.

    int n = std::pow(2,lev);
    T h = (Lmax-Lmin) / n;
    int dof_1D = deg*n;
    matrix2d<T> f(dof_1D,1);

    // Initial Condition for f

    std::vector<T> x(N,0);
    for (auto i=0; i<=n-1; ++i) {

        // Map quad_x from [-1,+1] to [Lmin,Lmax] physical domain.

        for (auto j=0; j<N; ++j) {
            x[j] = h * (quad_x[j]/2.0+1.0/2.0+i) + Lmin;
        }

        // Get the f(v) initial condition at the quadrature points.

        std::vector<T> fHere;
        matrix2d<T>this_quad(N,1);
        for (auto j=0; j<N; ++j) {
            fHere.push_back( foo(x[j],params) );
            this_quad(j,0) = (quad_w[j] * fHere[j]);
        }

        // Generate the coefficients for DG basis

        auto tmp = p_val * this_quad;

        f.set(i*deg,i*deg+deg-1,0,0,tmp);

    }

    f = f * (h * std::sqrt(1.0/h)/2.0);


    // Transfer to multi-DG bases

    auto tmp2 = FMWT * f;

    return( tmp2.data_container );

}
*/
