#include <cmath>
#include <cfloat>
#include <fenv.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <mpfr.h>

using std::cerr;
using std::cout;
using std::endl;
using std::function;
using std::get;
using std::isinf;
using std::isnan;
using std::make_pair;
using std::make_tuple;
using std::ostream;
using std::pair;
using std::setprecision;
using std::setw;
using std::string;
using std::tuple;
using std::vector;

template<typename F, typename MF>
pair<double, double> compute(double x, F std_f, MF mpfr_f, tuple<string, int, mpfr_rnd_t> rnd) {
    mpfr_t mpfr_tmp;
    mpfr_init2(mpfr_tmp, 256);
    double std_result, mpfr_result;
    int fe_rnd = get<1>(rnd);
    mpfr_rnd_t mpfr_rnd = get<2>(rnd);
    int fe_old_rnd = fegetround();
    fesetround(fe_rnd);
    std_result = std_f(x);
    fesetround(fe_old_rnd);
    mpfr_set_d(mpfr_tmp, x, mpfr_rnd);
    mpfr_f(mpfr_tmp, mpfr_tmp, mpfr_rnd);
    mpfr_result = mpfr_get_d(mpfr_tmp, mpfr_rnd);
    return make_pair(std_result, mpfr_result);
}

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

ostream & display(ostream & out, string f, double x, double r1, double r2, string rnd) {
    out << f << "(" << setprecision(15) << x << ") = " << setprecision(15) << r1 << "\t" << rnd << endl;
    out << setw(25) << " = " << setprecision(15) << r2 << endl;
    return out;
}

int main() {
    srand(time(NULL));
    // Rounding Modes
    vector<tuple<string, int, mpfr_rnd_t>> rnds = {make_tuple("NEAREST",  FE_TONEAREST,  MPFR_RNDN),
                                                   make_tuple("UPWARD",   FE_UPWARD,     MPFR_RNDU),
                                                   make_tuple("DOWNWARD", FE_DOWNWARD,   MPFR_RNDD),
                                                   make_tuple("TO_ZERO",  FE_TOWARDZERO, MPFR_RNDZ)};
    vector<tuple<string, function<double(double)>, function<int(mpfr_t, mpfr_t, mpfr_rnd_t)>,
                 double, double, double, double> > funcs =
    //       FuncName    C_FUN  MPFR_FUN   DOM_L  DOM_U  RANGE_L    RANGE_U
        {make_tuple("sin",   sin,   mpfr_sin,  -31.4, +31.4, -1.00,     +1.00),
         make_tuple("cos",   cos,   mpfr_cos,  -31.4, +31.4, -1.00,     +1.00),
         make_tuple("tan",   tan,   mpfr_tan,  -31.4, +31.4, -INFINITY, +INFINITY),
         make_tuple("acos",  acos,  mpfr_acos, -1.00, +1.00, -INFINITY, +INFINITY),
         make_tuple("asin",  asin,  mpfr_asin, -1.00, +1.00, -INFINITY, +INFINITY),
         make_tuple("atan",  atan,  mpfr_atan, -1.00, +1.00, -INFINITY, +INFINITY),
         make_tuple("cosh",  cosh,  mpfr_cosh, -31.4, +31.4, 1.0,       +INFINITY),
         make_tuple("sinh",  sinh,  mpfr_sinh, -31.4, +31.4, -INFINITY, +INFINITY),
         make_tuple("tanh",  tanh,  mpfr_tanh, -31.4, +31.4, -1.0,      +1.0),
         make_tuple("exp",   exp,   mpfr_exp,  -31.4, +31.4, 0.0,       +INFINITY),
         make_tuple("log",   log,   mpfr_log,   0.01, +31.4, -INFINITY, +INFINITY),
         make_tuple("log10", log10, mpfr_log10, 0.01, +31.4, -INFINITY, +INFINITY),
         make_tuple("sqrt",  sqrt,  mpfr_sqrt,  0.00, +31.4, 0,         +INFINITY),
        };
    double std_result, mpfr_result, eps, x;
    unsigned const max_iter = 100000;
    pair<double, double> result;
    cout << setw(15) << "Function"     << setw(15) << "Round" << setw(15) << "Imprecise"
         << setw(15) << "Unreasonable" << setw(15) << "Infty" << setw(15) << "Nan"
         << setw(15) << "Total"        << endl;
    for (auto const & func : funcs) {
        for (auto const & rnd : rnds) {
            unsigned int imprecise = 0, nan = 0, inf = 0, unreasonable = 0;
            for (unsigned i = 0; i < max_iter; i++) {
                x = fRand(get<3>(func), get<4>(func));
                result = compute(x, get<1>(func), get<2>(func), rnd);
                std_result  = result.first;
                mpfr_result = result.second;

                eps = 1e15 * fabs(mpfr_result - std::nextafter(mpfr_result, DBL_MAX));
                if (fabs(mpfr_result - std_result) > eps) {
                    imprecise++;
                    display(cerr, get<0>(func), x, std_result, mpfr_result, get<0>(rnd));
                    if (isnan(std_result)) { nan++; }
                    else if (isinf(std_result)) { inf++; }
                    else if (std_result < get<5>(func) || get<6>(func) < std_result) { unreasonable++; }
                }
            }
            unsigned int total = imprecise + unreasonable + inf + nan;
            if (total > 0) {
                cout << setw(15) << get<0>(func) << setw(15) << get<0>(rnd)
                     << setw(15) << imprecise    << setw(15) << unreasonable
                     << setw(15) << inf          << setw(15) << nan
                     << setw(15) << total        << endl;
            }
        }
    }
    return 0;
}
