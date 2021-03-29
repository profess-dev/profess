// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "kinetic_gga.hpp"
#include "fourier.hpp"

/*
The generalized gradient approximations for $T_s[n]$ may be written in the form
$$
T_{GGA}[n] = \int_\Omega \mathrm{d}\mathbf{r} \, t_{TF} f(s^2) \quad,
$$
where $\Omega$ is the volume, $t_{TF} = c_0 n^{5/3}$ with $c_0 = \tfrac{3}{10}(3\pi^2)^{2/3}$,  the function $f(x)$ defines the enhancement factor, and $s^2 = |\nabla n|^2 / (c_s^2 n^{8/3})$ with $c_s = 2(3\pi^2)^{1/3}$.

The kinetic potential is
$$
\begin{aligned}
\frac{\delta T_{GGA}}{\delta n(\mathbf{r})} =
\frac{5}{3} & c_0 n^{2/3} \left( f(s^2)
    - \frac{8}{5} f'(s^2) s^2 \right) -\nabla \cdot \left(\frac{2 t_{TF} f'(s^2)}{c_s^2 n^{8/3}} \nabla n\right)
\end{aligned}
$$
and the kinetic stress is
$$
\begin{aligned}
\sigma_{\alpha \beta} =
    - \frac{2 \delta_{\alpha \beta}}{3\Omega}  T_{GGA}
    + \frac{2 \delta_{\alpha \beta}}{3\Omega}
            \int_\Omega \mathrm{d}\mathbf{r} \, t_{TF} f'(s^2) s^2 \\
    - \frac{2}{\Omega} \int_\Omega \mathrm{d}\mathbf{r} \,  
            \frac{t_{TF} f'(s^2) }{c_s^2 n^{8/3}} \nabla_\alpha n \cdot \nabla_\beta n
\end{aligned} \quad .
$$
*/

namespace profess
{

KineticGGA::KineticGGA(
        Box box,
        std::function<std::tuple<Double3D,Double3D,Double3D>
                      (Double3D,Double3D)> integrand_and_derivatives,
        std::function<std::tuple<Double3D,Double3D>
                      (Double3D,Double3D)> second_derivatives)
    : Functional::Functional(box),
      _integrand_and_derivatives(integrand_and_derivatives),
      _second_derivatives(second_derivatives)
{
}

std::string KineticGGA::name()
{
    return "KineticGGA";
}

void KineticGGA::set_box(std::array<std::array<double,3>,3> box)
{
    _box = box;
}

double KineticGGA::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double,Double3D> KineticGGA::energy_potential(Double3D den)
{
    // assemble grad_dot_grad
    Double3D grad_x = gradient_x(den, _box);
    Double3D grad_y = gradient_y(den, _box);
    Double3D grad_z = gradient_z(den, _box);
    Double3D grad_dot_grad = grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;

    // compute integrand and derivatives
    Double3D integrand(den.shape());
    Double3D deriv_1(den.shape()), deriv_2(den.shape());
    std::tie(integrand, deriv_1, deriv_2) =
        _integrand_and_derivatives(den, grad_dot_grad);
    Double3D deriv_12(den.shape()), deriv_22(den.shape());
    std::tie(deriv_12, deriv_22) = _second_derivatives(den, grad_dot_grad);

    // return energy and potential
    double ene = integrate(integrand, _box);
    Double3D pot = deriv_1
        - 2.0 * deriv_2 * laplacian(den, _box)
        - 2.0 * deriv_12 * grad_dot_grad
        - 2.0 * deriv_22 * gradient_x(grad_dot_grad, _box) * grad_x
        - 2.0 * deriv_22 * gradient_y(grad_dot_grad, _box) * grad_y
        - 2.0 * deriv_22 * gradient_z(grad_dot_grad, _box) * grad_z;
    return {ene, pot};
}

std::array<std::array<double,3>,3> KineticGGA::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};

    // compute gradient and s2
    auto grad_x = gradient_x(den, _box);
    auto grad_y = gradient_y(den, _box);
    auto grad_z = gradient_z(den, _box);
    auto grad_dot_grad = grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;

    // compute integrand and derivatives
    Double3D integrand(den.shape());
    Double3D deriv_1(den.shape()), deriv_2(den.shape());
    std::tie(integrand, deriv_1, deriv_2) =
        _integrand_and_derivatives(den, grad_dot_grad);

    // compute diagonal-only contributions to stress tensor
    integrand -= (den*deriv_1 + 2.0*grad_dot_grad*deriv_2);
    double s = integrate(integrand, _box) / _box.volume();
    for (int i=0; i<3; ++i) stress[i][i] = s;

    // compute remaining parts of stress tensor
    integrand = grad_x * grad_x * deriv_2;
    stress[0][0] -= 2.0/_box.volume()*integrate(integrand, _box);
    integrand = grad_y * grad_y * deriv_2;
    stress[1][1] -= 2.0/_box.volume()*integrate(integrand, _box);
    integrand = grad_z * grad_z * deriv_2;
    stress[2][2] -= 2.0/_box.volume()*integrate(integrand, _box);
    integrand = grad_y * grad_z * deriv_2;
    stress[1][2] = -2.0/_box.volume()*integrate(integrand, _box);
    stress[2][1] = stress[1][2];
    integrand = grad_x * grad_z * deriv_2;
    stress[0][2] = -2.0/_box.volume()*integrate(integrand, _box);
    stress[2][0] = stress[0][2];
    integrand = grad_x * grad_y * deriv_2;
    stress[0][1] = -2.0/_box.volume()*integrate(integrand, _box);
    stress[1][0] = stress[0][1];

    return stress;
}

}
