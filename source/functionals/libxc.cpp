// Copyright (c) 2019-2020 William C. Witt
// SPDX-License-Identifier: MIT

#include "libxc.hpp"

#include "fourier.hpp"

/*
For functionals of the form
$$  
F[n] = \int_\Omega \mathrm{d}\mathbf{r} \, f\left(n,|\nabla n|^2\right),
$$  
where $\Omega$ is the volume, the first functional derivative is  
$$  
\frac{\delta F}{\delta n(\mathbf{r})} =  
\frac{\partial f}{\partial n}
    - 2 \nabla \cdot \frac{\partial f}{\partial |\nabla n|^2} \nabla n
$$  
and the stress tensor is
$$
\begin{aligned}
\sigma_{\alpha \beta} =
    & \delta_{\alpha\beta} \frac{1}{\Omega} F[n] \\
    & -\delta_{\alpha\beta} \frac{1}{\Omega} 
            \int_{\Omega} \mathrm{d}\mathbf{r} \, 
                n \frac{\partial f}{\partial n} \\
    & -\frac{2}{\Omega} \int_\Omega \mathrm{d}\mathbf{r} 
         \left[ \delta_{\alpha\beta} |\nabla n|^2
               + \nabla_\alpha n \nabla_\beta n \right] 
               \frac{\partial f}{\partial |\nabla n|^2}
\end{aligned}.
$$
*/

namespace profess
{

LibXC::LibXC(Box box, std::vector<int> xc_func_ids)
    : Functional::Functional(box)
{
    _init(xc_func_ids);
}

LibXC::~LibXC()
{
    for (size_t i=0; i<_xc_func.size(); ++i)
        xc_func_end(&_xc_func[i]);
}

std::string LibXC::name()
{
    return "libxc";
}

void LibXC::set_box(std::array<std::array<double,3>,3> box)
{
    _box = Box(box);
}

double LibXC::energy(Double3D den)
{
    return std::get<0>(energy_potential(den));
}

std::tuple<double, Double3D> LibXC::energy_potential(Double3D den)
{
    Double3D exc(den), vxc(den);

    if (_xc_all_lda) {
        // use libxc to obtain the...
        exc.fill(0.0);  // 1. energy density (per electron)
        vxc.fill(0.0);  // 2. derivative with respect to density
        for (size_t i=0; i<_xc_func.size(); ++i) {
            Double3D texc(den), tvxc(den);
            xc_lda_exc_vxc(&_xc_func[i], den.size(), den.data(),
                    texc.data(), tvxc.data());
            exc += texc;
            vxc += tvxc;
        }    
    } else {
        // compute grad(density) and grad_dot_grad
        Double3D grad_x = gradient_x(den, _box);
        Double3D grad_y = gradient_y(den, _box);
        Double3D grad_z = gradient_z(den, _box);
        Double3D gdg = grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;
        Double3D vgdg(den.shape());
        // use libxc to obtain the...
        exc.fill(0.0);  // 1. energy density (per electron)
        vxc.fill(0.0);  // 2. derivative with respect to density
        vgdg.fill(0.0); // 3. derivative with respect to grad_dot_grad
        for (size_t i=0; i<_xc_func.size(); ++i) {
            if (_xc_func[i].info->family == XC_FAMILY_LDA) {
                Double3D texc(den), tvxc(den);
                xc_lda_exc_vxc(&_xc_func[i], den.size(), den.data(),
                        texc.data(), tvxc.data());
                exc += texc;
                vxc += tvxc;
            } else {
                Double3D texc(den), tvxc(den), tvgdg(den);
                xc_gga_exc_vxc(&_xc_func[i], den.size(), den.data(), gdg.data(),
                        texc.data(), tvxc.data(), tvgdg.data());
                exc += texc;
                vxc += tvxc;
                vgdg += tvgdg;
            }
        }    
        // add final contribution to potential
        grad_x *= vgdg;
        grad_y *= vgdg;
        grad_z *= vgdg;
        grad_x = gradient_x(grad_x, _box);
        grad_y = gradient_y(grad_y, _box);
        grad_z = gradient_z(grad_z, _box);
        vxc = vxc - 2.0*(grad_x + grad_y + grad_z);
    }
    exc *= den;
    return {integrate(exc, _box), std::move(vxc)};
}

std::array<std::array<double,3>,3> LibXC::stress(Double3D den)
{
    std::array<std::array<double,3>,3> stress =
            {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    if (_xc_all_lda) {
        double s;
        Double3D pot(den);
        std::tie(s, pot) = energy_potential(den);
        pot *= den;
        s = (s-integrate(pot,_box)) / _box.volume();
        for (int i=0; i<3; ++i) stress[i][i] = s;
    } else {
        // compute grad(density) and grad_dot_grad
        Double3D grad_x = gradient_x(den, _box);
        Double3D grad_y = gradient_y(den, _box);
        Double3D grad_z = gradient_z(den, _box);
        Double3D gdg = grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;
        // use libxc to obtain the...
        Double3D exc(den), vxc(den), vgdg(den);
        exc.fill(0.0);  // 1. energy density (per electron)
        vxc.fill(0.0);  // 2. derivative with respect to density
        vgdg.fill(0.0); // 3. derivative with respect to grad_dot_grad
        for (size_t i=0; i<_xc_func.size(); ++i) {
            if (_xc_func[i].info->family == XC_FAMILY_LDA) {
                Double3D texc(den), tvxc(den);
                xc_lda_exc_vxc(&_xc_func[i], den.size(), den.data(),
                        texc.data(), tvxc.data());
                exc += texc;
                vxc += tvxc;
            } else {
                Double3D texc(den), tvxc(den), tvgdg(den);
                xc_gga_exc_vxc(&_xc_func[i], den.size(), den.data(), gdg.data(),
                        texc.data(), tvxc.data(), tvgdg.data());
                exc += texc;
                vxc += tvxc;
                vgdg += tvgdg;
            }
        }    
        // compute diagonal-only contributions to stress tensor
        exc = exc*den - den*vxc - 2.0*gdg*vgdg;
        double s = integrate(exc,_box) / _box.volume();
        for (int i=0; i<3; ++i) stress[i][i] = s;
        // compute remaining parts of stress tensor
        exc = grad_x * grad_x * vgdg;
        stress[0][0] -= 2.0/_box.volume()*integrate(exc,_box);
        exc = grad_y * grad_y * vgdg;
        stress[1][1] -= 2.0/_box.volume()*integrate(exc,_box);
        exc = grad_z * grad_z * vgdg;
        stress[2][2] -= 2.0/_box.volume()*integrate(exc,_box);
        exc = grad_y * grad_z * vgdg;
        stress[1][2] = -2.0/_box.volume()*integrate(exc,_box);
        stress[2][1] = stress[1][2];
        exc = grad_x * grad_z * vgdg;
        stress[0][2] = -2.0/_box.volume()*integrate(exc,_box);
        stress[2][0] = stress[0][2];
        exc = grad_x * grad_y * vgdg;
        stress[0][1] = -2.0/_box.volume()*integrate(exc,_box);
        stress[1][0] = stress[0][1];
    }
    return stress;
}

void LibXC::_init(std::vector<int> xc_func_ids)
{
    _xc_func = std::vector<xc_func_type>(xc_func_ids.size());
    _xc_all_lda = true;
    for (size_t i=0; i<_xc_func.size(); ++i) {
        xc_func_init(&_xc_func[i], xc_func_ids[i], XC_UNPOLARIZED);
        if (_xc_func[i].info->family == XC_FAMILY_GGA)
            _xc_all_lda = false;
    }
}

}
