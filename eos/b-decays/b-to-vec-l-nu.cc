/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2019, 2020 Danny van Dyk
 * Copyright (c) 2021 Christoph Bobeth
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/form-factors/form-factors.hh>
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/save.hh>

#include <cmath>
#include <functional>
#include <map>
#include <string>

namespace eos
{
    using std::norm;

    namespace b_to_dstar_l_nu
    {
        struct Amplitudes
        {
            complex<double> a_0;
            complex<double> a_0_T;
            complex<double> a_plus;
            complex<double> a_plus_T;
            complex<double> a_minus;
            complex<double> a_minus_T;
            complex<double> a_P;
            complex<double> a_t;
            complex<double> a_para;
            complex<double> a_para_T;
            complex<double> a_perp;
            complex<double> a_perp_T;
            double mlH;
            double NF;
        };

        // angular observables V's. cf. from [DSD2014], p. 16, redifined V's in order to include NF
        struct AngularObservables
        {
            std::array<double, 12> _vv;

            AngularObservables(const Amplitudes & a)
            {
                // charged lepton velocity in the dilepton rest frame
                const double mlH = a.mlH;
                const double mlH2 = mlH * mlH;
                const double NF = a.NF;

                _vv[0] = NF * 2.0 * (
                            (1.0 + mlH2) * (std::norm(a.a_0) + 16.0 * std::norm(a.a_0_T))
                          + 2.0 * mlH2 * std::norm(a.a_t)
                          + 2.0 * std::norm(a.a_P)
                          + 4.0 * mlH * std::real(a.a_t * std::conj(a.a_P))
                          - 16.0 * mlH * std::real(a.a_0_T * std::conj(a.a_0))
                         );

                _vv[1] = NF * 2.0 * (1.0 - mlH2) * ( - std::norm(a.a_0) + 16.0 * std::norm(a.a_0_T) );

                _vv[2] = - NF * 8.0 * std::real(
                            mlH * (mlH * a.a_t + a.a_P) * std::conj(a.a_0)
                          - 4.0 * (mlH * a.a_t + a.a_P) * std::conj(a.a_0_T)
                        );

                _vv[3] = NF * (
                            (3.0 + mlH2) * (std::norm(a.a_para) + std::norm(a.a_perp)) / 2.0
                          + 8.0 * (1.0 + 3.0 * mlH2) * (std::norm(a.a_para_T) + std::norm(a.a_perp_T))
                          - 16.0 * mlH * std::real(a.a_para_T * std::conj(a.a_para) + a.a_perp_T * std::conj(a.a_perp))
                        );

                _vv[4] = NF * (1.0 - mlH2) * (
                            (std::norm(a.a_para) + std::norm(a.a_perp)) / 2.0
                          - 8.0 * (std::norm(a.a_para_T) + std::norm(a.a_perp_T))
                        );

                _vv[5] = NF * 4.0 * std::real(
                          - a.a_para * std::conj(a.a_perp) 
                          - 16.0 * mlH2 * a.a_para_T * std::conj(a.a_perp_T)
                          + 4.0 * mlH * (a.a_perp_T * std::conj(a.a_para) + a.a_para_T * std::conj(a.a_perp))
                        );

                _vv[6] = NF * (1.0 - mlH2) * (
                          - (std::norm(a.a_para) - std::norm(a.a_perp))
                          + 16.0 * (std::norm(a.a_para_T) - std::norm(a.a_perp_T))
                        );

                _vv[7] = NF * 2.0 * (1.0 - mlH2) * std::imag( a.a_para * std::conj(a.a_perp));

                _vv[8] = NF * std::sqrt(2.0) * (1.0 - mlH2) * std::real(
                            a.a_para * std::conj(a.a_0)
                          - 16.0 * a.a_para_T * std::conj(a.a_0_T)
                        );

                _vv[9] = NF * 2.0 * std::sqrt(2.0) * std::real(
                          - a.a_perp * std::conj(a.a_0)
                          + a.a_para * mlH * std::conj(mlH * a.a_t + a.a_P) 
                          - 16.0 * mlH2 * a.a_perp_T * std::conj(a.a_0_T)
                          + 4.0 * mlH * (a.a_0_T * std::conj(a.a_perp) + a.a_perp_T * std::conj(a.a_0))
                          - 4.0 * a.a_para_T * std::conj(mlH * a.a_t + a.a_P)
                        );

                _vv[10] = NF * 2.0 * std::sqrt(2.0) * std::imag(
                          - a.a_para * std::conj(a.a_0) 
                          + mlH * a.a_perp * std::conj(mlH * a.a_t + a.a_P)
                          + 4.0 * mlH * ( a.a_0_T * std::conj(a.a_para) - a.a_para_T * std::conj(a.a_0))
                          + 4.0 * a.a_perp_T * std::conj(mlH * a.a_t + a.a_P) 
                        );

                _vv[11] = NF * std::sqrt(2.0) * (1.0 - mlH2) * std::imag( a.a_perp * std::conj(a.a_0));
            }

            AngularObservables(const std::array<double, 12> & vv) :
                _vv(vv)
            {
            }

            inline double vv10()    const  { return _vv[0]; }  // J_1c in B->K* ll literature
            inline double vv20()    const  { return _vv[1]; }  // J_2c
            inline double vv30()    const  { return _vv[2]; }  // J_6c
            inline double vv1T()    const  { return _vv[3]; }  // J_1s
            inline double vv2T()    const  { return _vv[4]; }  // J_2s
            inline double vv3T()    const  { return _vv[5]; }  // J_6s
            inline double vv4T()    const  { return _vv[6]; }  // J_3
            inline double vv5T()    const  { return _vv[7]; }  // J_9
            inline double vv10T()   const  { return _vv[8]; }  // J_4
            inline double vv20T()   const  { return _vv[9]; }  // J_5
            inline double vv30T()   const  { return _vv[10];}  // J_7
            inline double vv40T()   const  { return _vv[11];}  // J_8

            // longitudinal polarization amplitude
            inline double normalized_amplitude_polarization_L() const
            {
                return  vv10() - vv20() / 3.0;
            }

            // transverse polarization amplitude
            inline double normalized_amplitude_polarization_T() const
            {
                return  2.0 * (vv1T() - vv2T() / 3.0);
            }

            // redefined decay width
            inline double normalized_decay_width() const
            {
                return 3.0 / 4.0 * ( normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T() );
            }

            // polarization fraction
            inline double f_L() const
            {
                return  normalized_amplitude_polarization_L() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            // polarization fraction from cos(theta_l) distribution; identical to F_L in the SM and the limit m_l -> 0.
            inline double ftilde_L() const
            {
                // (1 - 3 Ftilde_L)  == 16/3 (S2s + S2c/2)
                return 1.0 / 3.0 - 16.0 / 9.0 * (vv2T() + vv20() / 2.0) / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            // a_fb leptonic
            inline double a_fb_leptonic() const
            {
                return  (vv3T() + vv30() / 2.0) / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            // transverse azimuthal asymmetries
            inline double a_c_1() const
            {
                return  4.0 * vv4T() / (3.0 * (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T()));
            }

            inline double a_c_2() const
            {
                return  vv20T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            inline double a_c_3() const
            {
                return  vv10T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            // T-odd CP asymmetries
            inline double a_t_1() const
            {
                return  4.0 * vv5T() / (3.0 * (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T()));
            }

            inline double a_t_2() const
            {
                return  vv30T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }

            inline double a_t_3() const
            {
                return  vv40T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
            }
        };
    }

    /**/
    template <> struct Implementation<BToVectorLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Parameters parameters;

        SwitchOption opt_q;

        UsedParameter hbar;

        UsedParameter tau_B;

        UsedParameter g_fermi;

        SwitchOption opt_l;

        UsedParameter m_l;

        UsedParameter m_B;

        UsedParameter m_V;

        UsedParameter mu;

        bool cp_conjugate;

        inline std::string _process() const
        {
            switch (opt_q.value()[0])
            {
                case 'd':
                case 'u':
                    return std::string("B->D^*");
                    break;

                case 's':
                    return std::string("B_s->D_s^*");
                    break;

                default:
                    throw InternalError("Should never reach this part, either!");
            }

            return "";
        }

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            opt_q(o, "q", { "u", "d", "s" }, "d"),
            hbar(p["hbar"], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            g_fermi(p["G_Fermi"], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            m_B(p["mass::B_" + opt_q.value()], u),
            m_V(p["mass::D_" + opt_q.value() + "^*"], u),
            mu(p["mu"], u),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            form_factors = FormFactorFactory<PToV>::create(_process() + "::" + o.get("form-factors", "BSZ2015"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);

        }

        // normalization cf. [DSD2014] eq. (7), p. 5
        double norm(const double & q2) const
        {
            const double lam  = lambda(m_B * m_B, m_V * m_V, q2);
            const double p    = (lam > 0.0) ? std::sqrt(lam) / (2.0 * m_B) : 0.0;
            
            // normalized prefactor (|Vcb|^2=1)
            return power_of<2>(g_fermi()) * p * q2 * power_of<2>(1.0 - m_l * m_l / q2) / (3.0 * 64.0 * power_of<3>(M_PI) * m_B * m_B);
        }

        b_to_dstar_l_nu::Amplitudes amplitudes(const double & q2) const
        {
            b_to_dstar_l_nu::Amplitudes result;

            // NP contributions in EFT including tensor operator cf. [DSD2014], p. 3
            const WilsonCoefficients<ChargedCurrent> wc = model->wilson_coefficients_b_to_c(opt_l.value(), cp_conjugate);
            const complex<double> gV_pl = wc.cvl() + wc.cvr();  // gV_pl = 1 + gV = 1 + VL + VR = cVL + cVR
            const complex<double> gV_mi = wc.cvl() - wc.cvr();  // gV_mi = 1 - gA = 1 + VL - VR = cVL - cVR
            const complex<double> gP = wc.csr() - wc.csl();
            const complex<double> TL = wc.ct();

            // form factors
            const double aff0  = form_factors->a_0(q2);
            const double aff1  = form_factors->a_1(q2);
            const double aff12 = form_factors->a_12(q2);
            const double vff   = form_factors->v(q2);
            const double tff1  = form_factors->t_1(q2);
            const double tff2  = form_factors->t_2(q2);
            const double tff3  = form_factors->t_3(q2);
            // running quark masses
            const double mbatmu = model->m_b_msbar(mu);
            const double mcatmu = model->m_c_msbar(mu);
            const double lam      = lambda(m_B * m_B, m_V * m_V, q2);
            const double sqrt_lam = (lam > 0.0) ? std::sqrt(lam) : 0.0;
            const double sqrtq2 = std::sqrt(q2);

            // transversity amplitudes A's. cf. [DSD2014], p.17
            result.a_0          = gV_mi * 8.0 * m_B() * m_V() / sqrtq2 * aff12;
            result.a_0_T        = TL / (2.0 * m_V) * ( (m_B * m_B + 3.0 * m_V * m_V - q2) * tff2 - lam * tff3 / (m_B * m_B - m_V * m_V) );
            result.a_plus       = (m_B + m_V) * aff1 * gV_mi - sqrt_lam * vff * gV_pl / (m_B + m_V);
            result.a_minus      = (m_B + m_V) * aff1 * gV_mi + sqrt_lam * vff * gV_pl / (m_B + m_V);
            result.a_plus_T     = TL / sqrtq2 * ( (m_B * m_B - m_V * m_V) * tff2 + sqrt_lam * tff1 );
            result.a_minus_T    = TL / sqrtq2 * ( (m_B * m_B - m_V * m_V) * tff2 - sqrt_lam * tff1 );
            result.a_t          =  sqrt_lam * aff0 * gV_mi / sqrtq2;
            result.a_P          =  sqrt_lam * aff0 * gP / (mbatmu + mcatmu);
            result.a_para       =  (result.a_plus + result.a_minus) / std::sqrt(2.0);
            result.a_para_T     =  (result.a_plus_T + result.a_minus_T) / std::sqrt(2.0);
            result.a_perp       =  (result.a_plus - result.a_minus) / std::sqrt(2.0);
            result.a_perp_T     =  (result.a_plus_T - result.a_minus_T) / std::sqrt(2.0);

            result.mlH = (m_l > 0.0) ? std::sqrt(m_l * m_l / q2) : 0.0;
            result.NF = norm(q2);

            return result;
        }


        std::array<double, 12> _differential_angular_observables(const double & q2) const
        {
            return b_to_dstar_l_nu::AngularObservables(this->amplitudes(q2))._vv;
        }

        // define below integrated observables in generic form
        std::array<double, 12> _integrated_angular_observables(const double & q2_min, const double & q2_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand(std::bind(&Implementation::_differential_angular_observables, this, std::placeholders::_1));
            // second argument of integrate1D is some power of 2
            // Christoph Bobeth: increased integration accuracy for
            //                   difference observables Delta F_L = F_L^mu - F_L^e                      
            return integrate1D(integrand, 4096/*256*/, q2_min, q2_max);
        }

        inline b_to_dstar_l_nu::AngularObservables differential_angular_observables(const double & q2) const
        {
            return b_to_dstar_l_nu::AngularObservables{ _differential_angular_observables(q2) };
        }

        inline b_to_dstar_l_nu::AngularObservables integrated_angular_observables(const double & q2_min, const double & q2_max) const
        {
            return b_to_dstar_l_nu::AngularObservables{ _integrated_angular_observables(q2_min, q2_max) };
        }

        double normalized_decay_width(const double & q2) const
        {
            return this->differential_angular_observables(q2).normalized_decay_width();
        }

        double integrated_pdf_q2(const double & q2_min, const double & q2_max) const
        {
            const double q2_abs_min = power_of<2>(m_l());
            const double q2_abs_max = power_of<2>(m_B() - m_V());

            std::function<double (const double &)> f = std::bind(&Implementation<BToVectorLeptonNeutrino>::normalized_decay_width, this, std::placeholders::_1);
            const double num   = integrate<GSL::QAGS>(f, q2_min,     q2_max);
            const double denom = integrate<GSL::QAGS>(f, q2_abs_min, q2_abs_max);

            return num / denom / (q2_max - q2_min);
        }

        double integrated_pdf_w(const double & w_min, const double & w_max) const
        {
            const double m_B    = this->m_B(), m_B2 = m_B * m_B;
            const double m_V    = this->m_V(), m_V2 = m_V * m_V;
            const double q2_max = m_B2 + m_V2 - 2.0 * m_B * m_V * w_min;
            const double q2_min = m_B2 + m_V2 - 2.0 * m_B * m_V * w_max;

            return integrated_pdf_q2(q2_min, q2_max) * (q2_max - q2_min);// / (w_max - w_min);
        }
    };

    BToVectorLeptonNeutrino::BToVectorLeptonNeutrino(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToVectorLeptonNeutrino>(new Implementation<BToVectorLeptonNeutrino>(p, o, *this))
    {
    }

    BToVectorLeptonNeutrino::~BToVectorLeptonNeutrino()
    {
    }

    /* q^2-differential observables */

    // |Vcb|=1
    double
    BToVectorLeptonNeutrino::normalized_differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).normalized_decay_width() * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::differential_branching_ratio(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).normalized_decay_width() * std::norm(_imp->model->ckm_cb()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::differential_a_fb_leptonic(const double & q2) const
    {
        return _imp->differential_angular_observables(q2).a_fb_leptonic();
    }

    double
    BToVectorLeptonNeutrino::differential_J1c(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv10();
    }

    double
    BToVectorLeptonNeutrino::differential_J1s(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv1T();
    }

    double
    BToVectorLeptonNeutrino::differential_J2c(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv20();
    }

    double
    BToVectorLeptonNeutrino::differential_J2s(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv2T();
    }

    double
    BToVectorLeptonNeutrino::differential_J3(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv4T();
    }

    double
    BToVectorLeptonNeutrino::differential_J4(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * -o.vv10T();
    }

    double
    BToVectorLeptonNeutrino::differential_J5(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv20T();
    }

    double
    BToVectorLeptonNeutrino::differential_J6c(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * -o.vv30();
    }

    double
    BToVectorLeptonNeutrino::differential_J6s(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * -o.vv3T();
    }

    double
    BToVectorLeptonNeutrino::differential_J7(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * -o.vv30T();
    }

    double
    BToVectorLeptonNeutrino::differential_J8(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * o.vv40T();
    }

    double
    BToVectorLeptonNeutrino::differential_J9(const double & q2) const
    {
        auto   o = _imp->differential_angular_observables(q2);
        return 3.0 / 4.0 * -o.vv5T();
    }

    /* q^2-integrated observables */

    // |Vcb|=1
        double
    BToVectorLeptonNeutrino::normalized_integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width() * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::integrated_decay_width(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width() * std::norm(_imp->model->ckm_cb()) / (1.0e-15);
    }

    double
    BToVectorLeptonNeutrino::integrated_branching_ratio(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).normalized_decay_width() * std::norm(_imp->model->ckm_cb()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::integrated_CPave_branching_ratio(const double & q2_min, const double & q2_max) const
    {
       Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return (o.normalized_decay_width() + o_c.normalized_decay_width()) / 2.0 * std::norm(_imp->model->ckm_cb()) * _imp->tau_B / _imp->hbar;
    }

    double
    BToVectorLeptonNeutrino::integrated_pdf_q2(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_pdf_q2(q2_min, q2_max);
    }

    double
    BToVectorLeptonNeutrino::integrated_pdf_w(const double & w_min, const double & w_max) const
    {
        return _imp->integrated_pdf_w(w_min, w_max);
    }

    double
    BToVectorLeptonNeutrino::integrated_a_fb_leptonic(const double & q2_min, const double & q2_max) const
    {
        return _imp->integrated_angular_observables(q2_min, q2_max).a_fb_leptonic();
    }

    double
    BToVectorLeptonNeutrino::integrated_amplitude_polarization_L(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.normalized_amplitude_polarization_L() * std::norm(_imp->model->ckm_cb());
    }

    double
    BToVectorLeptonNeutrino::integrated_amplitude_polarization_T(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.normalized_amplitude_polarization_T() * std::norm(_imp->model->ckm_cb());
    }

    double
    BToVectorLeptonNeutrino::integrated_f_L(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.f_L();
    }

    double
    BToVectorLeptonNeutrino::integrated_ftilde_L(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.ftilde_L();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_c_1(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.a_c_1();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_c_2(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.a_c_2();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_c_3(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.a_c_3();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_t_1(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.a_t_1();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_t_2(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.a_t_2();
    }

    double
    BToVectorLeptonNeutrino::integrated_a_t_3(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return o.a_t_3();
    }

    double
    BToVectorLeptonNeutrino::integrated_J1c(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * o.vv10();
    }

    double
    BToVectorLeptonNeutrino::integrated_J1s(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * o.vv1T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J2c(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * o.vv20();
    }

    double
    BToVectorLeptonNeutrino::integrated_J2s(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * o.vv2T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J3(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * o.vv4T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J4(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * -o.vv10T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J5(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * o.vv20T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J6c(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * -o.vv30();
    }

    double
    BToVectorLeptonNeutrino::integrated_J6s(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * -o.vv3T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J7(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * -o.vv30T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J8(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * o.vv40T();
    }

    double
    BToVectorLeptonNeutrino::integrated_J9(const double & q2_min, const double & q2_max) const
    {
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        return 3.0 / 4.0 * -o.vv5T();
    }

    double
    BToVectorLeptonNeutrino::integrated_CPave_a_fb_leptonic(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv3T() + o_c.vv3T() + o.vv30() / 2.0 + o_c.vv30() / 2.0) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_CPave_f_L(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv10() + o_c.vv10() - o.vv20() / 3.0 - o_c.vv20() / 3.0) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_CPave_ftilde_L(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);
           
        return 1.0 / 3.0 - 3.0 / 4.0 * 16.0 / 9.0 * (o.vv2T() + o_c.vv2T() + o.vv20() / 2.0 + o_c.vv20() / 2.0) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    /*  CP-averaged normalized observables
     *
     *    S_i ~ (J_i + barJ_i) / (Gam_i + barGam_i)
     */
    double
    BToVectorLeptonNeutrino::integrated_S1c(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv10() + o_c.vv10()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S1s(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv1T() + o_c.vv1T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S2c(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv20() + o_c.vv20()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S2s(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv2T() + o_c.vv2T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S3(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv4T() + o_c.vv4T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S4(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv10T() + o_c.vv10T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S5(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv20T() + o_c.vv20T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S6c(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv30() + o_c.vv30()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S6s(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv3T() + o_c.vv3T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S7(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);
        
        return -3.0 / 4.0 * (o.vv30T() + o_c.vv30T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S8(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv40T() + o_c.vv40T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_S9(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv5T() + o_c.vv5T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    /*  CP-asymmetric normalized observables
     *
     *    A_i ~ (J_i - barJ_i) / (Gam_i + barGam_i)
     */
    double
    BToVectorLeptonNeutrino::integrated_A1c(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv10() - o_c.vv10()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A1s(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv1T() - o_c.vv1T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A2c(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv20() - o_c.vv20()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A2s(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv2T() - o_c.vv2T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A3(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv4T() - o_c.vv4T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A4(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv10T() - o_c.vv10T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A5(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv20T() - o_c.vv20T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A6c(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv30() - o_c.vv30()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A6s(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv3T() - o_c.vv3T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A7(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);
        
        return -3.0 / 4.0 * (o.vv30T() - o_c.vv30T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A8(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return 3.0 / 4.0 * (o.vv40T() - o_c.vv40T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    double
    BToVectorLeptonNeutrino::integrated_A9(const double & q2_min, const double & q2_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        
        auto   o = _imp->integrated_angular_observables(q2_min, q2_max);
        _imp->cp_conjugate = true;
        auto   o_c = _imp->integrated_angular_observables(q2_min, q2_max);

        return -3.0 / 4.0 * (o.vv5T() - o_c.vv5T()) / (o.normalized_decay_width() + o_c.normalized_decay_width());
    }

    //* cf. [DSD2014], eq. (6), p. 5 - normalized(|Vcb|=1)
    double
    BToVectorLeptonNeutrino::normalized_four_differential_decay_width(const double & q2, const double & c_theta_l, const double & c_theta_d, const double & phi) const
    {
        // compute d^4 Gamma, cf. [DSD2014], p. 5, eq. (6)
        // Trigonometric identities: Cosine squared of the angles
        double c_theta_d_2 = c_theta_d * c_theta_d;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = std::cos(phi);
        // Sine squared of the angles
        double s_theta_d_2 = 1.0 - c_theta_d_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_d = std::sqrt(s_theta_d_2);
        double s_theta_l = std::sqrt(s_theta_l_2);
        double s_phi = std::sin(phi);
        // Cosine of twice the angle
        //double c_2_theta_d = 2.0 * c_theta_d_2 - 1.0;
        double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
        double c_2_phi = std::cos(2.0 * phi);
        // Sine of twice the angle
        double s_2_theta_d = 2.0 * s_theta_d * c_theta_d;
        double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;
        double s_2_phi = std::sin(2.0 * phi);

        b_to_dstar_l_nu::AngularObservables a_o = _imp->differential_angular_observables(q2);

        double result = 9.0 / 32.0 / M_PI * (
                                            (a_o.vv10() + a_o.vv20() * c_2_theta_l + a_o.vv30() * c_theta_l) * c_theta_d_2
                                            +  (a_o.vv1T() + a_o.vv2T() * c_2_theta_l + a_o.vv3T() * c_theta_l) * s_theta_d_2
                                            +  a_o.vv4T() * s_theta_d_2 * s_theta_l_2 * c_2_phi + a_o.vv10T() * s_2_theta_d * s_2_theta_l * c_phi
                                            +  a_o.vv20T() * s_2_theta_d * s_theta_l * c_phi + a_o.vv5T() * s_theta_d_2 * s_theta_l_2 * s_2_phi
                                            +  a_o.vv30T() * s_2_theta_d * s_theta_l * s_phi +  a_o.vv40T() * s_2_theta_d * s_2_theta_l * s_phi
                                            );

        return result;
    }

    const std::string
    BToVectorLeptonNeutrino::description = "\
    The decay B->D^* l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_q2 = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the charged lepton's helicity angle theta_l in the l-nubar rest frame.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_c_theta_d = "\
    The cosine of the D's helicity angle theta_d in the D-pi rest frame.";

    const std::string
    BToVectorLeptonNeutrino::kinematics_description_phi = "\
    The azimuthal angle between the D-pi plane and the l-nubar  plane.";

    const std::set<ReferenceName>
    BToVectorLeptonNeutrino::references
    {
    };
}
