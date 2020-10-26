/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Christoph Bobeth
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

#include <test/test.hh>
#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/form-factors/parametric-bgl1997-impl.hh>

#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BGL1997FormFactorsTest :
    public TestCase
{
    public:
        BGL1997FormFactorsTest() :
            TestCase("BGL1997_form_factor_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            //p["QCD::alpha_s(MZ)"] = 0.1176;

            /* Outer function */
            {
                BGL1997FormFactors<BToDstar> ff(p, Options{ });
                const double m_B(BToDstar::mB);
                const double m_V(BToDstar::mV);
                const double t_0((m_B + m_V)* power_of<2>(std::sqrt(m_B) - std::sqrt(m_V)));

                //TEST_CHECK_NEARLY_EQUAL( 0.0535063, ff._z(1.0, 2.0), eps);

                //_phi(s, t_0, K, a, b, c, chi);
                TEST_CHECK_NEARLY_EQUAL( 0.0331832,  ff._phi(-2.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0324458,  ff._phi(+1.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0316657,  ff._phi(+4.0, t_0, 48, 3, 3, 2, 3.1e-03), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.488275,   ff._phi(-2.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.470779,   ff._phi(+1.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.452784,   ff._phi(+4.0, t_0, 48, 3, 3, 1, 3.1e-03), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.00817026, ff._phi(-2.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00822179, ff._phi(+1.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.00827232, ff._phi(+4.0, t_0, 16, 1, 1, 1, 3.1e-03), eps);
            }

            /* B -> D^* FFs */
            {
                BGL1997FormFactors<BToDstar> ff(p, Options{ });

                p["B->D^*::a^g_0@BGL1997"] = 0.1e-02;
                p["B->D^*::a^g_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^g_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^g_3@BGL1997"] = 0.4e-02;
                
                p["B->D^*::a^f_0@BGL1997"] = 0.1e-02;
                p["B->D^*::a^f_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^f_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^f_3@BGL1997"] = 0.4e-02;
                
                p["B->D^*::a^F1_0@BGL1997"] = 0.1e-02;
                p["B->D^*::a^F1_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F1_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F1_3@BGL1997"] = 0.4e-02;
                
                p["B->D^*::a^F2_0@BGL1997"] = 0.1e-02;
                p["B->D^*::a^F2_1@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F2_2@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F2_3@BGL1997"] = 0.4e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.0120945, ff.g(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0131032, ff.g(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0143205, ff.g(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.598215, ff.f(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.620714, ff.f(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.647914, ff.f(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 12.448531, ff.F1(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 12.736923, ff.F1(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 13.101888, ff.F1(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.0463648, ff.F2(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0498538, ff.F2(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0540411, ff.F2(+4.0), eps);

                p["B->D^*::a^g_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^g_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^g_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^g_3@BGL1997"] = 0.1e-02;
                
                p["B->D^*::a^f_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^f_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^f_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^f_3@BGL1997"] = 0.1e-02;
                
                p["B->D^*::a^F1_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^F1_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F1_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F1_3@BGL1997"] = 0.1e-02;
                
                p["B->D^*::a^F2_0@BGL1997"] = 0.4e-02;
                p["B->D^*::a^F2_1@BGL1997"] = 0.3e-02;
                p["B->D^*::a^F2_2@BGL1997"] = 0.2e-02;
                p["B->D^*::a^F2_3@BGL1997"] = 0.1e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.0461237, ff.g(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0508857, ff.g(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0566738, ff.g(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 2.281365, ff.f(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.410521, ff.f(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 2.564134, ff.f(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 47.474000, ff.F1(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 49.463380, ff.F1(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 51.850996, ff.F1(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.176818, ff.F2(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.193605, ff.F2(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.213869, ff.F2(+4.0), eps);
             }

            /* B -> D FFs*/
            {
                BGL1997FormFactors<BToD> ff(p, Options{ });

                p["B->D::a^f+_0@BGL1997"] = 0.1e-02;
                p["B->D::a^f+_1@BGL1997"] = 0.2e-02;
                p["B->D::a^f+_2@BGL1997"] = 0.3e-02;
                p["B->D::a^f+_3@BGL1997"] = 0.4e-02;
                
                p["B->D::a^f0_0@BGL1997"] = 0.1e-02;
                p["B->D::a^f0_1@BGL1997"] = 0.2e-02;
                p["B->D::a^f0_2@BGL1997"] = 0.3e-02;
                p["B->D::a^f0_3@BGL1997"] = 0.4e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.0862919, ff.f_p(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0911714, ff.f_p(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.0970192, ff.f_p(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 0.439279, ff.f_0(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.435522, ff.f_0(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.432496, ff.f_0(+4.0), eps);

                p["B->D::a^f+_0@BGL1997"] = 0.4e-02;
                p["B->D::a^f+_1@BGL1997"] = 0.3e-02;
                p["B->D::a^f+_2@BGL1997"] = 0.2e-02;
                p["B->D::a^f+_3@BGL1997"] = 0.1e-02;
                
                p["B->D::a^f0_0@BGL1997"] = 0.4e-02;
                p["B->D::a^f0_1@BGL1997"] = 0.3e-02;
                p["B->D::a^f0_2@BGL1997"] = 0.2e-02;
                p["B->D::a^f0_3@BGL1997"] = 0.1e-02;

                TEST_CHECK_NEARLY_EQUAL( 0.327129, ff.f_p(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.352242, ff.f_p(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 0.382318, ff.f_p(+4.0), eps);

                TEST_CHECK_NEARLY_EQUAL( 1.66529, ff.f_0(-2.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.68264, ff.f_0(+1.0), eps);
                TEST_CHECK_NEARLY_EQUAL( 1.70431, ff.f_0(+4.0), eps);
            }

        }
} BGL1997_form_factor_test;
